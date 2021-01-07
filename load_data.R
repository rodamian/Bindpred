
directories <- list.dirs("Immunizations", recursive = F)[2:3]
labels <- "2020_09_labels/2019_09_label_data_immunizations.csv"

# Function that loads data and produces a feature matrix --------------------------------------------------

load_data <- function(VDJ.directory, clonotype.level, clone.strategy, filter.unique) {

  require(tidyverse)
  require(Rsamtools)
  require(RcppSimdJson)
  require(parallel)
  require(Seurat)
  require(seqinr)
  require(progress)
  require(Biostrings)

  # Progress bar
  pb <- progress_bar$new(total = 100)

  if(.Platform$OS.type == "unix") options(mc.cores = detectCores()-1) else options(mc.cores = 1)
  if (missing(clonotype.level)) {clonotype.level <- FALSE}

  # Reading input files -----------------------------------------------------

  if (all(str_detect(VDJ.directory, "mouse"))) {mouse_number <- c(substr(sub(".*mouse", "", VDJ.directory), 1, 1))}
  else {mouse_number <- seq(1, length(VDJ.directory))}

  print("Reading in input files")
  contig.list <- lapply(paste0(VDJ.directory, "/filtered_contig_annotations.csv"), read.csv, stringsAsFactors=F, sep=",", header=T)
  fasta.list <- lapply(paste0(VDJ.directory, "/filtered_contig.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  reference.list <- lapply(paste0(VDJ.directory, "/concat_ref.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  annotations.list <- lapply(paste0(VDJ.directory, "/all_contig_annotations.json"), fload)
  clonotype.list <- lapply(paste0(VDJ.directory, "/clonotypes.csv"), read.table, stringsAsFactors = F, sep=",", header=T)

  # Filter out clones that do not contain exactly one heavy/beta and one light/alpha chain
  clonotype.list <- lapply(clonotype.list, filter, str_count(cdr3s_aa, "IGH:|IGL:|IGK:|TRB:|TRA:") == 2)

  clonotype.list <- map(clonotype.list,
    ~mutate(.x, CDR3H_aa = gsub(".*IGH:(.+);.*", "\\1", cdr3s_aa), CDR3L_aa = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_aa),
              CDR3H_nt = gsub(".*IGH:(.+);.*", "\\1", cdr3s_nt), CDR3L_nt = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_nt)))

  # Filtering contigs by different measures (add as option)
  contig.list <- map(contig.list, ~filter(., is_cell == "True", high_confidence == "True", productive == "True") %>% droplevels %>%
    rename("raw_clonotype_id" = "clonotype_id") %>% group_by(clonotype_id) %>%
      mutate(barcodes = gsub(toString(unique(barcode)), pattern = ", ", replacement = ";")))

  # Counting the most common gene for the V,D,J segments for each clonotype and chain
  genes <- map(contig.list, ~group_by(., clonotype_id, chain, barcodes) %>%
            summarize_at(vars(ends_with("gene")), ~names(which.max(table(.)))) %>%
              left_join(filter(., chain == "IGH"), filter(., str_detect(chain, "IGK|IGL")), by = "clonotype_id", suffix = c("_HC", "_LC")))

  clonotype.list <- clonotype.list <- map2(clonotype.list, genes, ~left_join(.x, .y, by = "clonotype_id"))

  # Adding cigar string for position of mutations and counting number of mutations
  shm <- map_if(VDJ.directory, ~file.exists(paste0(., "/concat_ref.bam")), ~as.data.frame(scanBam(paste0(., "/concat_ref.bam"))) %>%
          rename("qname" = "contig_id") %>% mutate(., clonotype_id = str_extract(rname, "[^_]+")) %>%
            mutate(mutation_count = str_count(cigar, "[0-9]+(?=[X,D])")))

  if (clonotype.level == TRUE) {
    if (!missing(clone.strategy)) {clonotype.list <- VDJ_clonotype(clonotype.list, clone.strategy)}

    clonotype.list <- map2(clonotype.list, mouse_number, ~mutate(., mouse = .y))
    GEX.object <- map_if(VDJ.directory, ~dir.exists(paste0(., "/GEX")), ~Read10X(paste0(., "/GEX")) %>% CreateSeuratObject)


      pca_data <- rownames_to_column(as.data.frame(GEX.object[[1]]@reductions[["pca"]]@cell.embeddings), var = "barcode")
      all.cells <- bind_rows(VDJ_per_cell(vdj, VDJ.directory)[[1]])
      temp_pca_data <- left_join(pca_data, bind_rows(all.cells), by = "barcode", keep = FALSE)
      means <- temp_pca_data %>% group_by(raw_clonotype_id)  %>%
          summarize_at(vars(starts_with("PC_")), list(~ mean(., na.rm = TRUE)))

      # Adding mean of pca components
      vdj[[1]] <- left_join(vdj[[1]], means, by = c("clonotype_id" = "raw_clonotype_id"))
    }
    return(vdj[[1]])
  }

  if (clonotype.level == FALSE) {

    vdj.per.cell <- list(data.frame())
    for (i in 1:length(clonotype.list)) {

    # Reference and barcodes ---------------------------------------------------

    references_HC <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {
      temp_index <- gsub("_consensus_", "_concat_ref_",
                         contig.list[[i]]$raw_consensus_id[contig.list[[i]]$raw_clonotype_id ==
                     x & contig.list[[i]]$chain == "IGH" & contig.list[[i]]$is_cell == "True"])
      references_HC[[x]] <- as.character(reference.list[[i]][unique(temp_index[temp_index != "None"])])
    }

    references_LC <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {
      temp_index <- gsub("_consensus_", "_concat_ref_",
                         contig.list[[i]]$raw_consensus_id[contig.list[[i]]$raw_clonotype_id ==
                     x & contig.list[[i]]$chain == "IGK" & contig.list[[i]]$is_cell == "True"])
      references_LC[[x]] <- as.character(reference.list[[i]][unique(temp_index[temp_index != "None"])])
    }

    barcode <- list()
    for (x in clonotype.list[[i]]$clonotype_id) {barcode[[x]] <- filter(contig.list[[i]], clonotype_id == x) %>% pull(barcode)}

    # Trimmed sequences -------------------------------------------------------

    # Selecting only sequences with start and end annotations
    annotations.list[[i]] <- annotations.list[[i]] %>% filter(sapply(annotations.list[[i]]$annotations,
      function(x) sum(x$feature$region_type %in% c("L-REGION+V-REGION", "J-REGION")) == 2), high_confidence, productive)

    # Extracts trimmed sequence for each cell in the clonotype
    per_clone <- mclapply(barcode, function(cells) seq <- filter(annotations.list[[i]], barcode %in% cells))

    seq_trimmed <- mclapply(per_clone, function(x)
      data.frame(sequence = substr(x$sequence,
                 map(x$annotations, ~pull(filter(.x, feature$region_type == "L-REGION+V-REGION"), contig_match_start)+1),
                 map(x$annotations, ~pull(filter(.x, feature$region_type == "J-REGION"), contig_match_end)-1)),
                 contig_id = x$contig_name)
    )

    # Pasted cdr3 -------------------------------------------------------------

    pasted_cdr3 <- mclapply(barcode, function(x) {

      cdr3_HC <- filter(contig.list[[i]], barcode %in% x, productive=="True", chain=="IGH") %>% select (cdr3, barcode, chain)
      cdr3_LC <- filter(contig.list[[i]], barcode %in% x, productive=="True", chain %in% c("IGK", "IGL")) %>% select (cdr3, barcode, chain)

      cdr3 <- full_join(cdr3_HC, cdr3_LC, by="barcode", suffix = c("_HC", "_LC"))

      data.frame(pasted_cdr3 = paste0(cdr3$chain_HC, ":", cdr3$cdr3_HC, ";", cdr3$chain_LC, ":",cdr3$cdr3_LC),
                 CDR3H = cdr3$cdr3_HC, CDR3L = cdr3$cdr3_LC, barcode = cdr3$barcode)
    })

    # Extracting gene usage features from contig.list -------------------------

    extract_feature <- function(barcodes, trimmed_sequence, feature, chain) {
      matching <- contig.list[[i]]$barcode %in% barcodes & contig.list[[i]]$chain == chain
      data <- data.frame(barcode = contig.list[[i]]$barcode[matching], contig.list[[i]][matching, feature],
                         full_seq = as.character(fasta.list[[i]][contig.list[[i]]$contig_id[matching]]))
      data <- inner_join(data, trimmed_sequence, by = "contig_id")
    }

    data_HC <- mcMap(extract_feature, barcode, seq_trimmed, MoreArgs = list(c("umis", "contig_id", "c_gene", "v_gene", "j_gene", "clonotype_id"), "IGH"))
    data_LC <- mcMap(extract_feature, barcode, seq_trimmed, MoreArgs = list(c("umis", "contig_id", "c_gene", "v_gene", "j_gene"), "IGK"))

    vdj.per.cell[[i]] <- mcMap(inner_join, data_HC, data_LC, by = "barcode", suffix = list(c("_HC", "_LC")))
    vdj.per.cell[[i]] <- mcMap(left_join, vdj.per.cell[[i]], pasted_cdr3, by = "barcode", suffix = list(c("", "")))

    # Creating new column with trimmed sequence in amino acid sequence
    vdj.per.cell[[i]] <- map(vdj.per.cell[[i]], ~mutate(.x, sequence_HC_aa = as.character(Biostrings::translate(DNAStringSet(sequence_HC)))))
    vdj.per.cell[[i]] <- map(vdj.per.cell[[i]], ~mutate(.x, sequence_LC_aa = as.character(Biostrings::translate(DNAStringSet(sequence_LC)))))

    # Adding SHM quantification
    if (file.exists(paste0(VDJ.directory, "/concat_ref.bam"))) {
      vdj.per.cell[[i]] <- mclapply(vdj.per.cell[[i]], left_join, shm, by=c("contig_id_HC"="contig_id"), suffix=c("", "_HC"))
      vdj.per.cell[[i]] <- mclapply(vdj.per.cell[[i]], left_join, shm, by=c("contig_id_LC"="contig_id"), suffix=c("", "_LC"))
    }

    # Adding GEX data
    # if (dir.exists(paste0(VDJ.directory,"/GEX"))) {
    #   GEX.object <- automate_GEX(paste0(VDJ.directory,"/GEX"), integration.method = "scale.data")[[1]]
    #   pca_data <- rownames_to_column(as.data.frame(GEX.object@reductions[["pca"]]@cell.embeddings), var = "barcode")
    #   vdj.per.cell[[i]] <- mclapply(vdj.per.cell[[i]], left_join, pca_data, by = "barcode")
    # }

    vdj.per.cell[[i]] <- bind_rows(vdj.per.cell[[i]]) # %>% rename("clonotype_id"="raw_clonotype_id")
    vdj.per.cell[[i]] <- mapply(cbind, vdj.per.cell[[i]], "mouse_number" = mouse_number, SIMPLIFY = F)

    return(vdj.per.cell)
    }
  }
}


# function that labels cell or clonotype data -------------------------------------------------------------

label_data <- function(features, labels, share.specific, share.affinity) {

  if (missing(share.specific)) {share.specific <- TRUE}
  if (missing(share.affinity)) {share.affinity <-  FALSE}

  mouse_number <- unique(features$mouse)

  if (is.null(mouse_number)) mouse_number <- unique(features$mouse_number)

  labels <- read.csv(labels) %>% filter(mouse == mouse_number)

  if (share.specific == TRUE) {
    features <- left_join(features, select(labels, ELISA_bind, clonotype_id), by = "clonotype_id")
  }
  else {
    features <- left_join(features, select(labels, ELISA_bind, cdr3s_aa), by = c("cdr3s_aa"))
  }
  if (share.affinity == TRUE) {
    features <- left_join(features, select(labels, octet, clonotype_id), by = "clonotype_id")
  }
  else {
    features <- left_join(features, select(labels, octet, cdr3s_aa), by = c("cdr3s_aa"))
  }
  return(features)
}


# Writing features to csv -------------------------------------------------

features <- load_data(directories, clone.strategy = "Hvj.Lvj.CDR3length.CDRH3homology")
features_labeled <- lapply(features, label_data, labels)


if (clonotype_level) {level <- "clono"} else {level <-  "cell"}

lapply(features_labeled, function(x)
  write.csv(x, paste0("features/", level, "_features_mouse", unique(x$mouse)), row.names = F))

