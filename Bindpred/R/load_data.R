#' Analyzes and processes the repertoire sequencing data from cellranger vdj. This provides information on the single-cell level for each clone.
#' @param VDJ.directory Character vector with each element containing the path to the output of cellranger vdj runs. To integrate the repertoire with transcriptome data GEX data can by provided inside the corresponding folder. The name of the GEX folder must contain "gex" in his name (case insensitive). This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder.
#' @param filtered.prod Logical indicating if the contigs file should be filtered by productive cells. Default set to TRUE
#' @param filtered.confidence Logical indicating if the contigs file should be filtered by contigs with high_confidence. Default set to TRUE
#' @param filtered.is.cell Logical indicating if the contigs file should be filtered by sequences that only come from cells. Default set to TRUE
#' @param red.method Method of dimensionality reduction used. Options are "pca", "tsne", "umap".
#' @param n.dims Number of dimension to be used in the chosen reduction. Default is 20.
#' @return Returns a list of dataframes containing information on the cell or clonotype level.
#' @export
#' @examples
#' \dontrun{
#' check_load_data <- load_data(VDJ.directory = "path/to/cellranger/outs/", clonotype.level = FALSE)
#' }
#'

load_data <- function(VDJ.directory, clonotype.level, filter.prod, filter.confidence, filter.is.cell, red.method, n.dims) {

  require(tidyverse, quietly = T)
  require(Rsamtools)
  require(RcppSimdJson, quietly = T)
  require(Seurat, quietly = T)
  require(seqinr, quietly = T)
  require(Biostrings)

  if (!all(dir.exists(VDJ.directory))) print("Some directory supplied are not supported")
    VDJ.directory[dir.exists(VDJ.directory)]
  if (missing(clonotype.level)) clonotype.level <- FALSE
  if (missing(filter.prod)) filter.prod <- TRUE
  if (missing(filter.confidence)) filter.confidence <- TRUE
  if (missing(filter.is.cell)) filter.is.cell <- TRUE
  if (missing(red.method)) red.method <- "pca"
  if (missing(n.dims)) n.dims <- 1:10

  # Reading input files
  if (all(str_detect(VDJ.directory, "mouse"))) {sample_number <- c(substr(sub(".*mouse", "", VDJ.directory), 1, 1))}
  else {sample_number <- seq(1, length(VDJ.directory))}

  print("Reading input files")
  contig.list <- lapply(paste0(VDJ.directory, "/filtered_contig_annotations.csv"), read.csv, stringsAsFactors=F, sep=",", header=T)
  fasta.list <- lapply(paste0(VDJ.directory, "/filtered_contig.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  reference.list <- lapply(paste0(VDJ.directory, "/concat_ref.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  annotations.list <- lapply(paste0(VDJ.directory, "/all_contig_annotations.json"), fload)
  clonotype.list <- lapply(paste0(VDJ.directory, "/clonotypes.csv"), read.table, stringsAsFactors = F, sep=",", header=T)

  # Filter clones that contain exactly one heavy/beta and one light/alpha chain
  print("Filtering clones that have both HC and LC")
  clonotype.list <- map(clonotype.list, ~filter(., str_count(cdr3s_aa, "IGH:") == 1 & str_count(cdr3s_aa, "IGL:|IGK:") == 1))
  clonotype.list <- map2(clonotype.list, sample_number, ~mutate(., sample = .y))

  clonotype.list <- map(clonotype.list,
                    ~mutate(.x, CDR3H_aa = gsub(".*IGH:(.+);.*", "\\1", cdr3s_aa), CDR3L_aa = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_aa),
                                CDR3H_nt = gsub(".*IGH:(.+);.*", "\\1", cdr3s_nt), CDR3L_nt = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_nt)))

  # Filtering contigs by specified parameters
  contig.list <- map(contig.list, ~filter(., case_when(filter.is.cell ~ is_cell %>% str_detect("(?i)true")),
                                             case_when(filter.confidence ~ high_confidence %>% str_detect("(?i)true")),
                                             case_when(filter.prod ~productive %>% str_detect("(?i)true"))) %>%
                                             rename("raw_clonotype_id" = "clonotype_id") %>%
                                             dplyr::select(-c(is_cell, full_length, productive, reads, high_confidence)))

  # Filtering annotations by specified parameters
  annotations.list <- map(annotations.list, ~filter(., case_when(filter.is.cell ~ is_cell == T),
                                                       case_when(filter.confidence ~ high_confidence == T),
                                                       case_when(filter.prod ~ productive == T)) %>%
                                                       dplyr::select(-c(is_cell, productive, high_confidence)))

  # Adding gene expression information # -------------------------------------------------------------------------
  print("Processing gene expression information")

  GEX.list <- map_if(VDJ.directory,
    ~any(dir.exists(unique(str_extract(list.files(., recursive = T, "features.tsv.gz|barcodes.tsv.gz|matrix.mtx.gz", full.names = T), ".*/")))),
       ~Read10X(unique(str_extract(list.files(., recursive = T, "features.tsv.gz|barcodes.tsv.gz|matrix.mtx.gz", full.names = T), ".*/"))) %>%
          CreateSeuratObject)

  antibody_gene_names <- c("IGHA", "IGHG", "IGHM", "IGHD", "IGHE", "IGHJ", "IGK", "IGHV", "JCHAIN", "IGL", "TRAV", "TRAC", "TRBC",
                           "TRGC", "TRDC", "TRBD", "TRBJ", "TRGV", "TRGJ", "TRGJ", "TRDV", "TRDD", "TRDJ", "TRBV", "TRAJ")

  GEX.list <- lapply(GEX.list, function(x) x[!toupper(row.names(x)) %in% antibody_gene_names])
  GEX.list <- lapply(GEX.list, function(x) x$barcode <- {as.character(gsub(colnames(x), pattern = "-1", replacement = "")); return(x)})
  GEX.list <- lapply(GEX.list, function(x) {if (typeof(x) == "S4") x$percent.mt <- PercentageFeatureSet(x, pattern="^mt-"); return(x)})
  GEX.list <- lapply(GEX.list, function(x) {if (typeof(x) == "S4")
    RunPCA(ScaleData(FindVariableFeatures(NormalizeData(subset(x, percent.mt<20 & nFeature_RNA>100 & nFeature_RNA<2500), verbose=F), verbose=F), verbose=F), npcs = 20, verbose=F)})

  GEX.list <- lapply(GEX.list, function(x)
    if (!is.null(x)) switch(red.method, tsne = {RunTSNE(x, dims = 1:n.dims, verbose=F)},
                                        umap = {RunUMAP(x, dims = 1:n.dims, verbose=F)},
                                        pca = {x}))

  red_data <- lapply(GEX.list, function(x) {if (typeof(x) == "S4")
    rownames_to_column(as.data.frame(x@reductions[[red.method]]@cell.embeddings), var = "barcode")})

  red_data <- map2(red_data, contig.list, ~{if (!is.null(.x)) inner_join(.x, dplyr::select(.y, c(barcode, clonotype_id)), by = "barcode")})

  # -------------------------------------------------------------------------

  # Adding cigar string for position of mutations and counting number of mutations
  print("Processing somatic hypermutation data")
  shm <- map_if(VDJ.directory, ~file.exists(paste0(., "/concat_ref.bam")), ~as.data.frame(scanBam(paste0(., "/concat_ref.bam"))) %>%
    rename("qname" = "contig_id") %>%
      mutate(mutation_count = str_extract_all(cigar, "[0-9]+(?=[X])")) %>%
      dplyr::select(contig_id, mutation_count, cigar) %>% mutate(mutation_count = sapply(.$mutation_count, function(x) sum(as.numeric(x)))))


  if (clonotype.level == TRUE) {

    # Counting the most common gene for the V,D,J segments for each clonotype and chain
    print("Extracting gene usage information")
    genes <- map(contig.list, ~group_by(., clonotype_id, chain, barcode) %>%
               summarize_at(vars(ends_with("gene")), ~names(which.max(table(.)))) %>%
                inner_join(filter(., chain == "IGH"), filter(., str_detect(chain, "IGK|IGL")), by = "clonotype_id", suffix = c("_HC", "_LC")) %>%
                  distinct(clonotype_id, .keep_all = T) %>% dplyr::select(-c(barcode_LC, chain_HC, chain_LC, barcode_HC)))

    clonotype.list <- clonotype.list <- map2(clonotype.list, genes, ~left_join(.x, .y, by = "clonotype_id", suffix = c("", "")))

    # Adding mean of pca components
    clonotype.list <- map2(clonotype.list, red_data, ~{if (!is.null(.y))
      left_join(.x, group_by(.y, clonotype_id) %>% dplyr::select(-barcode) %>% summarise_each(funs(mean(., na.rm = T))), by = "clonotype_id")
      else .x} %>% rename("CDR3H_nt" = "cdr3_nt_HC", "CDR3L_nt" = "cdr3_nt_LC"))

    return(clonotype.list)
  }


  if (clonotype.level == FALSE) {

    # Fusing HC and LC contigs from same cell in one row
    vdj.per.cell <- map(contig.list,
      ~left_join(filter(., chain %in% "IGH"),
        filter(dplyr::select(., barcode, cdr3, cdr3_nt, chain, contig_id, umis, ends_with("gene")), chain %in% c("IGK", "IGL")), by = "barcode", suffix = c("_HC", "_LC")))

    # Pasted CDR3
    vdj.per.cell <- map(vdj.per.cell, ~unite(., cdr3_HC, "chain_HC", "cdr3_HC", sep = ":") %>%
                                       unite(., cdr3_LC, "chain_LC", "cdr3_LC", sep = ":") %>%
                                       unite(., cdr3s_aa, "cdr3_HC", "cdr3_LC", sep = ";", remove = F))

    # Adding SHM
    vdj.per.cell <- map2(vdj.per.cell, shm, ~{if (is.data.frame(.y)) left_join(.x, .y, by = c("contig_id_HC" = "contig_id"), suffix = c("_HC", "_LC")) else .x})
    vdj.per.cell <- map2(vdj.per.cell, shm, ~{if (is.data.frame(.y)) left_join(.x, .y, by = c("contig_id_LC" = "contig_id"), suffix = c("_HC", "_LC")) else .x})
    vdj.per.cell <- mapply(cbind, vdj.per.cell, "sample" = sample_number, SIMPLIFY = F)

    # Adding aa sequence
    vdj.per.cell <- map2(vdj.per.cell, annotations.list, ~left_join(.x, dplyr::select(.y, aa_sequence, contig_name), by = c("contig_id_HC" = "contig_name"), suffix=c("_HC", "_LC")))
    vdj.per.cell <- map2(vdj.per.cell, annotations.list, ~left_join(.x, dplyr::select(.y, aa_sequence, contig_name), by = c("contig_id_LC" = "contig_name"), suffix=c("_HC", "_LC")))

    # Adding nt sequence
    trimmed <- map(annotations.list, ~filter(., !is.na(start_codon_pos)) %>%
        mutate(., sequence = str_sub(sequence, start_codon_pos + 1)) %>% dplyr::select(sequence, contig_name))

    vdj.per.cell <- map2(vdj.per.cell, trimmed, ~left_join(.x, .y, by = c("contig_id_HC" = "contig_name"), suffix = c("_HC", "_LC")))
    vdj.per.cell <- map2(vdj.per.cell, trimmed, ~left_join(.x, .y, by = c("contig_id_LC" = "contig_name"), suffix = c("_HC", "_LC")))

    # Adding pca components
    vdj.per.cell <- map2(vdj.per.cell, red_data, ~{if (!is.null(.y)) left_join(.x, dplyr::select(.y, -clonotype_id), by = "barcode") else .x})

    # Removing contig_id
    vdj.per.cell <- map(vdj.per.cell, ~dplyr::select(., -c(contig_id_HC, contig_id_LC, raw_consensus_id)))

    return(vdj.per.cell)
  }
}

