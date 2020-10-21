
VDJ.directory <- "Immunizations/mouse1"
GEX.directory <- "GEX"
label_data <- "2020_09_labels/2019_09_label_data_immunizations.csv"
clone.strategy <- "CDRH3.homology"

# Function that loads data and produces a feature matrix --------------------------------------------------

load_data <- function(VDJ.directory, GEX.directory, clonotype.level, clone.strategy) {

  require(Platypus)
  require(dplyr)
  require(jsonlite)
  require(stringr)
  require(tibble)
  require(Rsamtools)
  if (missing(VDJ.directory)) {stop("Missing VDJ sequence information. Please provide directory containing VDJ data")}

  mouse_number <- grep("mouse.*", VDJ.directory)
  vdj <- VDJ_analyze(VDJ.directory, filter.1HC.1LC = TRUE)

  vdj[[1]]$mouse <- mouse_number
  all_annotations <- read_json(paste0(VDJ.directory, "/all_contig_annotations.json"))
  filtered_annotations <- read_json(paste0(VDJ.directory, "/consensus_annotations.json"))
  concat_fasta <- read.fasta(paste0(VDJ.directory, "/concat_ref.fasta"), as.string = TRUE, forceDNAtolower = FALSE)
  ref <- as.data.frame(scanBam("Immunizations/mouse1/concat_ref.bam")[[1]])
  shm <- data.frame(cbind(lapply(str_extract_all(ref$cigar, "[0-9]+(?=[X,D])"), function(x) {sum(as.numeric(x))}), ref$qname), check.names = T)
  colnames(shm) <- c("SHM", "contig_id")

  if (clonotype.level == TRUE) {
    # retrieving VDJ sequence trimmed
    seq <- c()
    clonotype <- c()
    for (i in 1:length(filtered_annotations)) {
      if (filtered_annotations[[i]]$productive) {
        first_region_index <- which(sapply(filtered_annotations[[i]]$annotations, function(x) {x$feature$region_type == "L-REGION+V-REGION"}))
        last_region_index <- which(sapply(filtered_annotations[[i]]$annotations, function(x) {x$feature$region_type == "J-REGION"}))

        vdj_start_index <- filtered_annotations[[i]]$annotations[[first_region_index]]$contig_match_start
        vdj_end_index <- filtered_annotations[[i]]$annotations[[last_region_index]]$contig_match_end
        seq[i] <- substr(filtered_annotations[[i]]$sequence, vdj_start_index, vdj_end_index)
        clonotype[i] <- filtered_annotations[[i]]$clonotype
      }
    }
    VDJ_sequence <- data.frame(cbind(sequences = seq, clonotype_id = clonotype))
    HC_subset <- sapply(filtered_annotations, function(x) {x$annotations[[1]]$feature$chain == "IGH"})
    LC_subset <- sapply(filtered_annotations, function(x) {x$annotations[[1]]$feature$chain == "IGK"})

    VDJ_sequence <- inner_join(VDJ_sequence[HC_subset, ], VDJ_sequence[LC_subset, ], by = "clonotype_id")
    colnames(VDJ_sequence) <- c("VDJ_HC_sequence", "clonotype_id", "VDJ_LC_sequence")

    # Adding VDJ trimmed sequence
    left_join(vdj[[1]], VDJ_sequence, by = "clonotype_id")

    if (!missing(GEX.directory)) {
      GEX.object <- automate_GEX(GEX.directory, integration.method = "scale.data")
      pca_data <- rownames_to_column(as.data.frame(GEX.object[[1]]@reductions[["pca"]]@cell.embeddings), var = "barcode")
      all.cells <- bind_rows(VDJ_per_clone(vdj, VDJ.directory))
      temp_pca_data <- left_join(pca_data, bind_rows(all.cells), by = "barcode", keep = FALSE)
      means <- temp_pca_data %>%
        group_by(clonotype_id)  %>%
          summarize_at(vars(starts_with("PC_")), list(~ mean(., na.rm = TRUE)))

      # Adding mean of pca components
      vdj[[1]] <- left_join(vdj[[1]], means, by = "clonotype_id")
    }
    return(vdj)
  }

  if (clonotype.level == FALSE) {
    # retrieving VDJ sequence trimmed
    seq <- c()
    barcode <- c()
    clonotype <- c()
    for (i in 1:length(all_annotations)) {
      if (all_annotations[[i]]$productive) {
        first_region_index <- which(sapply(all_annotations[[i]]$annotations, function(x) {x$feature$region_type == "L-REGION+V-REGION"}))
        last_region_index <- which(sapply(all_annotations[[i]]$annotations, function(x) {x$feature$region_type == "J-REGION"}))
        vdj_start_index <- all_annotations[[i]]$annotations[[first_region_index]]$contig_match_start
        vdj_end_index <- all_annotations[[i]]$annotations[[last_region_index]]$contig_match_end
        seq[i] <- substr(all_annotations[[i]]$sequence, vdj_start_index, vdj_end_index)
        barcode[i] <- all_annotations[[i]]$barcode
      }
    }
    VDJ_sequence <- data.frame(cbind(sequences = seq, barcode = barcode))
    HC_subset <- sapply(all_annotations, function(x) {x$annotations[[1]]$feature$chain == "IGH"})
    LC_subset <- sapply(all_annotations, function(x) {x$annotations[[1]]$feature$chain == "IGK"})

    VDJ_sequence <- inner_join(VDJ_sequence[HC_subset, ], VDJ_sequence[LC_subset, ], by = "barcode")
    colnames(VDJ_sequence) <- c("VDJ_HC_sequence", "barcode", "VDJ_LC_sequence")

    vdj.per.clone <- VDJ_per_clone(vdj, VDJ.directory)[[1]]
    for (i in 1:length(vdj.per.clone)) {
      vdj.per.clone[[i]]$cdr3s_aa <- vdj[[1]]$cdr3s_aa[[i]][[1]]
      vdj.per.clone[[i]]$mouse <- mouse_number

      # Adding VDJ trimmed sequence
      vdj.per.clone[[i]] <- left_join(vdj.per.clone[[i]], VDJ_sequence, by = "barcode")

      # Adding SHM quantification
      left_join(vdj.per.clone[[i]], shm, by = c("contig_id_lc", "contig_id"))
      left_join(vdj.per.clone[[i]], shm, by = c("contig_id_hc", "contig_id"))
      }
    if (!missing(GEX.directory)) {
      GEX.object <- automate_GEX(GEX.directory, integration.method = "scale.data")[[1]]
      pca_data <- rownames_to_column(as.data.frame(GEX.object@reductions[["pca"]]@cell.embeddings), var = "barcode")
      for (i in 1:length(vdj.per.clone)) {
        # Adding mean of pca components
        vdj.per.clone[[i]] <- left_join(vdj.per.clone[[i]], pca_data, by = "barcode")
      }
    }
    return(vdj.per.clone)
  }
}

features_CELL <- load_data(VDJ.directory, clonotype.level = FALSE)
features_CLONE <- load_data(VDJ.directory, clonotype.level = TRUE)



# function that labels cell or clonotype data -------------------------------------------------------------

label_clonotype_data <- function(features, labels, share.specific, share.affinity) {

  if (missing(share.specific)) {share.specific <- TRUE}
  if (missing(share.affinity)) {share.affinity <-  FALSE}

  mouse_number <- features[[1]]$mouse[1]
  labels <- read.csv(labels) %>% filter(mouse == mouse_number)
  depth <- function(x) ifelse(is.list(x), 1L + max(sapply(x, depth)), 0L)

  ### Cell level
  if (depth(features) == 2) {
    for (i in 1:length(features)) {
      if (share.specific == TRUE) {
        features[[i]] <- left_join(features[[i]], select(labels, ELISA_bind, clonotype_id), by = "clonotype_id")
      }
      else {
        features[[i]] <- left_join(features[[i]], select(labels, ELISA_bind, cdr3s_aa), by = "cdr3s_aa")
      }
      if (share.affinity == TRUE) {
        features[[i]] <- left_join(features[[i]], select(labels, octet, clonotype_id), by = "clonotype_id")
      }
      else {
        features[[i]] <- left_join(features[[i]], select(labels, octet, cdr3s_aa), by = "cdr3s_aa")
      }
    }
    return(features)
  }
  ### Clonotype level
  if (depth(features) == 1) {
    features <- left_join(features, select(labels, octet, ELISA_bind, clonotype_id), by = "clonotype_id")
    return(features)
  }
}

labeled_CELL <- label_clonotype_data(features_CELL, labels = label_data, share.specific = FALSE, share.affinity = FALSE)
labeled_CLONE <- label_clonotype_data(features_CLONE, labels = label_data, share.specific = FALSE, share.affinity = TRUE)





consensus <- read_json("Immunizations/mouse1/all_contig_annotations.json")

# retrieving VDJ sequence trimmed
seq <- c()
barcode <- c()
clonotype <- c()
for (i in 1:length(consensus)) {
  first_region_index <- which(sapply(consensus[[i]]$annotations, function(x) {x$feature$region_type == "L-REGION+V-REGION"}))
  last_region_index <- which(sapply(consensus[[i]]$annotations, function(x) {x$feature$region_type == "J-REGION"}))
  if (any(first_region_index) & any(last_region_index)) {
    vdj_start_index <- consensus[[i]]$annotations[[first_region_index]]$contig_match_start
    vdj_end_index <- consensus[[i]]$annotations[[last_region_index]]$contig_match_end
    seq[i] <- substr(consensus[[i]]$sequence, vdj_start_index, vdj_end_index)
    barcode[i] <- consensus[[i]]$barcode
  }
}
VDJ_sequence <- data.frame(cbind(sequences = seq, barcode = barcode))
HC_subset <- sapply(consensus, function(x) {x$annotations[[1]]$feature$chain == "IGH"})
LC_subset <- !HC_subset
VDJ_sequence <- inner_join(VDJ_sequence[HC_subset, ], VDJ_sequence[LC_subset, ], by = "barcode")
colnames(VDJ_sequence) <- c("VDJ_HC_sequence", "barcode", "VDJ_LC_sequence")

vdj.per.clone[[1]] <- left_join(vdj.per.clone[[1]], VDJ_sequence, by = "barcode")

try <- VDJ_sequence %>% count(sequences)
unique(VDJ_sequence$sequences)

for (i in 1:length(vdj.per.clone)) {
  for (j in 1:length(vdj.per.clone[[i]]$VDJ_HC_sequence)) {
    temp_start <- regexpr(vdj.per.clone[[i]]$VDJ_HC_sequence[j], vdj.per.clone[[i]]$full_HC_sequence[j])[1]
    temp_end <- temp_start + str_length(vdj.per.clone[[i]]$VDJ_HC_sequence[j])
    print(str_length(vdj.per.clone[[i]]$full_HC_sequence) == str_length(vdj.per.clone[[i]]$full_HC_germline))
    # try <- substring(vdj.per.clone[[i]]$VDJ_HC_sequence[j], temp_start, temp_end)
    # print(try)
  }
}

trimmed_germline <- lapply(dj.per.clone[[i]]$full_HC_germline, substr, temp_start) substr(vdj.per.clone[[i]]$full_HC_germline, )
try <- vdj.per.clone[[1]] %>% mutate(trimmed_germline = substr(full_HC_germline, ))

library(Rsamtools)
ref <- as.data.frame(scanBam("Immunizations/mouse1/concat_ref.bam")[[1]])
shm <- data.frame(cbind(lapply(str_extract_all(ref$cigar, "[0-9]+(?=[X,D])"), function(x) {sum(as.numeric(x))}), ref$qname))


for (i in 1:length(features_CELL)) {
  for (j in 1:length(features_CELL[[i]])) {
    cigar_hc[i] <- ref$cigar[ref$qname == features_CELL[[i]]$contig_id_hc[j]]
    cigar_lc[i] <- ref$cigar[ref$qname == features_CELL[[i]]$contig_id_lc[j]]
    str_extract_all(ref$cigar[i], "[0-9]+(?=[X,D])")
    temp_annotation_hc <- consensus[which(sapply(consensus, function(x) {x$info$raw_consensus_id == features_CELL[[i]]$contig_id_hc[j]}))]
    temp_annotation_lc <- consensus[which(sapply(consensus, function(x) {x$info$raw_consensus_id == features_CELL[[i]]$contig_id_lc[j]}))]
  }
}
    first_region_index <- which(sapply(temp_annotation_hc[[i]]$annotations, function(x) {x$feature$region_type == "L-REGION+V-REGION"}))
    last_region_index <- which(sapply(consensus[[i]]$annotations, function(x) {x$feature$region_type == "J-REGION"}))
    if (any(first_region_index) & any(last_region_index)) {
      vdj_start_index <- consensus[[i]]$annotations[[first_region_index]]$contig_match_start
      vdj_end_index <- consensus[[i]]$annotations[[last_region_index]]$contig_match_end
      seq[i] <- substr(consensus[[i]]$sequence, vdj_start_index, vdj_end_index)
      barcode[i] <- consensus[[i]]$barcode
  }
}
  ref$cigar[ref$rname == "clonotype1_concat_ref_1"]
  if (str_detect(ref$cigar[i], "X")) {
    mismatches[i] <- str_count(ref$cigar[i], "X(.*)")
    name[i] <- ref[i]$qname
  }
}

lapply(str_extract_all(ref$cigar[ref$rname == "clonotype1_concat_ref_1"][1], "[0-9]{1,4}"), function(x) {sum(as.numeric(x))})

mismatch_dataframe <- as.data.frame(mismatches, name)




