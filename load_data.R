#' Analyzes and processes the repertoire sequencing data from cellranger vdj. This provides information on the single-cell level for each clone.
#' @param VDJ.directory Character vector with each element containing the path to the output of cellranger vdj runs. To integrate the repertoire with transcriptome data GEX data can by provided inside the corresponding folder. This pipeline assumes that the output file names have not been changed from the default 10x settings in the /outs/ folder.
#' @param filtered.prod Logical indicating if the contigs file should be filtered by productive cells. Default set to TRUE
#' @param filtered.confidence Logical indicating if the contigs file should be filtered by contigs with high_confidence. Default set to TRUE
#' @param filtered.is.cell Logical indicating if the contigs file should be filtered by . Default set to TRUE
#' @return Returns a list of dataframes containing
#' @export
#' @examples
#' \dontrun{
#' check_load_data <- load_data(VDJ.directory = "path/to/cellranger/outs/", clonotype.level = FALSE, )
#' }
#'

directories <- list.dirs("Immunizations", recursive = F)
labels <- "2020_09_labels/2019_09_label_data_immunizations.csv"

# Function that loads data and produces a feature matrix --------------------------------------------------

load_data <- function(VDJ.directory, clonotype.level, clone.strategy, filter.prod, filter.confidence, filter.is.cell) {

  require(tidyverse)
  require(Rsamtools)
  # require(RcppSimdJson)
  require(Seurat)
  require(seqinr)
  require(Biostrings)

  if (missing(clonotype.level)) clonotype.level <- FALSE
  if (missing(filter.prod)) filter.prod <- TRUE
  if (missing(filter.confidence)) filter.confidence <- TRUE
  if (missing(filter.is.cell)) filter.is.cell <- TRUE

  # Reading input files
  if (all(str_detect(VDJ.directory, "mouse"))) {mouse_number <- c(substr(sub(".*mouse", "", VDJ.directory), 1, 1))}
  else {mouse_number <- seq(1, length(VDJ.directory))}

  print("Reading input files")
  contig.list <- lapply(paste0(VDJ.directory, "/filtered_contig_annotations.csv"), read.csv, stringsAsFactors=F, sep=",", header=T)
  fasta.list <- lapply(paste0(VDJ.directory, "/filtered_contig.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  reference.list <- lapply(paste0(VDJ.directory, "/concat_ref.fasta"), read.fasta, as.string=T, seqonly=F, forceDNAtolower=F)
  # annotations.list <- lapply(paste0(VDJ.directory, "/all_contig_annotations.json"), fload)
  clonotype.list <- lapply(paste0(VDJ.directory, "/clonotypes.csv"), read.table, stringsAsFactors = F, sep=",", header=T)

  # Filter out clones that do not contain exactly one heavy/beta and one light/alpha chain
  print("Filtering clones that have both HC and LC")
  clonotype.list <- map(clonotype.list, filter, str_count(cdr3s_aa, "IGH:|IGL:|IGK:|TRB:|TRA:") == 2)
  clonotype.list <- map2(clonotype.list, mouse_number, ~mutate(., mouse = .y))

  clonotype.list <- map(clonotype.list,
                        ~mutate(.x, CDR3H_aa = gsub(".*IGH:(.+);.*", "\\1", cdr3s_aa), CDR3L_aa = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_aa),
                                CDR3H_nt = gsub(".*IGH:(.+);.*", "\\1", cdr3s_nt), CDR3L_nt = gsub(".*IG[K-L]:(.+)", "\\1", cdr3s_nt)))

  # Filtering contigs by different measures
  contig.list <- map(contig.list, ~filter(., case_when(filter.is.cell ~ is_cell == "True"),
                                             case_when(filter.confidence ~ high_confidence == "True"),
                                             case_when(filter.prod ~productive == "True")) %>%
                                             rename("raw_clonotype_id" = "clonotype_id") %>%
                                             select(-c(is_cell, full_length, productive, reads, high_confidence)))

  # Adding gene expression information
  print("Processing gene expression information")
  GEX.list <- map_if(VDJ.directory, ~dir.exists(paste0(., "/GEX")), ~Read10X(paste0(., "/GEX")) %>% CreateSeuratObject)

  antibody_gene_names <- c("IGHA", "IGHG", "IGHM", "IGHD", "IGHE", "IGHJ", "IGK", "IGHV", "JCHAIN", "IGL", "TRAV", "TRAC", "TRBC",
                           "TRGC", "TRDC", "TRBD", "TRBJ", "TRGV", "TRGJ", "TRGJ", "TRDV", "TRDD", "TRDJ", "TRBV", "TRAJ")

  GEX.list <- lapply(GEX.list, function(x) x[!toupper(row.names(x)) %in% antibody_gene_names])
  GEX.list <- lapply(GEX.list, function(x) x$barcode <- {as.character(gsub(colnames(x), pattern = "-1", replacement = "")); return(x)})
  GEX.list <- lapply(GEX.list, function(x) {if (typeof(x) == "S4") x$percent.mt <- PercentageFeatureSet(x, pattern="^mt-"); return(x)})
  GEX.list <- lapply(GEX.list, function(x) {if (typeof(x) == "S4")
    RunPCA(ScaleData(FindVariableFeatures(NormalizeData(subset(x, percent.mt < 20 & nFeature_RNA > 100 & nFeature_RNA < 2500)))), npcs = 20, ndims.print = 1:4, nfeatures.print = 10)})

  pca_data <- lapply(GEX.list, function(x) {if (typeof(x) == "S4")
    rownames_to_column(as.data.frame(x@reductions$pca@cell.embeddings), var = "barcode")})

  pca_data <- map2(pca_data, contig.list, ~{if (!is.null(.x)) inner_join(.x, select(.y, c(barcode, clonotype_id)), by = "barcode")})


  if (clonotype.level == TRUE) {

    contig.list <- map(contig.list, ~group_by(., clonotype_id) %>%
                         mutate(barcodes = gsub(toString(unique(barcode)), pattern = ", ", replacement = ";")))

    if (!missing(clone.strategy)) {clonotype.list <- VDJ_clonotype(clonotype.list, clone.strategy)}

    # Counting the most common gene for the V,D,J segments for each clonotype and chain
    print("Extracting gene usage information")
    genes <- map(contig.list, ~group_by(., clonotype_id, chain, barcodes) %>%
               summarize_at(vars(ends_with("gene")), ~names(which.max(table(.)))) %>%
                inner_join(filter(., chain == "IGH"), filter(., str_detect(chain, "IGK|IGL")), by = "clonotype_id", suffix = c("_HC", "_LC")) %>%
                  distinct(clonotype_id, .keep_all = T) %>% select(-c(barcodes_LC, chain_HC, chain_LC)) %>% rename("barcodes_HC" = "barcodes"))

    clonotype.list <- clonotype.list <- map2(clonotype.list, genes, ~left_join(.x, .y, by = "clonotype_id", suffix = c("", "")))

    # Adding mean of pca components
    clonotype.list <- map2(clonotype.list, pca_data, ~{if (!is.null(.y))
      left_join(.x, group_by(.y, clonotype_id) %>% select(-barcode) %>% summarise_each(funs(mean(., na.rm = T))), by = "clonotype_id")
      else .x})

    return(clonotype.list)
  }


  if (clonotype.level == FALSE) {
    # Fusing HC and LC contigs from same cell in one row
    vdj.per.cell <- map(contig.list,
      ~left_join(filter(., chain %in% "IGH"),
        filter(select(., barcode, cdr3, cdr3_nt, chain, contig_id, umis, ends_with("gene")), chain %in% c("IGK", "IGL")), by = "barcode", suffix = c("_HC", "_LC")))

    # Pasted CDR3
    vdj.per.cell <- map(vdj.per.cell, ~unite(., cdr3_HC, "chain_HC", "cdr3_HC", sep = ":") %>%
                                       unite(., cdr3_LC, "chain_LC", "cdr3_LC", sep = ":") %>%
                                       unite(., cdr3s_aa, "cdr3_HC", "cdr3_LC", sep = ";", remove = F))

    # Adding cigar string for position of mutations and counting number of mutations
    print("Processing somatic hypermutation data")
    shm <- map_if(VDJ.directory, ~file.exists(paste0(., "/concat_ref.bam")), ~as.data.frame(scanBam(paste0(., "/concat_ref.bam"))) %>%
                    rename("qname" = "contig_id") %>% mutate(mutation_count = str_count(cigar, "[0-9]+(?=[X,D])")) %>%
                    select(contig_id, mutation_count, cigar))

    vdj.per.cell <- map2(vdj.per.cell, shm, ~if (is.data.frame(.y)) left_join(.x, .y, by = c("contig_id_HC" = "contig_id"), suffix=c("", "_HC")) else .x)
    vdj.per.cell <- map2(vdj.per.cell, shm, ~if (is.data.frame(.y)) left_join(.x, .y, by = c("contig_id_LC" = "contig_id"), suffix=c("", "_LC")) else .x)
    vdj.per.cell <- mapply(cbind, vdj.per.cell, "mouse_number" = mouse_number, SIMPLIFY = F)

    # # References
    # print("Processing reference sequence and barcodes")
    # reference.list <- map(reference.list, ~plyr::ldply(.) %>% rename(".id" = "clonotype_id", "V1" = "reference"))
    # vdj.per.cell <- map2(vdj.per.cell, reference.list,
    #                      ~left_join(.x, filter(.y, str_detect(clonotype_id, "concat_ref_1")), by = "clonotype_id", suffix = c("", "")) %>%
    #                        left_join(.x, filter(.y, str_detect(clonotype_id, "concat_ref_2")), by = "clonotype_id", suffix = c("", "")))

    return(vdj.per.cell)
  }

}

# function that labels cell or clonotype data -------------------------------------------------------------

label_data <- function(features, labels, share.specific, share.affinity) {

  if (missing(share.specific)) {share.specific <- TRUE}
  if (missing(share.affinity)) {share.affinity <-  FALSE}

  mouse_number <- map(features, ~unique(.$mouse))
  labels <- read.csv(labels) %>% filter(mouse %in% mouse_number)

  if (share.specific == TRUE) {
    features <- map(features, ~left_join(., select(labels, ELISA_bind, clonotype_id), by = "clonotype_id"))
  }
  else {
    features <- map(features, ~left_join(., select(labels, ELISA_bind, cdr3s_aa), by = "cdr3s_aa"))
  }
  if (share.affinity == TRUE) {
    features <- map(features, ~left_join(., select(labels, octet, clonotype_id), by = "clonotype_id"))
  }
  else {
    features <- map(features, ~left_join(., select(labels, octet, cdr3s_aa), by = "cdr3s_aa"))
  }
  return(features)
}

features <- label_data(load_data(directories, clonotype.level = F), labels)
