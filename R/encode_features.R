#' Encodes repertoire features with different methods. Gene usage features are encoded using one hot encoding irrespective of the encoding method used. Sequence features can be encoded using multiple methods.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer", "protr", "protr.cdr3". To set the kmer size to size 5 use encoding = "5mer" for example. If only "kmer" is given the default size is 3. The default encoding method overall is set to "onehot".
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param unique.sequences Level at which unique sequences are filtered. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain amino acid sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @return Returns encoded features
#' @export
#' @examples
#' \dontrun{
#' check_encode_features <- encode_features(features = output.load_data, encoding = "onehot", unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"),  to.use = NULL)
#' }
#'

encode_features <- function(features, encoding, unique.sequences, to.use) {

    require(tidyverse)
    require(Biostrings)
    require(protr)

    if (missing(unique.sequences)) unique.sequences <- c("aa_sequence_HC", "aa_sequence_LC")

    f_temp <- features %>% bind_rows %>% drop_na(matches("cdr3_nt_[H,L]C|sequence_[H,L]C")) %>%
        distinct(across(all_of(unique.sequences)), .keep_all = T)

    ELISA_bind <- f_temp %>% pull(ELISA_bind)
    octet <- f_temp %>% pull(octet)
    f_temp <- f_temp %>% select(-c(ELISA_bind, octet, barcode))

    # Selecting only some features if specified
    if (!missing(to.use)) f_temp <- f_temp %>% select(to.use, "cdr3s_aa")

    # One hot encoding gene information
    genes <- f_temp %>% select(matches("gene")) %>%
        sapply(., as.data.frame, stringsAsFactors = T) %>%
            sapply(., function(x) data.frame(sapply(x, as.numeric))) %>% as.data.frame

    colnames(genes) <- c("v_gene_HC", "d_gene_HC",  "j_gene_HC", "c_gene_HC",
                         "v_gene_LC", "d_gene_LC",  "j_gene_LC", "c_gene_LC")

    # Selecting sequence features
    f_temp <- f_temp %>% select(matches("cdr|sequence"))

    if (encoding == "onehot") {
        f_encoded <- lapply(f_temp, function(x)
            data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))
        f_encoded <- lapply(f_encoded, function(x) data.frame(sapply(x, as.numeric)))
        f_encoded <- cbind.data.frame(f_encoded)
        f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]
    }

    if (str_detect(encoding, "[0-9,k]mer")) {
        if (str_starts(encoding, "k")) num <- 3 else num <- as.numeric(substr(encoding, 1, 1))
        f_encoded <- f_temp %>% select(matches("cdr3_nt_[H,L]C|$sequence_[H,L]C")) %>%
            lapply(function(x) as.data.frame(oligonucleotideFrequency(DNAStringSet(x), num))) %>% bind_cols
    }

    if (encoding == "protr") {
        f_encoded_HC <- map(gsub("\\*", "", f_temp$aa_sequence_HC), extractCTDC) %>% as.data.frame %>% t
        f_encoded_LC <- map(gsub("\\*", "", f_temp$aa_sequence_LC), extractCTDC) %>% as.data.frame %>% t
        f_encoded <- cbind(f_encoded_HC, f_encoded_LC, by=0, all=T) %>% as.data.frame
        colnames(f_encoded) <- seq(ncol(f_encoded))
    }

    if (encoding == "protr.cdr3") {
        f_encoded <- t(sapply(gsub("IG[H,K,L]:|;|:|", "", f_temp$cdr3s_aa), extractTC))
    }

    f_encoded <- cbind(f_encoded, genes, ELISA_bind, octet)

    return(f_encoded)
}
