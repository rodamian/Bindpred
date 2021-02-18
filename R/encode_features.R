#' Encodes repertoire features with different methods. Gene usage features are encoded using one hot encoding irrespective of the encoding method used. Positions of insertions deletions and mutations are provided if not specified otherwise. Sequence features can be encoded using multiple methods.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param encoding Character indicating which encoding strategy to use. Options are "onehot", "kmer", "protr", "protr.cdr3". To set the kmer size to size 5 use encoding = "5mer" for example. If only "kmer" is given the default size is 3. The default encoding method overall is set to "onehot".
#' @param unique.sequences Level at which unique sequences are filtered. Default is c("aa_sequence_HC", "aa_sequence_LC") which keeps every cell with a unique combination of heavy and light chain amino acid sequences. Other options include "clonotype_id", "cdr3s_aa", "aa_sequence_HC".
#' @param to.use Character vector indicating which features to use. If not supplied all the features will be used
#' @param filter.corr Numeric indicating the threshold of correlation after which features are dropped. Default is 0.9.
#' @param filter.corr Numeric indicating the minimum number of unique values per feature. Important when calculating the tripeptide composition for short sequences. Default is 3.
#' @return Returns encoded features
#' @export
#' @examples
#' \dontrun{
#' check_encode_features <- encode_features(features = output.load_data, encoding = "onehot", unique.sequences = c("aa_sequence_HC", "aa_sequence_LC"),  to.use = NULL)
#' }
#'

encode_features <- function(features,
                            encoding = "onehot",
                            unique.sequences = "cdr3s_aa",
                            to.use = c("cdr3s_aa", "cdr3s_nt", "aa_sequence_HC", "aa_sequence_LC"),
                            filter.corr = 0.9,
                            filter.unique.values = 4) {

    require(protr)
    require(Biostrings)
    require(caret)
    require(tidyverse)
    library(corrplot)

    f_temp <- features %>% bind_rows %>% distinct(across(all_of(unique.sequences)), .keep_all = T) %>%
        filter(ELISA_bind %in% c("yes", "no") | !is.na(octet)) %>%
            drop_na(any_of(c("aa_sequence_HC", "aa_sequence_LC")))

    # Selecting only some features if specified
    # f_temp <- f_temp %>% dplyr::select(all_of(to.use))

    ELISA_bind <- f_temp %>% pull(ELISA_bind)
    octet <- f_temp %>% pull(octet)
    f_temp <- f_temp %>% dplyr::select(-c(ELISA_bind, octet, barcode, clonotype_id))

    # One hot encoding gene information
    genes <- f_temp %>% dplyr::select(matches("gene")) %>%
        sapply(., as.data.frame, stringsAsFactors = T) %>%
            sapply(., function(x) data.frame(sapply(x, as.numeric))) %>% as.data.frame

    colnames(genes) <- c("v_gene_HC", "d_gene_HC",  "j_gene_HC", "c_gene_HC",
                         "v_gene_LC", "d_gene_LC",  "j_gene_LC", "c_gene_LC")

    # Creating features with mutation, insertion and deletion positions
    find_positions <- function(cigar, letter)
        sapply(replace_na(lapply(cigar, function(x) head(str_split(x, letter)[[1]], -1)), ""),
               function(x) cumsum(sapply(parse(text = gsub("\\D", "+", x)), eval)))

    mutations <- f_temp %>% transmute(del_pos_HC = find_positions(cigar_HC, "D"),
                             mut_pos_HC = find_positions(cigar_HC, "X"),
                             ins_pos_HC = find_positions(cigar_HC, "I"),
                             del_pos_LC = find_positions(cigar_LC, "D"),
                             mut_pos_LC = find_positions(cigar_LC, "X"),
                             ins_pos_LC = find_positions(cigar_LC, "I")) %>% select(-matches("cigar|gene")) %>%
        lapply(function(x) t(data.frame(lapply(x, "length<-", max(lengths(x)))))) %>% do.call("cbind" ,.) %>%
            replace_na(-1) %>% as.data.frame

    sequences <- f_temp %>% select(matches("aa|_nt")) %>% transmute_all(str_replace_all, "IGH:|;IG[K,L]|NA|:|;|\\*", "") %>%
        select(all_of(to.use))

    if (encoding == "onehot") {
        f_encoded <- lapply(sequences, function(x)
            data.frame(str_split_fixed(x, "", max(sapply(x, str_length))), stringsAsFactors = T))
        f_encoded <- lapply(f_encoded, function(x) data.frame(sapply(x, as.numeric)))
        f_encoded <- cbind.data.frame(f_encoded)
        f_encoded <- f_encoded[vapply(f_encoded, function(x) length(unique(x)) > 1, logical(1L))]
    }

    if (str_detect(encoding, "[0-9,k]mer")) {
        if (str_starts(encoding, "k")) num <- 3 else num <- as.numeric(substr(encoding, 1, 1))
        f_encoded <- sequences %>% dplyr::select(matches("nt")) %>%
            lapply(function(x) as.data.frame(oligonucleotideFrequency(DNAStringSet(x), num))) %>% bind_cols
    }

    if (encoding == "blosum") {
        f_encoded <- sequences %>% dplyr::select(matches("aa")) %>%
        lapply(function(x) map(x, extractBLOSUM, k = 5, lag = 4) %>% bind_cols %>% t)
    }

    if (encoding == "dc") {
        f_encoded <- sequences %>% dplyr::select(matches("aa")) %>%
            lapply(function(x) map(x, extractDC) %>% bind_cols %>% t)
    }

    if (encoding == "tc.cdr3") {
        f_encoded <- sequences %>% dplyr::select(matches("cdr3s_aa")) %>%
            lapply(function(x) map(x, extractTC) %>% bind_cols %>% t)
    }

    f_encoded <- do.call("cbind", f_encoded) %>% as.data.frame %>% cbind(genes, mutations)

    # Filtering constant columns
    f_encoded <- Filter(var, f_encoded)

    # Filtering columns with low numbers of unique values
    f_encoded <- Filter(function(x) n_distinct(x) >= filter.unique.values, f_encoded)

    # Correlation filter
    corr <- cor(f_encoded)
    highlyCor <- findCorrelation(corr, filter.corr)
    f_encoded <- f_encoded[,-highlyCor]

    f_encoded <- cbind(f_encoded, ELISA_bind, octet)
    rownames(f_encoded) <-  NULL

    return(f_encoded)
}
