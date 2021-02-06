#' Plots repertoire features for binding versus non binding cells. This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param per.sample Logical indicating if the plots produced are separated by sample or aggregated.
#' @return This function plots the mean mutation for each cell (divided in Heavy chain and Light chain mutations), gene usage, CDR3 length, and the ratio between coding versus non-coding mutations.
#' @export
#' @examples
#' \dontrun{
#' check_explore_features <- explore_features(features = output.load_data)
#' }
#'

explore_features <- function(features, per.sample) {

  require(ggplot2); theme_set(theme_bw())
  require(tidyverse)
  update_geom_defaults("smooth", list(alpha = 0.25))

  f_temp <- bind_rows(features) %>% filter(ELISA_bind %in% c("yes", "no"))

  if (missing(per.sample)) per.sample <- T
  if (per.sample) f_temp <- f_temp %>% group_by(sample)

  output.plot <- list()

  if ("barcode" %in% colnames(f_temp)) {

    # Mutations data
    if (any(grepl("mutation", colnames(f_temp)))) {
    # Mutation mean for each group
    mut_matrix <- f_temp %>% filter(mutation_count_HC < 40) %>%  filter(mutation_count_LC < 20) %>%
        select(sample, ELISA_bind, mutation_count_HC, mutation_count_LC) %>%
            group_by(ELISA_bind) %>%
              mutate(mean_HC = mean(mutation_count_HC, na.rm=T), mean_LC = mean(mutation_count_LC, na.rm=T))

    # Plots histogram of mutations in the heavy chain (_HC)
    output.plot[["mutation_HC"]] <- ggplot(mut_matrix, aes(x=mutation_count_HC, fill=ELISA_bind)) +
      geom_histogram(binwidth=.5, alpha=.5, position="identity", stat="count") +
        geom_vline(aes(xintercept=mean_HC, colour=ELISA_bind), linetype="dashed", size=1) +
          geom_vline(aes(xintercept=mean_HC, colour=ELISA_bind), linetype="dashed", size=1) +
            {if (per.sample) facet_wrap(~sample)}

    # Plots histogram of mutations in the light chain (_LC)
    output.plot[["mutation_LC"]] <- ggplot(mut_matrix, aes(x=mutation_count_LC, fill=ELISA_bind)) +
      geom_histogram(binwidth=.5, alpha=.5, position="identity", stat="count") +
        geom_vline(aes(xintercept=mean_LC,  colour=ELISA_bind), linetype="dashed", size=1) +
          geom_vline(aes(xintercept=mean_LC,  colour=ELISA_bind), linetype="dashed", size=1) +
            {if (per.sample) facet_wrap(~sample)}
    }

  # gene usage  -------------------------------------------------------------
  gene_matrix <- f_temp %>% select(matches("_gene|sample|ELISA_bind")) %>% group_by(ELISA_bind)

  output.plot[["c_gene"]] <- ggplot(gene_matrix, aes(x=ELISA_bind, fill=c_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + {if (per.sample) facet_wrap(~sample)}

  output.plot[["v_gene"]] <- ggplot(gene_matrix, aes(x=ELISA_bind, fill=v_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + {if (per.sample) facet_wrap(~sample)}

  output.plot[["j_gene"]] <- ggplot(gene_matrix, aes(x=ELISA_bind, fill=j_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + {if (per.sample) facet_wrap(~sample)}


  # CDR3 length -------------------------------------------------------------
  clono_matrix <- f_temp %>% group_by(ELISA_bind) %>%
      mutate(cdr3_len = str_length(gsub("IG[H,K,L]:|;|:|NA", "", cdr3s_aa)),
             mean_len = mean(str_length(gsub("IG[H,K,L]:|;|:|NA", "", cdr3s_aa))))

  output.plot[["cdr3_length"]] <- ggplot(clono_matrix, aes(x=cdr3_len, fill=ELISA_bind)) +
    geom_histogram(alpha = .5, position = "identity", stat="count") +
      geom_vline(aes(xintercept=mean_len,  colour=ELISA_bind), linetype="dashed", size=1) +
        {if (per.sample) facet_wrap(~sample)}


  # Degenerate mutations analysis # --------------------------------------------
  deg_data <- f_temp %>% group_by(clonotype_id) %>%
    mutate(mean_umis_HC = mean(umis_HC), mean_umis_LC = mean(umis_LC)) %>%
      mutate(deg_ratio_HC = length(unique(aa_sequence_HC))/length(unique(sequence_HC))) %>%
        mutate(deg_ratio_LC = length(unique(aa_sequence_LC))/length(unique(sequence_LC))) %>%
          mutate(deg_ratio = deg_ratio_HC + deg_ratio_LC)

  output.plot[["deg_ratio_HC"]] <- ggplot(deg_data, aes(x=ELISA_bind, y=deg_ratio_HC)) + geom_violin() + geom_boxplot(width=0.06, color="red") + {if (per.sample) facet_wrap(~sample)}
  output.plot[["deg_ratio_LC"]] <- ggplot(deg_data, aes(x=ELISA_bind, y=deg_ratio_LC)) + geom_violin() + geom_boxplot(width=0.06, color="red") + {if (per.sample) facet_wrap(~sample)}
  output.plot[["deg_ratio"]] <- ggplot(deg_data, aes(x=ELISA_bind, y=deg_ratio)) + geom_violin() + geom_boxplot(width=0.06, color="red") + {if (per.sample) facet_wrap(~sample)}

  }

    # Affinity analysis
    f_temp <- features %>% bind_rows %>% filter(!is.na(octet))
    if (any(grep("mutation", colnames(f_temp)))) {f_temp <- f_temp %>% filter(mutation_count_HC < 40, mutation_count_LC < 20)
        if (any(grep("octet", colnames(f_temp)))) {

          output.plot[["affinity_mutations_HC"]] <- ggplot(f_temp, aes(x=mutation_count_HC, y=octet, color=sample)) +
            scale_y_log10() + scale_x_log10() + geom_point() + geom_smooth(method="lm", aes(fill=sample)) + {if (per.sample) facet_wrap(~sample)}

          output.plot[["affinity_mutations_LC"]] <- ggplot(f_temp, aes(x=mutation_count_LC, y=octet, color=sample)) +
            scale_y_log10() + scale_x_log10() + geom_point() + geom_smooth(method="lm", aes(fill=sample)) + {if (per.sample) facet_wrap(~sample)}

          output.plot[["affinity_length"]] <- ggplot(f_temp, aes(x=length, y=octet, color=sample)) +
            scale_y_log10() + scale_x_log10() + geom_point() + geom_smooth(method="lm", aes(fill=sample)) + {if (per.sample) facet_wrap(~sample)}

          n_variants <- f_temp %>% group_by(clonotype_id) %>%
            mutate(variants_HC = n_distinct(aa_sequence_HC), variants_LC = n_distinct(aa_sequence_LC))

          output.plot[["variants"]] <- ggplot(n_variants, aes(x=variants_HC, y=octet, color=sample)) +
            scale_y_log10() + scale_x_log10() +
              geom_point() + geom_smooth(method="lm", aes(fill=sample)) + {if (per.sample) facet_wrap(~sample)}
        }
    }
    return(output.plot)
}
