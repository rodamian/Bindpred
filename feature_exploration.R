#' Plots repertoire features for binding versus non binding cells. This can be done both at the clonotype level or at the single cell level.
#' @param features List of dataframes containing the extracted features. This is the output of load_data function.
#' @param clono.level Logical indicating if the classified sequences are at the single cell level or at the clonotype level. Default is set to FALSE
#' @return This function plots the mean mutation for each cell (divided in Heavy chain and Light chain mutations), gene usage, CDR3 length, and the ratio between coding versus non-coding mutations.
#' @export
#' @examples
#' \dontrun{
#' check_explore_features <- explore_features(features = output.load_data, clono.level = FALSE)
#' }
#'

explore_features <- function(features, clono.level) {

  require(ggplot2); theme_set(theme_bw())
  require(tidyverse)

  features <- bind_rows(features) %>% rename("mouse_number" = "mouse") %>%
    filter(ELISA_bind %in% c("yes", "no")) %>% select(which(colMeans(is.na(.)) < 0.5))

  # Mutations data ----------------------------------------------------------

  # Calculation mutation mean for each group
  mut_matrix <- features %>% filter(!is.na(ELISA_bind)) %>%
   filter(ELISA_bind %in% c("yes", "no")) %>%
    select(mouse, ELISA_bind, mutation_count, mutation_count_LC) %>%
      group_by(ELISA_bind, mouse) %>%
        mutate(mean_HC = mean(mutation_count, na.rm=T)) %>%
          mutate(mean_LC = mean(mutation_count_LC, na.rm=T))

  # Plots histogram of mutations for each mouse separately (HC)
  ggplot(mut_matrix, aes(x=mutation_count, fill=ELISA_bind)) +
    geom_histogram(binwidth=.5, alpha=.5, position="identity", stat="count") +
      geom_vline(aes(xintercept=mean_HC,  colour=ELISA_bind),
                 linetype="dashed", size=1) + facet_wrap(~mouse)

  # Plots histogram of mutations for each mouse separately (LC)
  ggplot(mut_matrix, aes(x=mutation_count_LC, fill=ELISA_bind)) +
     geom_histogram(binwidth=.5, alpha=.5, position="identity", stat="count") +
     geom_vline(aes(xintercept=mean_LC,  colour=ELISA_bind),
                linetype="dashed", size=1) + facet_wrap(~mouse)


  # gene usage  -------------------------------------------------------------
  gene_matrix <- features  %>% select(matches("_gene|mouse|ELISA_bind")) %>%
        group_by(ELISA_bind, mouse)

  # C gene
  ggplot(gene_matrix, aes(x=mouse, fill=c_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + facet_wrap(~ELISA_bind)

  # V gene
  ggplot(gene_matrix, aes(x=mouse, fill=v_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + facet_wrap(~ELISA_bind)

  # J gene
  ggplot(gene_matrix, aes(x=mouse, fill=j_gene_HC)) +
    geom_bar(alpha=.8, position="stack", stat="count") + facet_wrap(~ELISA_bind)


  # CDR3 length
  clono_matrix <- features %>% filter(!is.na(ELISA_bind)) %>%
    filter(ELISA_bind %in% c("yes", "no")) %>% group_by(ELISA_bind, mouse) %>%
      mutate(cdr3_len = str_length(cdr3s_aa), mean_len = mean(str_length(cdr3s_aa)))

  ggplot(clono_matrix, aes(x = cdr3_len, fill = ELISA_bind)) +
    geom_histogram(alpha = .5, position = "identity", stat="count") +
      geom_vline(aes(xintercept = mean_len, colour = ELISA_bind), linetype = "dashed", size = 1) +
                   facet_wrap(~mouse)

  # Degenerate mutations analysis
  features <- features %>% group_by(clonotype_id) %>%
    mutate(mean_umis_HC = mean(umis_HC), mean_umis_LC = mean(umis_LC)) %>%
    mutate(deg_ratio = length(unique(cdr3_HC)) + length(unique(cdr3_LC)) / length(unique(cdr3_nt_HC)) + length(unique(cdr3_nt_LC)))

  ggplot(features, aes(x=deg_ratio, y=umis_HC)) + geom_point() + geom_smooth()
  ggplot(features, aes(x=deg_ratio, y=umis_LC)) + geom_point() + geom_smooth()

  ggplot(features, aes(x=ELISA_bind, y=deg_ratio)) + geom_col()

}


explore_features(features)
