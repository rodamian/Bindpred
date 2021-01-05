
### Feature exploration and plotting

# Figure block 1. Experimental overview of OVA dataset and schematic of feature matrix extraction from single-cell
# immune repertoire sequencing
#
# Figure block 2. Repertoire features of binders vs non binders
# V gene usage, SHM, number of variants, mean isotype distribution within clonotype, cdr3 lengths.
# Here we can have summary statistics for clones - e.g. median, mean, standard deviation, max, min
#
# Figure block 3. Repertoire features of high affinity vs low affinity
# V gene usage, SHM, number of variants, mean isotype distribution within clonotype, cdr3 lengths.
# Here we can have summary statistics for clones - e.g. median, mean, standard deviation, max, min
#
# Figure block 4. Three encoding options in R
# (can also have accompanying python pipeline, but for me extremely preferable to have this in R as well
#   (or at least that I have the encoded matrix in R afterwards)
#    One hot encoding, kmer, blosum/pam
#    Figure block 5. Classification of OVA binder vs non-binder, Cell level, clonotype
#    KNN  Linear discriminant analysis CART  SVM - linear kernel  Random forests  Predict binding
#
#    Figure block 6. Affinity prediction
#    Linear regression
#    Logistic regression
#    Ridge regression
#    Lasso regression
#
#    Figure block 7. Scikit learn pipeline and more sophisticated ML models
#    —> Scikit learn pipeline, e.g. writing feature matrix, reading into python, training, predicting,
#    interpreting the results in either python or R.
#
#    Figure block 8. Adaptability to other datasets / framework/package overview

require(ggplot2)
library(ggplot2); theme_set(theme_bw())
require(tidyverse)

clono_level <- F

if (clono_level) level <- "^clono" else level <- "^cell"

temp_files <- list()
for (file in list.files("features", level, full.names = T)) {
 temp_files[[sub(".*/", "", file)]] <- read.csv(file)
}
feature_matrix <- bind_rows(temp_files)

# Mutations data ----------------------------------------------------------

# Calculation mutation mean for each group
mut_matrix <- feature_matrix %>% filter(!is.na(ELISA_bind)) %>%
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

gene_matrix <- feature_matrix %>% filter(!is.na(ELISA_bind)) %>%
  filter(ELISA_bind %in% c("yes", "no")) %>% select(-octet) %>%
     select(matches("_gene|mouse|ELISA_bind")) %>%
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

clono_matrix <- feature_matrix %>% filter(!is.na(ELISA_bind)) %>%
  filter(ELISA_bind %in% c("yes", "no")) %>% group_by(ELISA_bind, mouse) %>%
    mutate(cdr3_len = str_length(cdr3s_aa), mean_len = mean(str_length(cdr3s_aa)))

ggplot(clono_matrix, aes(x=cdr3_len, fill=ELISA_bind)) +
  geom_histogram(alpha=.5, position="identity", stat="count") +
    geom_vline(aes(xintercept=mean_len, colour=ELISA_bind), linetype="dashed", size=1) +
                 facet_wrap(~mouse)


# Degenerate codon usage --------------------------------------------------

feature_matrix

