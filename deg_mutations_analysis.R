
# Degenerate mutations analysis

temp_files <- list()
for (file in list.files("features", "^cell", full.names = T)) {temp_files[[sub(".*/", "", file)]] <- read.csv(file)}

feature_matrix <- bind_rows(temp_files) %>% filter(!is.na(ELISA_bind)) %>%
  filter(ELISA_bind %in% c("yes", "no")) %>% select(-octet) %>% select(which(colMeans(is.na(.)) < 0.5))

temp_feat <- temp_feat %>% group_by(clonotype_id) %>%
  mutate(mean_umis_HC = mean(umis_HC), mean_umis_LC = mean(umis_LC)) %>%
    mutate(deg_ratio = length(unique(CDR3H)) + length(unique(CDR3L)) / length(unique(sequence_HC)) + length(unique(sequence_LC)))

ggplot(temp_feat, aes(x=deg_ratio, y=umis_HC)) + geom_point() + geom_smooth()
ggplot(temp_feat, aes(x=deg_ratio, y=umis_LC)) + geom_point() + geom_smooth()
