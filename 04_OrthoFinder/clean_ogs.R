library(tidyverse)

data <- read_delim("run_ortho/results/Results_Nov09/Orthogroups/Orthogroups.GeneCount.tsv", delim = '\t')

data %>%
  filter(Total <= 25) %>%
  filter(hetAtr.translated != 0) %>%
  rowwise() %>%
  mutate(missing_spp = sum(c_across(anaPla_protein:tymCupPin.translated) == 0)) %>%
  filter(missing_spp <= 7) %>%
  select(Orthogroup) %>%
  write_delim(., "clean_ogs_geneCount.tsv", delim = '\t')
