library(tidyverse)

clean <- read_delim("clean_ogs.tsv", delim = '\t', col_names = c("Orthogroup"))
uniq <- read_delim("uniq", delim = '\t', col_names = c("Orthogroup"))

anti_join(clean, uniq, by = "Orthogroup") %>%
  write_delim(., "split", delim = '\t', col_names = F)
