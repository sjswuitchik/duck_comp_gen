library(tidyverse)

data <- read_delim("transGene.txt", delim = '\t', col_names = F) %>%
  select(X2, X1) %>% 
  write_delim(., "transGene_final.tsv", delim = '\t', col_names = F)
