library(tidyverse)

data <- read_delim("hetAtr_coverage_sites_clean_merged.bed", delim = '\t', col_names = F) %>%
  rename(chr = X1, start = X2, end = X3) %>%
  mutate(int_len = end-start)
sum(data$int_len)
