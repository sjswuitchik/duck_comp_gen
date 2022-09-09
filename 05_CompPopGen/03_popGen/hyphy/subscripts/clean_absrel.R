library(tidyverse)

abs <- read_delim("absrel_output_clean_parse_hetAtr_head.csv", delim = ',')
abs_clean <- abs %>%
  select(-c("baseline mean omega", "baseline median omega")) %>%
  rename("total_branches" = "num branches", "sig_branches" = "num branches pval less than alpha") %>%
  na.omit() %>%
  filter(sig_branches > 0, sig_branches <= 3) %>%
  separate('species:pval', into = c('branch01', 'branch02', 'branch03'), sep = ';', extra = 'drop')
write_delim(abs_clean, "abs_cleanR.csv", delim = ',')
