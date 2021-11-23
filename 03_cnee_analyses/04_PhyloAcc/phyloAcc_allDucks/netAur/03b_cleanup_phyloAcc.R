setwd("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/netAur")

library(tidyverse)

#### Top1 ####

col_elem <- c("num", "cnee", "lik_null", "lik_acc", "lik_full", "logBF1", "logBF2", "lik_maxM0", "lik_maxM1", "lik_maxM2")
top1_lik <- read_delim("top1_combined_elem_lik.txt.gz", col_names = col_elem, delim = "\t")

bed <- read_delim("galGal6_final_merged_CNEEs_named.bed", col_names = c("chr", "start", "end", "name"), delim="\t")

all_lik <- top1_lik %>%
  full_join(bed, by = c("cnee" = "name")) %>%
  mutate(length = end - start, num = as.character(num))

table(all_lik$lik_null == 0)

all_lik <- all_lik %>%
  filter(length >= 50, lik_null !=0)

write_tsv(all_lik, "top1_all_element_lik.txt.gz")

# get M2 results

all_cnee <- all_lik %>%
  select(num, cnee)

col_Z <- scan("rate_postZ_M2_top1.header", character(), quote = "")
col_Z[1] = "cnee_num"

top1Z <- read_delim("top1_combined_postZ_M2.txt.gz", col_names = col_Z, delim = "\t") %>%
  mutate(cnee_num = as.character(cnee_num)) %>%
  right_join(all_cnee, by = c("cnee_num" = "num"))

top1Z %>% write_tsv("top1_Zmat.txt.gz")

#### Top 2 ####

top2_lik <- read_delim("top2_combined_elem_lik.txt.gz", col_names = col_elem, delim = "\t")

all_lik <- top2_lik %>%
  full_join(bed, by = c("cnee" = "name")) %>%
  mutate(length = end - start, num = as.character(num))

table(all_lik$lik_null == 0)

all_lik <- all_lik %>%
  filter(length >= 50, lik_null !=0)

write_tsv(all_lik, "top2_all_element_lik.txt.gz")

# get M2 results

all_cnee <- all_lik %>%
  select(num, cnee)

col_Z <- scan("rate_postZ_M2_top2.header", character(), quote = "")
col_Z[1] = "cnee_num"

top2Z <- read_delim("top2_combined_postZ_M2.txt.gz", col_names = col_Z, delim = "\t") %>%
  mutate(cnee_num = as.character(cnee_num)) %>%
  right_join(all_cnee, by = c("cnee_num" = "num"))

top2Z %>% write_tsv("top2_Zmat.txt.gz")

#### Top 3 ####

top3_lik <- read_delim("top3_combined_elem_lik.txt.gz", col_names = col_elem, delim = "\t")

all_lik <- top3_lik %>%
  full_join(bed, by = c("cnee" = "name")) %>%
  mutate(length = end - start, num = as.character(num))

table(all_lik$lik_null == 0)

all_lik <- all_lik %>%
  filter(length >= 50, lik_null !=0)

write_tsv(all_lik, "top3_all_element_lik.txt.gz")

# get M2 results

all_cnee <- all_lik %>%
  select(num, cnee)

col_Z <- scan("rate_postZ_M2_top3.header", character(), quote = "")
col_Z[1] = "cnee_num"

top3Z <- read_delim("top3_combined_postZ_M2.txt.gz", col_names = col_Z, delim = "\t") %>%
  mutate(cnee_num = as.character(cnee_num)) %>%
  right_join(all_cnee, by = c("cnee_num" = "num"))

top3Z %>% write_tsv("top3_Zmat.txt.gz")
