## Prelim phyloAcc analyses - NCBI run ## 

library(tidyverse)
library(ggthemes)
library(GGMM)
setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run")

# functions
trunc_to_one <- function(x) {
  ifelse(x > 1, 1, x)
}

discretize <- function(x, cutoff = 0.90) {
  ifelse(x >= cutoff, 1, 0)
}

#### top 1 ####

# read in and clean
likTop1 <- read_tsv("top1_combined_elem_lik.txt", col_names = c("key", "cnee", "loglik.null", "loglik.target", "loglik.full", "logBF1", "logBF2", "loglik_Max_M0", "loglik_Max_M1", "loglik_Max_M2")) %>%
  arrange(loglik.full) %>%
  distinct(cnee, .keep_all = T) %>%
  mutate(bf1 = as.numeric(loglik.target) - as.numeric(loglik.null), bf2 = as.numeric(loglik.target) - as.numeric(loglik.full)) %>%
  mutate(key = as.numeric(key))

zpostTop1 <- read_tsv("top1_combined_postZ_M2.txt", col_types = cols(.default = "d"))

# posterior prob acceleration
postaccTop1 <- zpostTop1 %>%
  rename(key = No.) %>%
  select(key, contains("_3"))
postaccTop1$key = zpostTop1$No.
postaccTop1$acc_rate = zpostTop1$n_rate
postaccTop1$cons_rate = zpostTop1$c_rate

## convert to 1/0 matrix
#postaccmatTop1 <- postaccTop1[1:32]
#postaccmatTop1[postaccmatTop1 >= 0.95] = 1
#postaccmatTop1[postaccmatTop1 < 0.95] = 0
#postaccmatTop1$key = postaccTop1$key

# fix names
names(postaccTop1) = gsub("-", "_", names(postaccTop1))

# create analysis data
cneeTop1 <- inner_join(postaccTop1, likTop1, by = c("key" = "key")) %>%
  filter(cons_rate <= 0.60) %>%
  ungroup

# add bed info
cneeBED <- read_tsv("galGal6_final_merged_CNEEs_named.bed", col_names = c("chr", "start", "end", "cnee"))
cneeTop1 <- inner_join(cneeTop1, cneeBED, by = c("cnee" = "cnee"))

# write out
write_csv(cneeTop1, "top1_cnees_analysis.csv")

cneeTop1 <- cneeTop1 %>%
  mutate(Accel = bf1 > 10, Spec = bf1 > 10 & bf2 > 1)

table(cneeTop1$Accel)["TRUE"] 
# 330 accel

#### top 2 ####

# read in and clean
likTop2 <- read_tsv("top2_combined_elem_lik.txt", col_names = c("key", "cnee", "loglik.null", "loglik.target", "loglik.full", "logBF1", "logBF2", "loglik_Max_M0", "loglik_Max_M1", "loglik_Max_M2")) %>%
  arrange(loglik.full) %>%
  distinct(cnee, .keep_all = T) %>%
  mutate(bf1 = as.numeric(loglik.target) - as.numeric(loglik.null), bf2 = as.numeric(loglik.target) - as.numeric(loglik.full)) %>%
  mutate(key = as.numeric(key))

zpostTop2 <- read_tsv("top2_combined_postZ_M2.txt", col_types = cols(.default = "d"))

# posterior prob acceleration
postaccTop2 <- zpostTop2 %>%
  rename(key = No.) %>%
  select(key, contains("_3"))
postaccTop2$key = zpostTop2$No.
postaccTop2$acc_rate = zpostTop2$n_rate
postaccTop2$cons_rate = zpostTop2$c_rate

# fix names
names(postaccTop2) = gsub("-", "_", names(postaccTop2))

# create analysis data
cneeTop2 <- inner_join(postaccTop2, likTop2, by = c("key" = "key")) %>%
  filter(cons_rate <= 0.60) %>%
  ungroup

# add bed info
cneeTop2 <- inner_join(cneeTop2, cneeBED, by = c("cnee" = "cnee"))

# write out
write_csv(cneeTop2, "top2_cnees_analysis.csv")

cneeTop2 <- cneeTop2 %>%
  mutate(Accel = bf1 > 10, Spec = bf1 > 10 & bf2 > 1)

table(cneeTop2$Accel)["TRUE"]
# 358 accel

#### top 3 ####

# read in and clean
likTop3 <- read_tsv("top3_combined_elem_lik.txt", col_names = c("key", "cnee", "loglik.null", "loglik.target", "loglik.full", "logBF1", "logBF2", "loglik_Max_M0", "loglik_Max_M1", "loglik_Max_M2")) %>%
  arrange(loglik.full) %>%
  distinct(cnee, .keep_all = T) %>%
  mutate(bf1 = as.numeric(loglik.target) - as.numeric(loglik.null), bf2 = as.numeric(loglik.target) - as.numeric(loglik.full)) %>%
  mutate(key = as.numeric(key))

zpostTop3 <- read_tsv("top3_combined_postZ_M2.txt", col_types = cols(.default = "d"))

# posterior prob acceleration
postaccTop3 <- zpostTop3 %>%
  rename(key = No.) %>%
  select(key, contains("_3"))
postaccTop3$key = zpostTop3$No.
postaccTop3$acc_rate = zpostTop3$n_rate
postaccTop3$cons_rate = zpostTop3$c_rate

# fix names
names(postaccTop3) = gsub("-", "_", names(postaccTop3))

# create analysis data
cneeTop3 <- inner_join(postaccTop3, likTop3, by = c("key" = "key")) %>%
  filter(cons_rate <= 0.60) %>%
  ungroup

# add bed info
cneeTop3 <- inner_join(cneeTop3, cneeBED, by = c("cnee" = "cnee"))

# write out
write_csv(cneeTop3, "top3_cnees_analysis.csv")

cneeTop3 <- cneeTop3 %>%
  mutate(Accel = bf1 > 10, Spec = bf1 > 10 & bf2 > 1)

table(cneeTop3$Accel)["TRUE"]
# 356 accel

#### explorations #### 
bf.1 <- cneeTop1 %>%
  select(cnee, chr, start, end, bf1, bf2, Accel)
bf.2 <- cneeTop2 %>%
  select(cnee, chr, start, end, bf1, bf2, Accel)
bf.3 <- cneeTop3 %>%
  select(cnee, chr, start, end, bf1, bf2, Accel)

all_bf <- inner_join(bf.3, bf.2, by = c("cnee" = "cnee"), suffix = c("3", "2")) %>%
  inner_join(bf.1, by = c("cnee")) %>%
  select(cnee, chr3, start3, end3, bf1, bf2, Accel, bf12, bf22, Accel2, bf13, bf23, Accel3) %>%
  rename(chr = chr3, start = start3, end = end3, bf1.1 = bf1, bf2.1 = bf2, acc.1 = Accel, bf1.2 = bf12, bf2.2 = bf22, acc.2 = Accel2, bf1.3 = bf13, bf2.3 = bf23, acc.3 = Accel3)

ggplot(all_bf, aes(x = bf1.3, y = bf1.2)) + 
  geom_point() + 
  labs(x = "Topology3", y = "Topology2")
ggplot(all_bf, aes(x = bf1.3, y = bf1.1)) + 
  geom_point() + 
  labs(x = "Topology3", y = "Topology1")
ggplot(all_bf, aes(x = bf1.1, y = bf1.2)) + 
  geom_point() +
  labs(x = "Topology1", y = "Topology2")

accelerated <- all_bf %>%
  filter(acc.1 == T, acc.2 == T, acc.3 == T)
# 295 CNEEs that are accel in all three tops
write_csv(accelerated, "acc_cnees.csv")

bed <- accelerated %>%
  mutate(newstart = start - 1) %>%
  select(chr, newstart, end, cnee) %>% 
  rename('#chrom' = chr, chromStart = newstart, chromEnd = end, name = cnee)
write_delim(bed, "acc_cnees.bed", delim = '\t')



