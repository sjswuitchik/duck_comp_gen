# remove all CNEEs that don't show strong evidence for conservation in Galloanserae from phyloP under all topologies after FDR correction

library(tidyverse)
setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run") #I know, sorry!

# read in all phyloP outputs for topologies, calculate FDR, create list of CNEEs to exclude
top1 <- read_delim("cnees_ncbi_phyloP_top1.out", delim = "\t", col_names = T) %>%
  rename(cnee = '#chr') %>%
  mutate(fdr = p.adjust(pval, method = "fdr"),
         excl = fdr < 0.10)
top1.excl <- top1 %>% 
  filter(excl == F) %>%
  select(cnee)

top2 <- read_delim("cnees_ncbi_phyloP_top2.out", delim = "\t", col_names = T) %>%
  rename(cnee = '#chr') %>%
  mutate(fdr = p.adjust(pval, method = "fdr"),
         excl = fdr < 0.10)
top2.excl <- top2 %>% 
  filter(excl == F) %>%
  select(cnee)

top3 <- read_delim("cnees_ncbi_phyloP_top3.out", delim = "\t", col_names = T) %>%
  rename(cnee = '#chr') %>%
  mutate(fdr = p.adjust(pval, method = "fdr"),
         excl = fdr < 0.10)
top3.excl <- top3 %>% 
  filter(excl == F) %>%
  select(cnee)

# intersect exclusionary CNEE lists for all topologies
cnee.excl.alltop <- intersect(top1.excl, top2.excl, top3.excl)

# read in BED of accelerated CNEEs
acc.cnee <- read_delim("acc_cnees.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4)
acc <- acc.cnee %>%
  select(cnee)

# intersect acc cnees with excl cnees to see if any need to be removed 
final.excl <- intersect(acc, cnee.excl.alltop) #zfCNEE300947

# filter out the one CNEE - would need to filter/setdiff based on final.excl object if more than a couple CNEEs
final.cnee <- acc.cnee %>%
  filter(cnee != "zfCNEE300947")

# write out new acc.cnees BED file
write_delim(final.cnee, "acc.cnees.final.bed", delim = "\t", col_names = F)


