library(tidyverse)
setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run")

#### spatial enrichment

# load in accelerated CNEEs by window
acc.cnees <- read_delim("window.acc.cnees.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, acc = X4)
# convert the CNEE name to a binary measure of occurrence (1)
acc.convert <- str_replace_all(acc.cnees$acc, regex("zfCNEE[0-9]{0,6}"), "1")
# summarise number of accelerated CNEEs by window
acc.clean <- bind_cols(acc.cnees, acc.convert) %>% 
  select(-c(acc)) %>% 
  rename(acc = ...5) %>%
  mutate(cnee = as.numeric(acc)) %>% 
  group_by(chr, start, end) %>%
  summarise(acc.cnee = sum(cnee))

# load in total list of CNEEs by window
all.cnees <- read_delim("window.cnees.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4)
# convert CNEE names to a binary measure of occurrence (1)
all.convert <- str_replace_all(all.cnees$cnee, regex("zfCNEE[0-9]{0,6}"), "1")
# summarise total number of CNEEs by window
all.clean <- bind_cols(all.cnees, all.convert) %>% 
  select(-c(cnee)) %>% 
  rename(cnee = ...5) %>%
  mutate(cnee = as.numeric(cnee)) %>%
  group_by(chr, start, end) %>%
  summarise(total.cnee = sum(cnee))

# combine & clean
data <- full_join(all.clean, acc.clean, by = c("chr", "start", "end"))
final.data <- data %>%
  unite(window, start, end, sep = "-", remove = T) %>%
  filter(total.cnee > 0)

# binomial test function
bt <- function(x, n, p = 294/375606) {
  binom.test(x, n, 294/375606, alternative = "greater", conf.level = 0.95)$p.value
}

# add binomial p-values to table
final.data$pVal <- mapply(bt, final.data$acc.cnee, final.data$total.cnee)

# any windows with sig more acc CNEEs than would be expected? 359!
final.data %>% 
  filter(pVal < 0.05) %>%
  write_delim("final.cnees.sig.txt", delim = "\t", col_names = T)

# collapses 359 significant windows into 21 chromosomes
final.data %>%
  filter(pVal < 0.05) %>%
  group_by(chr) %>% 
  summarise(ntotal = sum(total.cnee),
            nacc = sum(acc.cnee)) %>%
  write_delim("final.cnees.sigbychr.txt", delim = "\t", col_names = T)


#### GO term enrichment

library(tidyverse)
library(permute)
library(biomaRt)
library(clusterProfiler)

setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run")

genes <- read_delim("galGal_cnees_genes.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4, gene = X5)

acc <- read_delim("acc.cnees.final.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4)

data <- left_join(genes, acc, by = "cnee") %>%
  replace_na(list(chr.y = 0, start.y = 0, end.y = 0)) %>%
  select(-c(start.y, end.y)) %>%
  rename(chr = chr.x, start = start.x, end = end.x, temp = chr.y)
data.convert <- gsub(regex("NC_[0-9]{0,6}.[0-9]{1}"), "1", data$temp)
data.clean <- bind_cols(data, data.convert) %>%
  select(-c(temp)) %>%
  rename(binary = ...7) %>%
  mutate(binary = as.numeric(binary)) %>%
  group_by(gene) %>%
  summarise(accel = sum(binary)) %>%
  slice(-1) # there are 11 CNEE entries that have . as the gene name? Need to investigate 
data.sub <- gsub(regex("[2-9]{1,2}"), "1", data.clean$accel)
data.final <- bind_cols(data.clean, data.sub) %>%
  select(-c(accel)) %>%
  rename(accel = ...3)

# background permutations 
bg.perms <- lapply(1:100, function(x){
  set.seed(42)
  cbind.data.frame(data.final[,1], sample(data.final[,2]))
})

# proportion of genes with accelerated CNEEs = 0.02108621 
data.final %>% 
  filter(accel == 1) %>% 
  count()
data.final %>% 
  filter(accel == 0) %>% 
  count()

# gene IDs
mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
geneList <- unique(sort(data.final$gene))
test <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = geneList, bmHeader = T, mart = mart)


