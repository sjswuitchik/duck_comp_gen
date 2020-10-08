## uploading for Oct 8 meeting

library(tidyverse)
library(coin)
library(permute)
library(biomaRt)
library(clusterProfiler)

setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run")

genes.perms <- read_delim("galGal_cnees_genes.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4, gene = X5)

acc.perms <- read_delim("acc.cnees.final.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4) %>%
  mutate(accel = "1")

data <- left_join(genes.perms, acc.perms, by = "cnee") %>%
  replace_na(list(chr.y = 0, start.y = 0, end.y = 0)) %>%
  select(-c(start.y, end.y, chr.y)) %>%
  rename(chr = chr.x, start = start.x, end = end.x)

data.convert <- replace_na(data$accel, 0)

data.clean <- bind_cols(data, data.convert) %>%
  select(-c(accel)) %>%
  rename(accel = ...7) %>%
  mutate(accel = as.numeric(accel),
         total = 1) %>%
  mutate(binary = accel/total) %>%
  group_by(gene) %>%
  summarise(accel = sum(accel),
            total = sum(total),
            prop = accel/total) %>%
  slice(-1) # there are 11 CNEE entries that have . as the gene name? Need to investigate 

prop <- function(x, accel, total) {
  (x[accel])/(x[total])
}

store <- data.frame(matrix(nrow = nrow(data.clean), ncol = 100))
nsim <- 100
set.seed(42)
for (i in 1:nsim) {
  perm <- shuffle(N)
  df <- data.clean %>%
    group_by(gene) %>%
    summarise(accel = sum(accel),
              accel.perm = sum(sample(accel)),
              total = sum(total),
              prop.perm = accel.perm/total)
  store[i] <- df$prop.perm
}

perm.data <- bind_cols(data.clean$prop, store) %>%
  rename(obs.prop = '...1') %>% 
  #bind_cols(data.clean$accel, .) %>%
  #rename(accel.cnee = '...1') %>%
  bind_cols(data.clean$gene, .) %>%
  rename(gene = '...1')

hist.data <- bind_cols(data.clean$prop, store) %>%
  rename(obs.prop = '...1')
d <- melt(hist.data, id.vars = "obs.prop")

ggplot(d, aes(x = value)) + 
  geom_histogram(position = 'identity', binwidth = 0.05) 


# permutation test of independence 
independence_test(obs.prop ~ value, data = d, alternative = "greater", distribution = "approximate")




#### ####

# proportion of genes with accelerated CNEEs 
data.clean %>% 
  count(accel >= 1)

prop = 238/11287 # 0.02108621

# gene IDs
mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
geneList <- unique(sort(data.clean$gene))
martList <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = geneList, bmHeader = T, mart = mart)

collapse <- martList %>% 
  rename(gene = 'Gene name', goid = 'GO term accession') %>%
  select(gene, goid) %>%
  group_by(gene) %>%
  summarise(goid = paste(sort(unique(goid)), collapse = ", "))

# remove genes without GO ids
sub1 <- collapse[!(is.na(collapse$goid) | collapse$goid == ""), ] 
clean.mart <- sub1[-1,]

# add in accel CNEEs
accel <- data.frame(data.clean$gene, data.clean$accel) %>%
  rename(gene = data.clean.gene, accel = data.clean.accel)
accel.mart <- left_join(accel, clean.mart, by = "gene") %>%
  na.omit()





