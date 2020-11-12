library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)


bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, accel = X6, total = X7) %>%
  separate(combo, into = c(NA, "pass1"), sep = "GeneID:", remove = T) %>%
  separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) 

target <- bg %>% filter(accel >= 1) %>% select(-c(chr, start, end, gene, accel, total))

calc_enrich <- function(targetset, background, ont) { 
  enrichGO(targetset$ncbi,'org.Gg.eg.db',
           pvalueCutoff=1.5,
           qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
           pAdjustMethod="BH",
           universe=background$ncbi,
           keyType="ENTREZID",
           ont=ont) 
}

bp.perms <- data.frame()

for (i in target[,2:1001]) {
  perm.target <- filter(target, i >= 1) %>%
    select(ncbi)
  bp <- calc_enrich(perm.target, bg, ont = "BP")
  bp.clean <- bp@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/') %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/') %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, geneID, logp, target_frac, bg_frac, enrich,pvalue, qvalue) %>%
    arrange(ID) 
  bp.perms <- as.data.frame(rbind(bp.perms, bp.clean))
  
}

bp.perms.clean <- bp.perms %>%
  select(-geneID) %>%
  distinct() %>%
  arrange()

bp.real <- calc_enrich(target, bg, "BP")

bp.real.clean <- bp.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, geneID, logp, target_frac, bg_frac, enrich, pvalue, qvalue) %>%
  distinct() %>%
  arrange(ID)

bp.merge <- bp.real.clean %>% 
  full_join(., bp.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
  distinct() %>%
  arrange(ID)

  
  
  
  
  
  
