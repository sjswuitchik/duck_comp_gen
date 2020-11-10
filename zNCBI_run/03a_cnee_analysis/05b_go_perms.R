#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)


bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, accel = X6, total = X7) %>%
  separate(combo, into = c(NA, "pass1"), sep = "GeneID:", remove = T) %>%
  separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
  select(-gene)

target <- bg %>% filter(accel >= 1) %>% select(-c(chr, start, end, accel, total))

calc_enrich <- function(targetset, background, ont) { 
  enrichGO(targetset$ncbi,'org.Gg.eg.db',
           pvalueCutoff=1.5,
           qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
           pAdjustMethod="none",
           universe=background$ncbi,
           keyType="ENTREZID",
           ont=ont) 
}

bp.perms <- data.frame()
mf.perms <- data.frame()
cc.perms <- data.frame()

for (i in 2:ncol(target)) {
  perm.target <- filter(target, i >= 1) %>%
    select(ncbi)
  bp <- calc_enrich(perm.target, bg, ont = "BP")
  mf <- calc_enrich(perm.target, bg, ont = "MF")
  cc <- calc_enrich(perm.target, bg, ont = "CC")
  bp.clean <- bp@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/') %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/') %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, geneID, logp, target_frac, bg_frac, OR, enrich) %>%
    arrange(ID) 
  mf.clean <- mf@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, geneID, logp, target_frac, bg_frac, OR, enrich) %>%
    arrange(ID) 
  cc.clean <- cc@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, geneID, logp, target_frac, bg_frac, OR, enrich) %>%
    arrange(ID) 
  bp.perms <- as.data.frame(rbind(bp.perms, bp.clean))
  mf.perms <- as.data.frame(rbind(mf.perms, mf.clean))
  cc.perms <- as.data.frame(rbind(cc.perms, cc.clean))
  
}

write_tsv(bp.perms, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/bp.perms.tsv", append = F, col_names = T)
write_tsv(mf.perms, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/mf.perms.tsv", append = F, col_names = T)
write_tsv(cc.perms, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cc.perms.tsv", append = F, col_names = T)

#bp.perms.load <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/bp.perms#.tsv", delim = "\t", )
#mf.perms.load <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/mf.perms#.tsv", delim = "\t", )
#cc.perms.load <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cc.perms#.tsv", delim = "\t", )

bp.real <- calc_enrich(target, bg, "BP")
#mf.real <- calc_enrich(target, bg, "MF")
#cc.real <- calc_enrich(target, bg, "CC")

bp.real.clean <- bp.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, geneID, logp, target_frac, bg_frac, OR, enrich) %>%
  filter(target_frac < 1, bg_frac < 1) %>%
  arrange(ID)

#mf.real.clean <- mf.real@result %>% 
#  separate(GeneRatio, into = c("target_in", "target_total")) %>%
#  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
#  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#         OR = target_frac/bg_frac,
#         enrich = log2(target_frac/bg_frac),
#         newpval = ifelse(is.na(pvalue), 1, pvalue),
#         logp = -log10(newpval)) %>%
#  dplyr::select(ID, logp, target_frac, bg_frac, OR, enrich) %>%
#  arrange(ID)
#
#cc.real.clean <- cc.real@result %>% 
#  separate(GeneRatio, into = c("target_in", "target_total")) %>%
#  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
#  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
#         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
#         OR = target_frac/bg_frac,
#         enrich = log2(target_frac/bg_frac),
#         newpval = ifelse(is.na(pvalue), 1, pvalue),
#         logp = -log10(newpval)) %>%
#  dplyr::select(ID, logp, target_frac, bg_frac, OR, enrich) %>%
#  arrange(ID)

bp.perms.merge <- bp.perms %>% 
  select(-c(geneID, OR)) %>%
  group_by(ID) %>%
  summarise(ecdf_frac = list(ecdf(target_frac)),
            ecdf_enrich = list(ecdf(enrich)),
            ecdf_logp = list(ecdf(logp)))


bp.real.merge <- bp.real.clean %>% select(-c(geneID, OR))

bp.merge <- left_join(bp.real.merge, bp.perms.merge, by = "ID")


# p values from this are all the same ... 
bp.merge.clean <-
  bp.merge %>% 
  rowwise %>% 
  mutate(pval_frac = max(1-ecdf_frac(target_frac), 0.0002), 
         pval_logp = max(1-ecdf_logp(logp), 0.0002), 
         pval_enrich = max(1-ecdf_enrich(enrich), 0.0002)) 

# adjust p values
pval_frac_adj <- data.frame(pval_frac = bp.merge.clean$pval_frac, qval_frac = p.adjust(bp.merge.clean$pval_frac, "BH"))
pval_logp_adj <- data.frame(pval_logp = bp.merge.clean$pval_logp, qval_logp = p.adjust(bp.merge.clean$pval_logp, "BH"))
pval_enrich_adj <- data.frame(pval_enrich = bp.merge.clean$pval_enrich, qval_enrich = p.adjust(bp.merge.clean$pval_enrich, "BH"))

plot(pval_frac_adj$pval_frac, pval_frac_adj$qval_frac)






