#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)
library(enrichplot)


bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, accel = X6, total = X7) %>%
  separate(combo, into = c(NA, "pass1"), sep = "GeneID:", remove = T) %>%
  separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
  select(-c(chr, start, end, gene, total))

calc_enrich <- function(targetset, background, ont) { 
  enrichGO(targetset$ncbi,'org.Gg.eg.db',
           pvalueCutoff=1.5,
           qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
           pAdjustMethod="BH",
           universe=background$ncbi,
           keyType="ENTREZID",
           ont=ont) 
}


bp.perms <- list()
mf.perms <- list()
cc.perms <- list()

for (i in 1:1000) {
  col_name <- paste0("X", i + 7)
  perm.target <- bg %>% filter(.data[[col_name]] >= 1) %>%
    select(ncbi, .data[[col_name]])
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
    dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  mf.clean <- mf@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/') %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/') %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  cc.clean <- cc@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/') %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/') %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  bp.perms[[i]] <- bp.clean
  mf.perms[[i]] <- mf.clean
  cc.perms[[i]] <- cc.clean
  
}

bp.perms.clean <- bind_rows(list(bp.perms), .id = "perm")
mf.perms.clean <- bind_rows(list(mf.perms), .id = "perm")
cc.perms.clean <- bind_rows(list(cc.perms), .id = "perm")

target <- bg %>% filter(accel >= 1)

bp.real <- calc_enrich(target, bg, "BP")
mf.real <- calc_enrich(target, bg, "MF")
cc.real <- calc_enrich(target, bg, "CC")

bp.real.clean <- bp.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
  arrange(ID)
mf.real.clean <- mf.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
  arrange(ID)
cc.real.clean <- cc.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, logp, target_frac, bg_frac, enrich) %>%
  arrange(ID)

bp.merge <- bp.real.clean %>% 
  left_join(., bp.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
  distinct() %>%
  arrange(ID) %>%
  select(-perm)
mf.merge <- mf.real.clean %>% 
  left_join(., mf.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
  distinct() %>%
  arrange(ID) %>%
  select(-perm)
cc.merge <- cc.real.clean %>% 
  left_join(., cc.perms.clean, by = c("ID" = "ID"), suffix = c(".real", ".perm")) %>% 
  distinct() %>%
  arrange(ID) %>%
  select(-perm)

bp.p <- bp.merge %>%
  group_by(ID) %>%
  mutate(gt_bg = bg_frac.perm >= bg_frac.real,
         gt_target = target_frac.perm >= target_frac.real,
         gt_enrich = enrich.perm >= enrich.real)
bp.cols <- sapply(bp.p[,10:12], as.numeric)
bp.pval <- bp.p %>%
  group_by(ID) %>%
  mutate(sum_bg = sum(gt_bg) + 1,
         pVal_bg = sum_bg/1001,
         sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  filter(pVal_bg <= 0.05 | pVal_target <= 0.05 | pVal_enrich <= 0.05) %>%
  select(ID, pVal_bg, pVal_target, pVal_enrich) %>%
  distinct()

mf.p <- mf.merge %>%
  group_by(ID) %>%
  mutate(gt_bg = bg_frac.perm >= bg_frac.real,
         gt_target = target_frac.perm >= target_frac.real,
         gt_enrich = enrich.perm >= enrich.real)
mf.cols <- sapply(mf.p[,10:12], as.numeric)
mf.pval <- mf.p %>%
  group_by(ID) %>%
  mutate(sum_bg = sum(gt_bg) + 1,
         pVal_bg = sum_bg/1001,
         sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  filter(pVal_bg <= 0.05 | pVal_target <= 0.05 | pVal_enrich <= 0.05) %>%
  select(ID, pVal_bg, pVal_target, pVal_enrich) %>%
  distinct()

cc.p <- cc.merge %>%
  group_by(ID) %>%
  mutate(gt_bg = bg_frac.perm >= bg_frac.real,
         gt_target = target_frac.perm >= target_frac.real,
         gt_enrich = enrich.perm >= enrich.real)
cc.cols <- sapply(cc.p[,10:12], as.numeric)
cc.pval <- cc.p %>%
  group_by(ID) %>%
  mutate(sum_bg = sum(gt_bg) + 1,
         pVal_bg = sum_bg/1001,
         sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  filter(pVal_bg <= 0.05 | pVal_target <= 0.05 | pVal_enrich <= 0.05) %>%
  select(ID, pVal_bg, pVal_target, pVal_enrich) %>%
  distinct()

write_delim(bp.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/sigterms_BP.tsv", delim = "\t", col_names = T)
write_delim(mf.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/sigterms_MF.tsv", delim = "\t", col_names = T)
write_delim(cc.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/sigterms_CC.tsv", delim = "\t", col_names = T)
