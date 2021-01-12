#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)

# load in background data and clean up NCBI IDs
bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, accel = X6, total = X7) %>%
  separate(combo, into = c(NA, "pass1"), sep = "GeneID:", remove = T) %>%
  separate(pass1, into = c("ncbi", NA), sep = "miRBase:", remove = T) %>%
  select(-c(chr, start, end, gene, total))

# set up enrichment calculation
calc_enrich <- function(targetset, background, ont) { 
  enrichGO(targetset$ncbi,'org.Gg.eg.db',
           pvalueCutoff=1.5,
           qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
           pAdjustMethod="BH",
           universe=background$ncbi,
           keyType="ENTREZID",
           ont=ont) 
}


# initialize empty lists for storing permutation results
bp.perms <- list()
mf.perms <- list()
cc.perms <- list()

# loop through every column (each column is a permutation), select target set for given permutation, run through calc_enrich for each subontology, calculate odds ratio and enrichment scores from gene and bg ratios, and store results in list
for (i in 1:1000) {
  col_name <- paste0("X", i + 7)
  perm.target <- bg %>% filter(.data[[col_name]] >= 1) %>%
    select(ncbi, .data[[col_name]])
  bp <- calc_enrich(perm.target, bg, ont = "BP")
  mf <- calc_enrich(perm.target, bg, ont = "MF")
  cc <- calc_enrich(perm.target, bg, ont = "CC")
  bp.clean <- bp@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  mf.clean <- mf@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  cc.clean <- cc@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
    separate(BgRatio, into = c("bg_in", "bg_total"), sep = '/', remove= F) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac,
           enrich = log2(target_frac/bg_frac),
           newpval = ifelse(is.na(pvalue), 1, pvalue),
           logp = -log10(newpval)) %>%
    dplyr::select(ID, GeneRatio, BgRatio, logp, target_frac, bg_frac, enrich) %>%
    arrange(ID)
  bp.perms[[i]] <- bp.clean
  mf.perms[[i]] <- mf.clean
  cc.perms[[i]] <- cc.clean
  
}


# bind permutation results
bp.perms.clean <- bind_rows(list(bp.perms), .id = "perm")
mf.perms.clean <- bind_rows(list(mf.perms), .id = "perm")
cc.perms.clean <- bind_rows(list(cc.perms), .id = "perm")

# write out perms so you don't have to re-run when your session inevitably crashes
write_delim(bp.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/bp.perms.clean.tsv", delim = "\t", col_names = T)
write_delim(mf.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/mf.perms.clean.tsv", delim = "\t", col_names = T)
write_delim(cc.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/cc.perms.clean.tsv", delim = "\t", col_names = T)

# read in permutation results if working in a new session
bp.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/bp.perms.clean.tsv", col_names = T)
mf.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/mf.perms.clean.tsv", col_names = T)
cc.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/cc.perms.clean.tsv", col_names = T)

# filter out observed target set from background data
target <- bg %>% filter(accel >= 1)

# calculate enrichment, OR, and enrichment score for real target set for each subontology
bp.real <- calc_enrich(target, bg, "BP")
mf.real <- calc_enrich(target, bg, "MF")
cc.real <- calc_enrich(target, bg, "CC")

bp.real.clean <- bp.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, logp, target_frac, target_total, bg_frac, enrich) %>%
  arrange(ID)
mf.real.clean <- mf.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total"), remove = F) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac,
         enrich = log2(target_frac/bg_frac),
         newpval = ifelse(is.na(pvalue), 1, pvalue),
         logp = -log10(newpval)) %>%
  dplyr::select(ID, logp, target_frac, target_total, bg_frac, enrich) %>%
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

# merge observed results with permutation results
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

# calculate empirical pvalues (number of instances where real data is greater than or equal to permutation)
bp.p <- bp.merge %>%
  group_by(ID) %>%
  mutate(gt_target = target_frac.perm <= target_frac.real,
         gt_enrich = enrich.perm <= enrich.real)
bp.cols <- sapply(bp.p[,10:11], as.numeric)
bp.pval <- bp.p %>%
  group_by(ID) %>%
  mutate(sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  select(ID, pVal_target, pVal_enrich) %>%
  distinct()

mf.p <- mf.merge %>%
  group_by(ID) %>%
  mutate(gt_target = target_frac.perm <= target_frac.real,
         gt_enrich = enrich.perm <= enrich.real)
mf.cols <- sapply(mf.p[,10:11], as.numeric)
mf.pval <- mf.p %>%
  group_by(ID) %>%
  mutate(sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  select(ID, pVal_target, pVal_enrich) %>%
  distinct()

cc.p <- cc.merge %>%
  group_by(ID) %>%
  mutate(gt_target = target_frac.perm <= target_frac.real,
         gt_enrich = enrich.perm <= enrich.real)
cc.cols <- sapply(cc.p[,10:11], as.numeric)
cc.pval <- cc.p %>%
  group_by(ID) %>%
  mutate(sum_target = sum(gt_target) + 1,
         pVal_target = sum_target/1001,
         sum_enrich = sum(gt_enrich) + 1,
         pVal_enrich = sum_enrich/1001) %>%
  select(ID, pVal_target, pVal_enrich) %>%
  distinct()

# write out GO terms
write_delim(bp.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/GOterms_BP.tsv", delim = "\t", col_names = T)
write_delim(mf.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/GOterms_MF.tsv", delim = "\t", col_names = T)
write_delim(cc.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/GOterms_CC.tsv", delim = "\t", col_names = T)

# write out sig GO terms
filter(bp.pval, pVal_enrich <= 0.05) %>% dplyr::select(-c(pVal_target)) %>% mutate(subontology = "BP") %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/sigGOterms_BP.tsv", delim = "\t", col_names = T)
filter(mf.pval, pVal_enrich <= 0.05) %>% dplyr::select(-c(pVal_target)) %>% mutate(subontology = "MF") %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/sigGOterms_MF.tsv", delim = "\t", col_names = T)
filter(cc.pval, pVal_enrich <= 0.05) %>% dplyr::select(-c(pVal_target)) %>% mutate(subontology = "CC") %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/sigGOterms_CC.tsv", delim = "\t", col_names = T)

