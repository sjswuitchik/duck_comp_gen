#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler)

# load in background data and clean up NCBI IDs
bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
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
write_delim(bp.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/bp.perms.clean.tsv", delim = "\t", col_names = T)
write_delim(mf.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/mf.perms.clean.tsv", delim = "\t", col_names = T)
write_delim(cc.perms.clean, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/cc.perms.clean.tsv", delim = "\t", col_names = T)

# read in permutation results if working in a new session
bp.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/bp.perms.clean.tsv", col_names = T)
mf.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/mf.perms.clean.tsv", col_names = T)
cc.perms.clean <- read_tsv("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/cc.perms.clean.tsv", col_names = T)

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
write_delim(bp.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/GOterms_BP.tsv", delim = "\t", col_names = T)
write_delim(mf.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/GOterms_MF.tsv", delim = "\t", col_names = T)
write_delim(cc.pval, "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/GOterms_CC.tsv", delim = "\t", col_names = T)

# write out sig GO terms
filter(bp.pval, pVal_enrich <= 0.05) %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/sigGOterms_BP.tsv", delim = "\t", col_names = T)
filter(mf.pval, pVal_enrich <= 0.05) %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/sigGOterms_MF.tsv", delim = "\t", col_names = T)
filter(cc.pval, pVal_enrich <= 0.05) %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/sigGOterms_CC.tsv", delim = "\t", col_names = T)

#### QC #### 
# randomly sample IDs to make sure distributions of target_frac and and enrich look okay i.e. normal distribution (just with BP for starters)
ids <- bp.perms.clean %>% group_by(ID) %>% distinct(ID)
sample(1:1587, 15, replace=FALSE) 
# 757 867 164 1138 981 662 133 231 846 95 892 1456 1458 694 1147
ids[1147,]
# GO:0034220 GO:0042886 GO:0006281 GO:0048762 GO:0045216 GO:0031396 GO:0003013 GO:0006939 GO:0042255 GO:0002429 GO:0043270 GO:0071383 GO:0071407 GO:0032332 GO:0048869
test <- bp.perms.clean %>% group_by(ID) %>% filter(ID == 'GO:0048869')
ggplot(data = test, aes(x = target_frac)) + geom_histogram()
ggplot(data = test, aes(x = enrich)) + geom_histogram()

# make sure real data is well represented by the permutations (i.e. does the real data fall within permutation or is it on the tails?)
test <- bp.perms.clean %>% 
  separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
  mutate(target_total = as.numeric(target_total)) %>%
  group_by(ID)
filt <- test %>% filter(ID == 'GO:0048869')
ggplot(data = filt, aes(x = target_total)) + 
  geom_histogram() + 
  geom_text(x = 85, y = 0, colour = "red", size = 8, label = '*')

## now for MF, just to make sure it's consistent across subontologies
ids <- mf.perms.clean %>% group_by(ID) %>% distinct(ID)
sample(1:528, 15, replace=FALSE) 
# 521 207 321  29  42 357  60 481 204 446 271 155 169 448  18
ids[18,]
# GO:0005549 GO:0043177 GO:0030594 GO:0003727 GO:0004675 GO:0098531 GO:0005126 GO:0034212 GO:0043167 GO:0034979 GO:0004721 GO:0019212 GO:0022804 GO:0042165 GO:0003674
test <- mf.perms.clean %>% group_by(ID) %>% filter(ID == 'GO:0003674')
ggplot(data = test, aes(x = target_frac)) + geom_histogram()
ggplot(data = test, aes(x = enrich)) + geom_histogram()

# make sure real data is well represented by the permutations (i.e. does the real data fall within permutation or is it on the tails?)
test <- mf.perms.clean %>% 
  separate(GeneRatio, into = c("target_in", "target_total"), sep = '/', remove= F) %>%
  mutate(target_total = as.numeric(target_total)) %>%
  group_by(ID)
filt <- test %>% filter(ID == 'GO:0003674')
ggplot(data = filt, aes(x = target_total)) + 
  geom_histogram() + 
  geom_text(x = 82, y = 0, colour = "red", size = 8, label = '*')
