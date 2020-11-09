#### run GO analysis ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(clusterProfiler)
bg <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/cnee_goperms.counts.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, gene = X4, combo = X5, accel = X6, total = X7) %>%
  separate(combo, into = c(NA, "ncbi"), sep = "GeneID:", remove = T) %>%
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
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac) %>%
    dplyr::select(ID, target_frac, bg_frac, OR) %>%
    arrange(ID) 
  mf.clean <- mf@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac) %>%
    dplyr::select(ID, target_frac, bg_frac, OR) %>%
    arrange(ID) 
  cc.clean <- cc@result %>% 
    separate(GeneRatio, into = c("target_in", "target_total")) %>%
    separate(BgRatio, into = c("bg_in", "bg_total")) %>%
    mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
           bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
           OR = target_frac/bg_frac) %>%
    dplyr::select(ID, target_frac, bg_frac, OR) %>%
    arrange(ID) 
  bp.perms <- as.data.frame(rbind(bp.perms, bp.clean))
  mf.perms <- as.data.frame(rbind(mf.perms, mf.clean))
  cc.perms <- as.data.frame(rbind(cc.perms, cc.clean))
  
}

bp.real <- calc_enrich(target, bg, "BP")
mf.real <- calc_enrich(target, bg, "MF")
cc.real <- calc_enrich(target, bg, "CC")

bp.real.clean <- bp.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac) %>%
  dplyr::select(ID, target_frac, bg_frac, OR) %>%
  arrange(ID)

mf.real.clean <- mf.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac) %>%
  dplyr::select(ID, target_frac, bg_frac, OR) %>%
  arrange(ID)

cc.real.clean <- cc.real@result %>% 
  separate(GeneRatio, into = c("target_in", "target_total")) %>%
  separate(BgRatio, into = c("bg_in", "bg_total")) %>%
  mutate(target_frac = as.numeric(target_in)/as.numeric(target_total),
         bg_frac = as.numeric(bg_in)/as.numeric(bg_total),
         OR = target_frac/bg_frac) %>%
  dplyr::select(ID, target_frac, bg_frac, OR) %>%
  arrange(ID)
  
left_join(bp.real.clean, bp.perms, by = "ID", suffix = c(".real", ".perm"))
