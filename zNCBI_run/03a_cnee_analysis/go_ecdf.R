library(tidyverse) 

wd <- "/n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/03a_cnee_analysis/postPhyloAcc/go_perms"

for (whichont in c("BP", "MF")) {
  spec_patt <- glob2rx(paste0("GO_top1_*_", whichont, "perm.tsv"))
  files <- list.files(path = wd, pattern = spec_patt, full.names = T)
  results <- list()
  for (file in files) {
    results[[file]] <- read_tsv(file) %>% 
      filter(!is.na(target_frac), target_frac < 1) %>% 
      mutate(enrich = log2(target_frac/bg_frac))
  }
  bind_rows(results) %>% 
    group_by(set, ID) %>% 
    summarize(ecdf_frac = list(ecdf(target_frac)),
              ecdf_enrich = list(ecdf(enrich)),
              ecdf_logp = list(ecdf(logp.perm))) %>%
    saveRDS(file = paste0(wdir, "GO_top1_", whichont, ".robj"))
}


for (whichont in c("BP", "MF")) {
  spec_patt <- glob2rx(paste0("GO_top2_*_", whichont, "perm.tsv"))
  files <- list.files(path = wd, pattern = spec_patt, full.names = T)
  results <- list()
  for (file in files) {
    results[[file]] <- read_tsv(file) %>% 
      filter(!is.na(target_frac), target_frac < 1) %>% 
      mutate(enrich = log2(target_frac/bg_frac))
  }
  bind_rows(results) %>% 
    group_by(set, ID) %>% 
    summarize(ecdf_frac = list(ecdf(target_frac)),
              ecdf_enrich = list(ecdf(enrich)),
              ecdf_logp = list(ecdf(logp.perm))) %>%
    saveRDS(file = paste0(wdir, "GO_top2_", whichont, ".robj"))
}


for (whichont in c("BP", "MF")) {
  spec_patt <- glob2rx(paste0("GO_top3_*_", whichont, "perm.tsv"))
  files <- list.files(path = wd, pattern = spec_patt, full.names = T)
  results <- list()
  for (file in files) {
    results[[file]] <- read_tsv(file) %>% 
      filter(!is.na(target_frac), target_frac < 1) %>% 
      mutate(enrich = log2(target_frac/bg_frac))
  }
  bind_rows(results) %>% 
    group_by(set, ID) %>% 
    summarize(ecdf_frac = list(ecdf(target_frac)),
              ecdf_enrich = list(ecdf(enrich)),
              ecdf_logp = list(ecdf(logp.perm))) %>%
    saveRDS(file = paste0(wdir, "GO_top3_", whichont, ".robj"))
}
