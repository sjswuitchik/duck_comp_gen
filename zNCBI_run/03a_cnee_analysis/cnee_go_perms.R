#### GO enrichment perms ####
library(org.Gg.eg.db) # BiocManager::install("org.Gg.eg.db")
library(tidyverse)
library(clusterProfiler) #  BiocManager::install("clusterProfiler")
library(parallel)
library(rlist)

compute_go_results <- function(DF, outname, CORES, PERMS) {
  
  ##INTERNAL FUNCTIONS##
  calc_enrich <- function(targetset, background, ont) { 
    enrichGO(targetset$ncbi,'org.Gg.eg.db',
             pvalueCutoff=1.5,
             qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
             pAdjustMethod="none",
             universe=background,
             keyType="ENTREZID",
             ont=ont) 
  }
  
  get_go_perm <- function(DF, samples, golist, ont) {
    rand <- DF %>% 
      sample_n(samples) %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    background <- DF %>% 
      filter(gene != ".") %>% 
      dplyr::select(gene) %>% 
      separate(gene, into=c("ncbi", "sym"), sep=":") %>% 
      distinct(ncbi)
    
    rand_go <- calc_enrich(targetset=rand, background=background$ncbi, ont=ont)
    
    golist %>% 
      left_join(rand_go@result, by=c("ID" = "ID")) %>% 
      separate(GeneRatio, into=c("target_in", "target_total")) %>% 
      separate(BgRatio, into=c("bg_in", "bg_total")) %>%
      mutate(newpval = ifelse(is.na(pvalue), 1, pvalue), 
             logp.perm = -log10(newpval),
             target_frac = as.numeric(target_in)/as.numeric(target_total), 
             bg_frac = as.numeric(bg_in)/as.numeric(bg_total)) %>%
      dplyr::select(ID, logp.perm, target_frac, bg_frac) %>% 
      arrange(ID)
  }
  
  #to do one permutation of the full set
  get_one_perm_set <- function(perm, input, DF, golist, ont) {
    try(lapply(input, get_go_perm, DF=DF, golist=golist, ont=ont) %>%
          dplyr::bind_rows(.id="    set"), TRUE)
  }
  
}


### REAL WORK ###

args <- commandArgs(trailingOnly = TRUE)
#args are 1 permutation index ID for slurm batch processing, 2 number of cores, 3 number of permutations, 4 is data path

path_to_data <- args[4]

gene_gg <- read_tsv(paste0(path_to_data, "/galGal_background.txt"), col_names = c("cnee", "gene"))

#top1
cnee_top1 <- read_tsv(paste0(path_to_data, "/cnee_final_top1.tsv")) %>% 
  dplyr::select(ID, logBF1, logBF2) %>%
  dplyr::rename(cnee = ID) %>%
  full_join(gene_gg, by = c("cnee" = "cnee")) %>%
  mutate(accel = ifelse(logBF1 >= 10 & logBF2 >= 1), TRUE, FALSE) %>%
  distinct(ID, .keep_all=TRUE) %>%
  dplyr::select(ID, accel, gene)

#top2
cnee_top2 <- read_tsv(paste0(path_to_data, "/cnee_final_top2.tsv")) %>% 
  dplyr::select(ID, logBF1, logBF2) %>%
  dplyr::rename(cnee = ID) %>%
  full_join(gene_gg, by = c("cnee" = "cnee")) %>%
  mutate(accel = ifelse(logBF1 >= 10 & logBF2 >= 1), TRUE, FALSE) %>%
  distinct(ID, .keep_all=TRUE) %>%
  dplyr::select(ID, accel, gene)

#top3
cnee_top3 <- read_tsv(paste0(path_to_data, "/cnee_final_top3.tsv")) %>% 
  dplyr::select(ID, logBF1, logBF2) %>%
  dplyr::rename(cnee = ID) %>%
  full_join(gene_gg, by = c("cnee" = "cnee")) %>%
  mutate(accel = ifelse(logBF1 >= 10 & logBF2 >= 1), TRUE, FALSE) %>%
  distinct(ID, .keep_all=TRUE) %>%
  dplyr::select(ID, accel, gene)


compute_go_results(cnee_top1, paste0(path_to_data, "/original_GO_top1_", "_run", args[1]), args[2], args[3])
compute_go_results(cnee_top2, paste0(path_to_data, "/original_GO_top2_", "_run", args[1]), args[2], args[3])
compute_go_results(cnee_top3, paste0(path_to_data, "/original_GO_top3_", "_run", args[1]), args[2], args[3])

