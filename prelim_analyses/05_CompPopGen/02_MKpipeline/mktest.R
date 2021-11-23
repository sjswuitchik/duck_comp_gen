#!/usr/bin/Rscript

library(tidyverse)
library(lme4)
library(arm)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

cds <- read_tsv("onlyCDS.genes.bed", col_names = c("chr", "start", "end", "gene")) %>%
  mutate(cds.temp = end - start) %>%
  group_by(gene) %>% 
  summarise(cds.len = sum(cds.temp))

callable <- read_tsv("callable.cds.bed", col_names = c("chr", "start", "end", "gene")) %>%
  mutate(call.temp = end - start) %>%
  group_by(gene) %>% 
  summarise(call.len = sum(call.temp))

ingroup <- read_tsv(args[1], col_names = c("chr", "start", "end", "effect", "gene"))
mis.in <- ingroup %>% 
  group_by(gene) %>% 
  tally(effect ==  "missense_variant") %>%
  set_names(c("gene","PR"))
syn.in <- ingroup %>%
  group_by(gene) %>%
  tally(effect == "synonymous_variant") %>%
  set_names(c("gene", "PS"))

outgroup <- read_tsv(args[2], col_names = c("chr", "start", "end", "effect", "gene"))
mis.out <- outgroup %>% 
  group_by(gene) %>% 
  tally(effect ==  "missense_variant") %>%
  set_names(c("gene","FR"))
syn.out <- outgroup %>%
  group_by(gene) %>%
  tally(effect == "synonymous_variant") %>%
  set_names(c("gene", "FS"))

# create full table, replacing any NAs with 0s
MKtable <- full_join(mis.out, syn.out, by = "gene") %>% full_join(mis.in, by = "gene") %>% full_join(syn.in, by = "gene") %>%
  mutate(FR = replace_na(FR, 0)) %>%
  mutate(FS = replace_na(FS, 0)) %>%
  mutate(PR = replace_na(PR, 0)) %>%
  mutate(PS = replace_na(PS, 0)) 

# qc: check for complete NA conversion in each column
MKtable %>% summarise(count = sum(is.na(FR)))
MKtable %>% summarise(count = sum(is.na(FS)))
MKtable %>% summarise(count = sum(is.na(PR)))
MKtable %>% summarise(count = sum(is.na(PS)))

# qc: make sure no ratios > 1
check <- full_join(callable, cds, by = "gene") %>%
  mutate(check = call.len/cds.len) %>%
  mutate(call.len = replace_na(call.len, 0)) %>%
  mutate(cds.len = replace_na(cds.len, 0)) %>%
  mutate(check = replace_na(check, 0))
max(check$check)
check %>% filter(check > 1)

# Prep for SnIPRE
in.miss <- read_delim(args[3], delim = "\t") %>% 
  mutate(indv = N_DATA - N_MISS) %>%
  dplyr::select(indv)
out.miss <- read_delim(args[4], delim = "\t") %>% 
  mutate(indv = N_DATA - N_MISS) %>%
  dplyr::select(indv)

snipre_data <- full_join(MKtable, callable, by = "gene") %>%
  mutate(Tsil = call.len * 0.33) %>%
  mutate(Trepl = call.len * 0.66) %>%
  add_column(npop = median(in.miss$indv)/2) %>%
  add_column(nout = median(out.miss$indv)/2) %>%
  dplyr::select(-c(call.len)) %>%
  filter((Trepl/Tsil)<5) %>%
  filter((PR+FR+PS+FS)>1) 

# functions for MK calculations and DOS 
mk_test <- function(dn=dn,ds=ds,pn=pn,ps=ps){
  dnds_mat <- matrix(data=c(ds,dn,ps,pn),nrow=2,byrow = F)
  pval = fisher.test(dnds_mat)$p.value
  alpha = 1-((ds*pn)/(dn*ps))
  return(list(pval,alpha))
} 

mk_tibble_calc <- function(snipre_res_obj){
  snipre_res_obj_new <- snipre_res_obj %>%
    rowwise %>%
    mutate(mk_pval = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[1]],
           alpha = mk_test(dn=FR,ds=FS,pn=PR,ps=PS)[[2]]) %>%
    ungroup %>%
    dplyr::mutate(mk_pval_fdr = p.adjust(mk_pval,method="BH"))
  return(snipre_res_obj_new)
}

MKtest <- mk_tibble_calc(snipre_data) %>%
  mutate(dos = FR/(FR + FS) - PR/(PR + PS),
         total_poly = PR + PS,
         total_div = FR + FS) 

write.table(MKtest, "MK_output.txt", sep = "\t", quote = F, row.names = F)

# Run SnIPRE 
source("SnIPRE_source.R")
source("my.jags2.R")

snipre.res <- SnIPRE(MKtest)
snipre.qres <- snipre.res$new.dataset
snipre.model <- snipre.res$model
snipre.table <- table(snipre.qres$SnIPRE.class)

write.table(snipre.qres, "snipre_output.txt", sep = "\t", quote = F, row.names = F)
