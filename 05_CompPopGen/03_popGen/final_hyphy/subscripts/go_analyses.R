library(tidyverse)

# load in Entrez IDs from NCBI batch retrieve
target <- read_delim("entrezIDs.tsv", delim = '\t', col_names = c("prot", "gene_sym", "gene_id")) %>%
  separate(gene_id, into = c(NA, "ncbi"), sep = "GeneID:", remove = T)
bg <- read_delim("bg_entrezIDs.tsv", delim = '\t', col_names = c("prot", "gene_sym", "gene_id")) %>%
  separate(gene_id, into = c(NA, "ncbi"), sep = "GeneID:", remove = T)

calc_enrich <- function(targetset, background, ont) { 
  enrichGO(targetset$ncbi,'org.Gg.eg.db',
           pvalueCutoff=1.5,
           qvalueCutoff = 1.5, minGSSize=5, maxGSSize=2000, 
           pAdjustMethod="BH",
           universe=background$ncbi,
           keyType="ENTREZID",
           ont=ont) 
}

bp <- calc_enrich(target, bg, ont = "BP")
mf <- calc_enrich(target, bg, ont = "MF")
cc <- calc_enrich(target, bg, ont = "CC")

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
cc.clean <- cc@result %>% 
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

bp.id <- bp.clean %>%
  select(ID)
mf.id <- mf.clean %>%
  select(ID)
cc.id <- cc.clean %>%
  select(ID)
id.df <- bind_rows(bp.id, mf.id, cc.id)
goList <- remove_rownames(id.df)

mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
martList <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = goList, bmHeader = T, mart = mart)

#### pull name, not gene 
collapse <- martList %>% 
  dplyr::rename(goID= `GO term accession`, goTerm = `GO term name`) %>%
  group_by(goID) %>%
  summarise(goTerm = paste(sort(unique(goTerm)), collapse = ", "))

# remove genes without GO ids
sub1 <- collapse[!(is.na(collapse$goID) | collapse$goID == ""), ] 
clean.mart <- sub1[-1,]

left_join(goList, clean.mart, by = c("ID" = "goID")) %>% 
  na.omit() %>% 
  write_delim(., "martList.tsv", delim = "\t", col_names = T)
