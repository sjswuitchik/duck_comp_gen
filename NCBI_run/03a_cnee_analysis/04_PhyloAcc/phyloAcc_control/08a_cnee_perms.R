#### new perms using shuf #### 
library(RColorBrewer)
library(ggrepel)
library(tidyverse)

data <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/cnee_perms.counts.bed", delim = "\t", col_names = F)
data.convert <- str_replace_all(data$X4, regex("gene-"), "")
data.clean <- bind_cols(data.convert, data) %>%
  select(-c(X4)) %>%
  rename(X4 = ...1) %>%
  select(X1:X3, X4, X5:X1006)
# calculating pVal this way is off by an order of magnitude .. sum(X7:X1006) isn't working the way I want it to, but count() and n() don't work either and sumRows() will add the actual numbers (e.g. if perm = 2 and is gt obs, then will be summed as 2 perms gt obs instead of one perm)
pVal <- data.clean %>%  
  rowwise() %>%
  mutate(pVal = sum((X7:X1006 >= X5)+1)/(1000+1)) %>%
  select(X1:X6, pVal) %>%
  rename(chr = X1, start= X2, end = X3, gene = X4, accel = X5, total = X6)

# try calculating pVal this way instead - set logical TRUE for count > obs, convert logical to numerical, then sum per gene 
t <- data.clean %>% 
  select(X4, X5, X6, X7:X1006) %>% 
  rename(gene = X4, obs = X5, total = X6) %>% 
  pivot_longer(cols = starts_with("X"), names_to = "perms", values_to = "count") %>%
  mutate(gt = count >= obs)
cols <- sapply(t, is.logical) 
t[,cols] <- lapply(t[,cols], as.numeric)
tpval <- t %>%
  group_by(gene) %>%
  mutate(sum = sum(gt) + 1,
         pVal = sum/1001)
tidyt <- tpval %>%
  pivot_wider(names_from = perms, values_from = count) %>%
  select(gene, obs, total, pVal) %>%
  filter(obs > 0) %>%
  distinct() 

# Adjust p old school
pv <- data.frame(adjP = p.adjust(tidyt$pVal, method = "fdr"), pVal = tidyt$pVal)
# plot to make sure adjustment worked (ie/ not a 1:1 line)
plot(-log10(pv$adjP), -log10(pv$pVal))
# only keep the adjusted p-values
pv <- pv %>% select(-c(pVal))
adjP <- bind_cols(tidyt, pv) %>%
  mutate(sig_class = case_when(
    adjP <= 0.05 ~ "Enriched for accelerated regions (5% FDR)",
    adjP > 0.05 ~ "Not enriched"))


ggplot(adjP, aes(x = obs, y = total, col = sig_class, label = gene)) +
  theme_classic() + 
  scale_y_log10() + 
  geom_jitter(shape = 16) +
  scale_colour_brewer(palette = "Dark2") +
  labs(x = "Number of accelerated CNEEs near gene", y = "Total number of CNEEs near gene", color = "Significance") + 
  scale_x_continuous(breaks = c(1,2,3,4,5)) #+
  geom_label_repel(data = subset(adjP, adjP <= 0.05),
                   aes(label = gene),
                   box.padding = 1.5,
                   point.padding = 0.5,
                   segment.size = 0.2,
                   force = 100,
                   segment.colour = 'grey50')

  
geneList <- adjP %>% filter(adjP < 0.05) %>% select(gene) %>% write_delim(., "~/Desktop/PDF/CNEEs/PhyloAcc_control/geneList.txt", delim = "\t", col_names = T)

# two significant CNEEs with only 1 obs
adjP %>% filter(adjP < 0.05) %>% filter(obs == 1) %>% View()



library(biomaRt) # installed dev version with BiocManager::install('grimbough/biomaRt')
mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
geneList <- adjP %>% filter(adjP < 0.05) %>% dplyr::select(gene)
martList <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = geneList, bmHeader = T, mart = mart)

collapse <- martList %>% 
  dplyr::rename(gene = 'Gene name', goID= `GO term accession`, goTerm = `GO term name`) %>%
  group_by(gene) %>%
  summarise(goID = paste(sort(unique(goID)), collapse = ", "),
            goTerm = paste(sort(unique(goTerm)), collapse = ", "))

# remove genes without GO ids
sub1 <- collapse[!(is.na(collapse$goID) | collapse$goID == ""), ] 
clean.mart <- sub1[-1,]

left_join(geneList, clean.mart, by = "gene") %>% write_delim(., "~/Desktop/PDF/CNEEs/PhyloAcc_control/martList.tsv", delim = "\t", col_names = T)
