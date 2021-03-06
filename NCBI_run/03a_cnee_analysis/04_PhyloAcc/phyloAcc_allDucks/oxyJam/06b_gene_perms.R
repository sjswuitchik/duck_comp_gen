#### new perms using shuf #### 
library(RColorBrewer)
library(ggrepel)
library(tidyverse)

data <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/oxyJam/cnee_perms.counts.bed", delim = "\t", col_names = F)
data.convert <- str_replace_all(data$X4, regex("gene-"), "")
data.clean <- bind_cols(data.convert, data) %>%
  dplyr::select(-c(X4)) %>%
  rename(X4 = ...1) %>%
  select(X1:X3, X4, X5:X1006)

# calculating pVals - set logical TRUE for count > obs, convert logical to numerical, then sum per gene 
t <- data.clean %>% 
  dplyr::select(X4, X5, X6, X7:X1006) %>% 
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
# plot to make sure adjustment worked (ie/ not an exact 1:1 line)
plot(-log10(pv$adjP), -log10(pv$pVal))
# only keep the adjusted p-values
pv <- pv %>% select(-c(pVal))
adjP <- bind_cols(tidyt, pv) %>%
  mutate(sig_class = case_when(
    adjP <= 0.05 ~ "Enriched for accelerated regions (5% FDR)",
    adjP > 0.05 ~ "Not enriched")) 

# a few organization checks  
adjP %>% arrange(desc(obs)) %>% View()
adjP %>% arrange(gene) %>% View()
adjP %>% filter(adjP <= 0.05) %>% count() # 0 genes

# add in 0 data to make plotting x axis consistent between spp
zero <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/oxyJam/dummy_data_cnee_acc.csv", delim = ',', col_names = T)
plotdf <- bind_rows(adjP, zero)

ggplot(plotdf, aes(x = obs, y = total, col = sig_class, label = gene, alpha = obs > 2)) +
  theme_classic() + 
  scale_y_log10() + 
  geom_jitter(shape = 16, width = 0.2, size = 2, show.legend = F, colour = c("#D95F02")) +
  scale_colour_brewer(palette = "Dark2") +
  scale_alpha_discrete(range = c(0.5, 0)) + 
  labs(x = "Number of accelerated CNEEs near gene", y = "Total number of CNEEs near gene", color = "Significance") + 
  scale_x_continuous(breaks = c(1,2,3,4,5,6,7)) 

geneList <- adjP %>% filter(adjP <= 0.05) %>% select(gene) %>% write_delim(., "~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/oxyJam/geneList.txt", delim = "\t", col_names = T)
