library(tidyverse)

abs <- read_delim("abs_clean_final.csv", delim = ',') %>% 
  separate("hetAtr_gene", into = c("spp_gene", "fdr_pval"), sep = ':') %>%
  mutate(fdr_pval = as.numeric(fdr_pval))
bus <- read_delim("~/Desktop/busted_output_clean.csv", delim = ',') %>%
  select(file, lrt, pval, omega2 = 'unconstrained omega 2', prop2 = 'proportion 2')
pv <- data.frame(adjP = p.adjust(bus$pval, method = "fdr"), pval = bus$pval)
plot(-log10(pv$adjP), -log10(pv$pval))
pv <- pv %>% select(adjP)
bus_adjP <- bind_cols(bus, pv) %>%
  filter(adjP < 0.05)

df <- inner_join(abs, bus_adjP, by = "file") %>%
  rename(abs_total_branches = total_branches, abs_sig_branches = sig_branches, abs_spp_gene = spp_gene, abs_adjP = fdr_pval, bus_lrt = lrt, bus_pval = pval, bus_omega2 = omega2, bus_prop2 = prop2, bus_adjP = adjP)

write_delim(df, "hyphy_final_ogs.csv", delim = ',')

hist(df$bus_adjP)
cor(df$bus_adjP, df$abs_adjP)
plot(df$bus_adjP, df$abs_adjP)
plot(df$abs_sig_branches, df$abs_total_branches)
plot(df$bus_omega2, df$bus_prop2)
plot(df$bus_prop2, df$bus_adjP)

abs_stripped <- abs_clean %>%
  select(-c(total_branches, sig_branches))
test <- left_join(df, abs_stripped, by = "file")
write_delim(test, "hyphy_final_branches.csv", delim = ',')

ogs <- read_delim("Orthogroups.tsv", delim = '\t') %>%
  rename(file = Orthogroup) %>%
  select(-c('...2'))
sig_ogs <- df %>%
  select(file)
orthos <- left_join(sig_ogs, ogs, by = "file") %>%
  select(file, galGal_protein) %>%
  na.omit()

orthos %>%
  select(galGal_protein) %>%
  write_delim(., "galGal_protIDs.tsv", delim = '\t')

ogs %>%
  select(galGal_protein) %>%
  na.omit() %>%
  write_delim("all_galGal_protIDs.tsv", delim = '\t')
