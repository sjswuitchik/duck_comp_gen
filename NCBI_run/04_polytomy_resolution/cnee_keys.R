#!/bin/R

data <- read_delim("galGal6_final_merged_CNEEs_named.bed", delim = '\t', col_names = F) %>%
rename(chr = X1, start = X2, end = X3, cnee = X4) %>%
mutate(int_len = end - start) %>%
filter(int_len > 500) %>%
write_delim(., "cnees_big.bed", delim = '\t', col_names = F)

cnee <- data %>% select(cnee)

tab <- read_delim("all_cnees.tab", delim = '\t', col_names = F) %>%
rename(spp = X1, cnee = X2, seq = X3)

combo <- inner_join(data, tab, by = "cnee")

length(unique(combo$cnee)) # there are 7499 unique CNEEs 

combo %>% group_by(cnee) %>% arrange(cnee) %>% tally() %>% write_delim(., "spp_by_cnee.tsv", delim = '\t', col_names = T) # count number of spp that has each CNEE

ducks <- filter(combo, spp == "hetAtr" | spp == "netAur" | spp == "oxyJam" | spp == "stiNae")

length(unique(ducks$cnee)) # checks number of unique CNEEs (7388)
length(unique(ducks$spp)) # checks that there are only 4 focal spp (4)

ducks %>% group_by(cnee) %>% arrange(cnee) %>% tally() %>% write_delim(., "spp_by_cnee_ducksOnly.tsv", delim = '\t', col_names = T)

combo %>% group_by(cnee) %>% arrange(cnee) %>% tally() %>% filter(n == 15) %>% write_delim(., "spp_by_cnee_allSpp.tsv", delim = '\t', col_names = T)
