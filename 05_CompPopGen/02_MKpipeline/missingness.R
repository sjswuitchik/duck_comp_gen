#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

# read in data and calculate missingness
ingroup <- read.delim(args[1]) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(ingroup$missing)
remove <- ingroup %>% 
  filter(missing >= threshold)
write.csv(remove %>% select(INDV), "ingroup.remove.indv", row.names = F, quote = F)

outgroup <- read.delim(args[2]) %>%
  group_by(INDV) %>%
  summarise_each(funs(sum)) %>%
  mutate(missing = N_MISS/N_DATA)
threshold <- 3*median(outgroup$missing)
remove <- outgroup %>% 
  filter(missing >= threshold)
write.csv(remove %>% select(INDV), "outgroup.remove.indv", row.names = F, quote = F)

# NB if only one individual in outgroup, it will be output as an individual to remove 

