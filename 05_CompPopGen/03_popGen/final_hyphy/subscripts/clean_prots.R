#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

spp <- read_delim(args[1], delim = '\t', col_names = F) %>%
  separate(X1, into = c(NA, "X1"), sep = '=') %>%
  write_delim(., args[2], delim = '\t', col_names = F)
