#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
} 

read_delim(args[1], col_names = F) %>%
  select(X9) %>%
  separate(X9, into = c("trans", NA), sep = ';') %>%
  separate(trans, into = c(NA, "id"), sep = ' ') %>%
  write_delim(args[2], col_names = F, quote = NULL)
