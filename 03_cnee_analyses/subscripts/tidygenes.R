#!/usr/bin/Rscript

library(tidyverse)

args <- commandArgs(trailingOnly = T)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).\n", call.=FALSE)
}

data <- read_delim(args[1], delim = "\t", col_names = c("chr", "start", "end", "gene"))

t <- data %>%
  group_by(gene) %>%
  summarise(chrom = chr,
            start = min(start),
            end = max(end)) %>%
  distinct() %>%
  select(chrom, start, end, gene) %>%
  write_delim(., "galGal.tidygenes.bed", delim = "\t", col_names = F)
