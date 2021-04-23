library(tidyverse)

#read in data, rename col names
data <- read_delim("hetAtr.frq.count", delim = '\t', col_names = T) %>%
  rename(chrom = CHROM, pos = POS, n_alleles = N_ALLELES, n_chr = N_CHR, major = MAJOR, AC = MINOR) 

# calc folded SFS
fold <- data %>%
  mutate(fold = ifelse(AC > 15, 30 - AC, AC)) %>%
  filter(fold > 0)
# put SFS into its own var
sfs <- table(fold$fold)
# check SFS
pdf("sfs.pdf")
plot(sfs)
dev.off()
# write out so you can manually put into the Stairway blueprint file 
sfs %>% 
  as.data.frame() %>%
  write_delim(., "sfs.tsv", delim = '\t', col_names = T)
