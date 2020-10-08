## uploading for Oct 8 meeting

library(tidyverse)
library(qqman)
setwd("~/Desktop/PDF/duck_assemblies/PhyloAcc_out/NCBI_run")

# load in accelerated CNEEs by window
acc.cnees.window <- read_delim("window.acc.cnees.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, acc = X4)
# convert remaining CNEE names to a binary measure of occurrence (1)
acc.convert.window <- str_replace_all(acc.cnees.window$acc, regex("zfCNEE[0-9]{0,6}"), "1")
# summarise by window
acc.clean.window <- bind_cols(acc.cnees.window, acc.convert.window) %>%
  select(-c(acc)) %>%
  rename(acc = ...5) %>%
  mutate(cnee = as.numeric(acc)) %>%
  group_by(chr, start, end) %>%
  summarise(acc.cnee = sum(cnee))

# load in total list of CNEEs by window
all.cnees.window <- read_delim("window.cnees.bed", delim = "\t", col_names = F) %>%
  rename(chr = X1, start = X2, end = X3, cnee = X4)
# convert CNEE names to a binary measure of occurrence (1)
all.convert.window <- str_replace_all(all.cnees.window$cnee, regex("CNEE[0-9]{0,6}"), "1")
# summarise total number of CNEEs by window
all.clean.window <- bind_cols(all.cnees.window, all.convert.window) %>% 
  select(-c(cnee)) %>% 
  rename(cnee = ...5) %>%
  mutate(cnee = as.numeric(cnee)) %>%
  group_by(chr, start, end) %>%
  summarise(total.cnee = sum(cnee))

# combine, clean, and calculate proportion of accel'd CNEEs per window
data.window <- full_join(all.clean.window, acc.clean.window, by = c("chr", "start", "end")) %>%
  filter(total.cnee > 0)
final.data.window <- data.window %>%
  mutate(prop = acc.cnee/total.cnee) 

# binomial test function
bt <- function(x, n, p = 294/375591) {
  binom.test(x, n, 294/375591, alternative = "greater", conf.level = 0.95)$p.value
}

# add binomial p-values to table
final.data.window$pVal <- mapply(bt, final.data.window$acc.cnee, final.data.window$total.cnee)

# adjust p-values for multiple comparisons
# don't use mutate for this, do it old school
pv <- data.frame(adjustP = p.adjust(final.data.window$pVal, method = "fdr"), pVal = final.data.window$pVal)
# plot to make sure adjustment worked (ie/ not a 1:1 line)
plot(-log10(pv$adjustP), -log10(pv$pVal))
# only keep the adjusted p-values
pv <- pv %>% select(-c(pVal))

adjP <- bind_cols(final.data.window, pv)

adjP %>% filter(adjustP < 0.05) %>%
  write_delim("final.cnees.window.sig.txt", delim = "\t", col_names = T)

test <- left_join(acc.cnees.window, adjP, by = c("chr" = "chr", "start" = "start", "end" = "end")) %>%
  na.omit()
write_delim(test, "~/Desktop/test.txt", delim = "\t", col_names = F)

## bash repl chr:
## wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/315/GCF_000002315.6_GRCg6a/GCF_000002315.6_GRCg6a_assembly_report.txt
## sed 's/\r$//g' GCF_000002315.6_GRCg6a_assembly_report.txt | grep -v "^#" | cut -f1,5,7,10 > galGal6_chr_key
## awk '{print $3, $1}' galGal6_chr_key > acckey
## ./replace_chrs.pl acckey test.txt > test.rep.txt

# also replaced W with 34, Z with 35

testrep <- read_delim("~/Desktop/test.rep.txt", delim = "\t", col_names = F, col_types = "dddcddddd") %>%
  rename(chr = X1, start = X2, end = X3, acc = X4, total = X5, acc.cnee = X6, prop = X7, pVal = X8, pAdj = X9) %>%
  select(-c(acc.cnee)) %>%
  na.omit()

qq(testrep$pAdj)

manhattan(testrep, chr="chr", bp = "start", snp = "acc", p = "pAdj", col = c("grey", "skyblue"), annotatePval = 0.05) # note: W = 34, Z = 35


