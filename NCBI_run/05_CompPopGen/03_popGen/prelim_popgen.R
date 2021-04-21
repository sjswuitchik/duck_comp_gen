library(tidyverse)
library(ggrepel)
library(tidygenomics)

# load in pi output from vcftools - chrom, start, end, nvar, pi, and add window size
pi <- read.table("hetAtr.pi.bial.windowed.pi", header = T) %>%
  rename(chrom = CHROM, start = BIN_START, end = BIN_END, nvar = N_VARIANTS, pi = PI) %>%
  mutate(window = end - start + 1)
# hist of initial data
hist(pi$pi, br = 20)
# load in clean callable sites - chrom, start, end, and add interval length
call.sites <- read_delim("hetAtr_coverage_sites_clean_merged.bed", delim = '\t', col_names = F) %>%
  rename(chrom = X1, start = X2, end = X3) %>%
  mutate(int_len = end - start) 
# intersect pi and callable (akin to bedtools intersect but in R)
combo <- genome_intersect(pi, call.sites, by=c("chrom", "start", "end"), mode = "both") 
# calculate pi for each window by accounting for callable sites 
pibywindow <- combo %>%
  select(-c(start, end, nvar)) %>% # get rid of variables not needed for this 
  group_by(chrom, pi) %>% # group first by chromosome, then by value of pi to keep windows consistent
  #summarise(call = sum(int_len), .groups = 'keep',) # sum the length of all the callable intervals in each window
  mutate(call = sum(int_len)) %>% # add number of callable sites per window
  select(-c(int_len)) %>% # get rid of interval length now, no longer needed 
  group_by(chrom, pi, call) %>% # group by chrom, then pi, then number of call sites to keep windows consistent
  mutate(call.pi = ((pi*window)/call))  %>% # standardize pi by callable sites within a window
  group_by(chrom) %>% # group by chrom
  mutate(meanPi = mean(call.pi)) %>% # calculate mean pi per chromosome
  unique() %>% # keep unique entries only (removes duplicated entries that result from ungrouping)
  summarise(meanPi = mean(call.pi), .groups = 'keep') # summarise mean pi by chromosome/scaffold
hist(pibywindow$meanPi, br = 200, xlim = c(0,0.1)) # check dist 
mean(pibywindow$meanPi) # check mean - implies Ne 10k-100k, depending on mutation rate

plot <- pibywindow %>%
  mutate(piplot = ifelse(meanPi > 0.05, 0.05, meanPi))
hist(plot$piplot, br = 20) # much better, now make nicer

ggplot(plot, aes(x = piplot)) +
  geom_histogram(bins = 100) + 
  theme_classic() + 
  labs(y = "Frequency", x = "Nucleotide diversity (pi)")
mean(plot$piplot) # 0.003931766

# plot Tajima's D quickly
taj <- read.table("hetAtr.bial.100kb.Tajima.D", header = T)
hist(taj$TajimaD, br = 20)

# plot PCA, but nicer
pca1 <- read.table("hetAtr.eigenvec", sep = " ", header = F)
plot(data=pca1, V3~V4)

ggplot(pca1, aes(x = V3, y = V4)) +
  geom_point() + 
  theme_classic() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  labs(x = "PC1", y = "PC2") + 
  geom_label_repel(aes(label = ifelse(V4 > 0.75, as.character(V2), ''))) + 
  geom_label_repel(aes(label = ifelse(V3 > 0.75, as.character(V2), '')))

# plot IBC quickly
ibc <- read.table("hetAtr.ibc",header=T)
plot(ibc$Fhat1,ibc$Fhat2,xlab="Fhat1",ylab="Fhat2")

# load in missingness data
miss <- read.table("missing_data_per_ind.txt", header = T) %>% 
  mutate(miss = N_MISS/N_DATA) %>%
  rename(sample = INDV) %>%
  select(sample, miss)
# load in bam stats
cov <- read.table("bam_sumstats.txt", header = T) %>%
  rename(meanSD = MeanSeqDepth) %>%
  select(sample, meanSD)
# combine
data <- full_join(miss, cov, "sample" = "sample") %>%
  filter(sample != "hetAtr_female",
         sample != "stiNae_male") %>%
  as_tibble()

ggplot(data, aes(x = miss, y = meanSD)) +
  geom_point() + 
  theme_classic() + 
  theme(plot.margin = unit(c(1,1,1,1), "cm")) + 
  labs(x = "Relative Missingness", y = "Mean Sequence Depth") + 
  geom_label_repel(aes(label = ifelse(miss > 0.003, as.character(sample), ''))) + 
  geom_label_repel(aes(label = ifelse(meanSD > 30, as.character(sample), ''))) 
