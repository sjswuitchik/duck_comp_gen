setwd("~/Desktop/PDF/duck_assemblies")
library(tidyverse)

celen <- read_tsv("ce.lengths", col_names = c("length", "set"))
allcelen <- read_tsv("galGal6_allce.lengths", col_names = c("length"))

celen %>% group_by(set) %>% summarize(min = min(length), max=max(length), mean=mean(length), median=median(length))

celen %>% filter(set == "Craig", length > 1) %>% mutate(length = pmin(length, 1000)) %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth=1)

celen %>% filter(set == "Sackton", length > 1) %>% mutate(length = pmin(length, 1000)) %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth=1)

celen %>% filter(set == "UCSC", length > 1) %>% mutate(length = pmin(length, 1000)) %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth=1)

celen %>% filter(set == "Lowe", length > 1) %>% mutate(length = pmin(length, 1000)) %>% ggplot(aes(x=length)) +
  geom_histogram(binwidth=1)

celen %>% filter(length > 1) %>% mutate(length = pmin(length, 300)) %>% 
  ggplot(aes(x=length, group=set, color=set)) + geom_density()
ggsave("CNEE_comp.pdf")

celen %>% filter(length > 1000) %>% group_by(set) %>% summarize(count = n())

allcelen %>% mutate(length = pmin(length,300)) %>% ggplot(aes(x=length)) + geom_density()
