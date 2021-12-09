library(tidyverse)

single_orthos <- read_delim("Orthogroups_SingleCopyOrthologues.txt", delim = '\t', col_names = c("Orthogroup"))
all_orthos <- read_delim("Orthogroups.tsv", delim = '\t')

clean_orthos <- left_join(single_orthos, all_orthos, by = "Orthogroup")

clean_orthos %>%
  select(Orthogroup, anaPla = anaPla_protein) %>%
  write_delim(., "anaPla_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, ansBra = ansBra.translated) %>%
  write_delim(., "ansBra_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, ansCyg = ansCyg_protein) %>%
  write_delim(., "ansCyg_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, ansInd = ansInd.translated) %>%
  write_delim(., "ansInd_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, braCan = braCan.translated) %>%
  write_delim(., "braCan_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, colVir = colVir_protein) %>%
  write_delim(., "colVir_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, cotJap = cotJap_protein) %>%
  write_delim(., "cotJap_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, galGal = galGal_protein) %>%
  write_delim(., "galGal_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, hetAtr = hetAtr.translated) %>%
  write_delim(., "hetAtr_ogs.tsv", delim = '\t')

clean_orthos %>%
  select(Orthogroup, netAur = netAur.translated) %>%
  write_delim(., "netAur_ogs.tsv", delim = '\t')

 clean_orthos %>%
  select(Orthogroup, numMel = numMel_protein) %>%
  write_delim(., "numMel_ogs.tsv", delim = '\t')
 
 clean_orthos %>%
  select(Orthogroup, oxyJam = oxyJam.translated) %>%
  write_delim(., "oxyJam_ogs.tsv", delim = '\t')
 
clean_orthos %>%
  select(Orthogroup, stiNae = stiNae.translated) %>%
  write_delim(., "stiNae_ogs.tsv", delim = '\t') 
 
clean_orthos %>%
  select(Orthogroup, syrMik = syrMik.translated) %>%
  write_delim(., "syrMik_ogs.tsv", delim = '\t') 
 
clean_orthos %>%
  select(Orthogroup, tymCupPin = tymCupPin.translated) %>%
  write_delim(., "tymCupPin_ogs.tsv", delim = '\t')
 

write_delim(clean_orthos, "clean_single_orthos.tsv", delim = '\t')
