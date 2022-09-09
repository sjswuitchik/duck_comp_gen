library(tidyverse)

counts <- read_delim("run_ortho/results/Results_Nov09/Orthogroups/Orthogroups.GeneCount.tsv", delim = '\t')
orthos <- read_delim("run_ortho/results/Results_Nov09/Orthogroups/Orthogroups.tsv", delim = '\t')

clean <- counts %>%
  filter(Total <= 25) %>%
  filter(hetAtr.translated != 0) %>%
  rowwise() %>%
  mutate(missing_spp = sum(c_across(anaPla_protein:tymCupPin.translated) == 0)) %>%
  filter(missing_spp <= 7) %>%
  select(Orthogroup)
write_delim(clean, "clean_ogs_geneCount.tsv", delim = '\t')

clean_orthos <- left_join(clean, orthos, by = "Orthogroup")

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
 

write_delim(clean_orthos, "clean_orthos.tsv", delim = '\t')
