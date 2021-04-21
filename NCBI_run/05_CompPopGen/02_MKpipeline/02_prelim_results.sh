library(tidyverse)
library(ggridges)

snip_data <- read_delim("snipre_output.tsv", delim = '\t', col_names = T)

# gamma 
gam <- snip_data %>%
  rename(gamma= SnIPRE.gamma) %>%
  filter(FS != 0) %>%
  filter(FR != 0) %>%
  filter(PR > 0 & PS > 0) %>%
  filter(gamma >= -2 & gamma <= 2) %>%
  select(gamma)

ggplot(gam, aes(x = gamma)) + 
  geom_density(aes(), show.legend = F) + 
  theme_classic() + 
  labs(y = "", x = "Gamma") 

# direction of selection
dos <- snip_data %>%
  dplyr::select(gene, dos)

ggplot(dos, aes(x = dos)) + 
  geom_density(aes(), show.legend = F) + 
  theme_classic() + 
  labs(y = "", x = "Direction of Selection")

# alpha
alpha <- snip_data %>%
  filter(FS != 0) %>%
  filter(FR != 0) %>%
  filter(PR > 0 & PS > 0) %>%
  dplyr::select(gene,alpha) %>%
  na.omit() %>%
  filter(alpha != '-Inf') %>%
  filter(alpha >= -2 & 2 >= alpha)
ggplot(alpha, aes(x = alpha)) + 
  geom_density(aes(), show.legend =F) + 
  theme_classic() + 
  labs(y = "", x = "Alpha")
min(alpha$alpha)

# cumulative alpha
allFR <- sum(snip_data$FR)
allFS <- sum(snip_data$FS)
allPR <- sum(snip_data$PR)
allPS <- sum(snip_data$PS)
adjAlpha = 1 - ((allFS*allPR)/(allFR*allPS)) #hint at little evidence for positive selection and more evidence perhaps for small pop sizes and segregating deleterious mutations

# summary of MK data
mk_summary <- snip_data %>%
  summarize(n_sp = mean(npop),
            n_out=mean(nout),
            n_loci_tested = n(),
            n_sig_mk = sum(mk_pval<0.05),
            n_sig_mk_fdr = sum(mk_pval_fdr<0.1),
            n_pos_snipre_ml = sum(SnIPRE.class == "pos"),
            n_neg_snipre_ml = sum(SnIPRE.class == "neg"), 
            mean_SnIPRE.est = mean(SnIPRE.est))

# sig without test correction
sig <- snip_data %>%
  filter(mk_pval < 0.05) # 715 genes 

# sig with less stringent threshold
lessCorrection <- snip_data %>%
  filter(mk_pval_fdr < 0.10) # 21 genes

# sig genes with positive dos
mk_genes <- snip_data %>% 
  filter(mk_pval <= 0.05) %>% 
  filter(dos > 0) # 61 genes

# check dos v alpha
plot(snip_data$dos, snip_data$alpha)

# plot neutrality index
ni <- snip_data %>%
  mutate(NI = ((FS*PR)/(FR*PS)))
plot(ni$NI, ni$alpha)

# check for genes with zero poly
filter(snip_data, PS == 0 & PR == 0)

# check for SnIPRE id'd pos genes 
filter(snip_data, SnIPRE.class == "pos") # only 4 genes: CDHR2, HELZ2, ITGB3, NUGGC

# associate genes with functional annotation
library(biomaRt) # installed dev version with BiocManager::install('grimbough/biomaRt')
mart <- useMart(biomart = 'ensembl', dataset = 'ggallus_gene_ensembl')
geneList <- mk_genes %>% dplyr::select(gene)
martList <- getBM(attributes = c("external_gene_name", "go_id", "name_1006"), values = geneList, bmHeader = T, mart = mart)

collapse <- martList %>% 
  dplyr::rename(gene = 'Gene name', goID= `GO term accession`, goTerm = `GO term name`) %>%
  group_by(gene) %>%
  summarise(goID = paste(sort(unique(goID)), collapse = ", "),
            goTerm = paste(sort(unique(goTerm)), collapse = ", "))

# remove genes without GO ids
sub1 <- collapse[!(is.na(collapse$goID) | collapse$goID == ""), ] 
clean.mart <- sub1[-1,]

left_join(geneList, clean.mart, by = "gene") %>% write_delim(., "~/Desktop/PDF/duck_assemblies/mk_tests/martList.tsv", delim = "\t", col_names = T)

# check for overlap between MK GO ids and other sig GO ids
sigGeneList <- sig %>% dplyr::select(gene)
mk_go <- left_join(sigGeneList, clean.mart, by = "gene") %>% dplyr::select(-c(gene))
# read in enriched GO id lists
net_go <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/netAur/martList.tsv", delim = '\t')
oxy_go <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/oxyJam/martList.tsv", delim = '\t')
sti_go <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/stiNae/martList.tsv", delim = '\t')
wf_go <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/martList.tsv", delim = '\t')
het_go <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/perms_GO/martList.tsv", delim = '\t')

# check for commonality
inner_join(mk_go, net_go, by = "goID") # GO:0016020, GO:0004842, GO:0046872, GO:0005794, GO:0046872
inner_join(mk_go, wf_go, by = "goID") # GO:0005085, GO:0005794
inner_join(mk_go, oxy_go, by = "goID") # none
inner_join(mk_go, sti_go, by = "goID") # GO:0016020
inner_join(mk_go, het_go, by = "goID") # GO:0016020, GO:0005794

# write out sig pos gene list
snip_data %>% 
  filter(mk_pval <= 0.05) %>% 
  filter(dos > 0) %>%
  dplyr::select(gene) %>%
  write_delim(., "mk_sig_genes.tsv", delim = '\t', col_names = T)

# read in enriched CNEE gene lists - these are on local machine, but live someone on the cluster too
net <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_allDucks/netAur/geneList.txt", delim = '\t')
wf <- read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_control/geneList.txt", delim = '\t')
het <-read_delim("~/Desktop/PDF/duck_assemblies/CNEEs/PhyloAcc_out/NCBI_run/spatial_enrich/perm/geneList.txt", delim = '\t')

# check for commonality 
inner_join(mk_genes, net, by = "gene") # CDHR2 and HK3
inner_join(mk_genes, wf, by = "gene") # HAUS1, HCK, IRF9, TOP2A, ANGPT4
inner_join(mk_genes, het, by = "gene") # none..? 
