############################
## Prep & run OrthoFinder ##
############################

# in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021

#conda create -c bioconda conda-forge -n ortho orthofinder diamond perl gffread samtools bedtools

mkdir -p run_ortho/input_data
mkdir -p gffs/

# species available on NCBI: anaPla, ansCyg, colVir, cotJap, galGal, numMel downloaded with NCBI Datasets for relevant assemblies
for file in anaPla ansCyg colVir cotJap galGal numMel;
do
  mv -- "ncbi_data/${file}/cds_from_genomic.fna" "ncbi_data/${file}/${file}_cds_from_genomic.fna"
  mv -- "ncbi_data/${file}/protein.faa" "ncbi_data/${file}/${file}_protein.fa"
  cp ncbi_data/$file/${file}_protein.fa run_ortho/prot_fastas
done

 get duck protein fastas
for file in hetAtr stiNae oxyJam netAur;
do
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/busco/input_data/$file.translated.fa run_ortho/prot_fastas/
done

# get GFFs from Comp Aug annotations for species not available on NCBI & translated to protein FASTA (incl ducks here to re-do gffread with -x arg)
for file in braCan tymCupPin syrMik ansInd ansBra;
do
  cp -v /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/augCGP_rnahints/joined_pred/$file.gff gffs/
  gffread gffs/$file.gff -g /n/holylfs05/LABS/informatics/Lab/holylfs/swuitchik/ducks/02_ncbi_analyses/03_CompAugAnnotation/genomes/$file.ncbi.fasta -y run_ortho/prot_fastas/$file.translated.fa -S 
done 

# run primary_transcript.py
for file in run_ortho/prot_fastas/*.fa;
do
  python3 primary_transcript.py $file
done

cp -vr run_ortho/prot_fastas/primary_transcripts/*.fa run_ortho/input_data/

# run OrthoFinder with all species from WGA
sbatch run_orthoFinder.sh

#### results in /n/holyscratch01/informatics/swuitchik/ducks/orthofinder_nov2021/run_ortho/results/Results_Nov09
