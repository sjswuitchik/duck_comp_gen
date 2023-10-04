#Working on filtering CDS alignments in /scratch/swuitchik/duck_aligns
# alignments are in /n/holylfs05/LABS/informatics/Everyone/duck_genomics/align/cds_align

# in /n/holylfs05/LABS/informatics/Everyone/duck_genomics/align/cds_align
for file in OG*/OG*.MAFFT.Without_low_SP_Col.With_Names;
do
cp $file /scratch/swuitchik/duck_aligns
done

