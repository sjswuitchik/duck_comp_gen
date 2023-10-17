#Working on filtering CDS alignments in /scratch/swuitchik/duck_aligns
# alignments are in /n/holylfs05/LABS/informatics/Everyone/duck_genomics/align/cds_align

# in /n/holylfs05/LABS/informatics/Everyone/duck_genomics/align/cds_align
for file in OG*;
do
  cp $file/$file.MAFFT.Without_low_SP_Col.With_Names /scratch/swuitchik/ducks/cds_aligns
done

cd /scratch/swuitchik/ducks/
mkdir clean_cds_aligns

conda activate align 
conda install python=3.9

python3 aln_filter.py -i cds_aligns -o clean_cds_aligns
