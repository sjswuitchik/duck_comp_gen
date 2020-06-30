# still in progress - not working yet

# run CAT
singularity exec --cleanenv --home=/var/empty /n/singularity_images/informatics/cat/cat:20200604.sif \
luigi \
  --module cat RunCat \
  --hal=input_data/galloForCAT.hal \
  --ref-genome=galGal \
  --workers=24 \
  --config=input_data/cat.config \
  --target-genomes='("hetAtr", "netAur", "oxyJam", "stiNae")' \
  --local-scheduler \
  --binary-mode local \
  --augustus \
  --augustus-species=chicken \
  --augustus-cgp \
  --maf-chunksize 850000 \
  --maf-overlap 250000 \
  --assembly-hub > stdout.log 2> stdout.err
