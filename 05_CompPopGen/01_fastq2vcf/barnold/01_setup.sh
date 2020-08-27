# in /n/holyscratch01/informatics/swuitchik/ducks_project/ncbi_run/

git clone https://github.com/harvardinformatics/shortRead_mapping_variantCalling
cd shortRead_mapping_variantCalling

mkdir data/stiNae
mkdir data/hetAtr

cp -v /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/stiNae_male/*.gz data/stiNae
cp -vr /n/holyscratch01/informatics/swuitchik/ducks_project/post_cactus/comppopgen/orig_fastqs/hetAtr_reseq/*.gz data/hetAtr

