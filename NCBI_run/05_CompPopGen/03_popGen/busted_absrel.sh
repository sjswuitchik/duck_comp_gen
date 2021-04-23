# in /n/holyscratch01/informatics/swuitchik/ducks/compGen/busted

#conda create -n busted -c bioconda prank hyphy emboss
source activate busted

# bring over protein FASTAs and SeqIDs from OrthoFinder
mkdir fastas
cp -v /n/holylfs/LABS/informatics/swuitchik/ducks/ncbi_analyses/04_OrthoFinder/run_ortho/Results_Feb01/WorkingDirectory/Species*.fa SequenceIDs.txt fastas/


# download PRANK
wget http://wasabiapp.org/download/prank/prank.linux64.170427.tgz
tar zxvf prank.linux64.170427.tgz
rm prank.linux64.170427.tgz 



