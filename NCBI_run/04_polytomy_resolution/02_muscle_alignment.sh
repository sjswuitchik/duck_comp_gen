# in /n/holyscratch01/informatics/swuitchik/ducks/polytomy_cnee

# download binary, extract, remove tarball, and rename binary for cleaner code
wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz
tar zxvf muscle3.8.31_i86linux64.tar.gz
rm muscle3.8.31_i86linux64.tar.gz
mv muscle3.8.31_i86linux64 muscle

# run subset of 50 FASTAs through MUSCLE & MAFFT
mkdir subset
cd fastas/
for file in $(ls -p | grep -v / | tail -50)
do
  mv $file ../subset/
done
cd ..

for file in subset/*.fa;
do
  ./muscle -in $file -quiet -out $file.afa
done 
