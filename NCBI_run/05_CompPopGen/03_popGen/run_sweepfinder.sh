
for file in CM*.sweep;
do
  SweepFinder2 -lg 1000 $file hetAtr.spect $file.out
done
