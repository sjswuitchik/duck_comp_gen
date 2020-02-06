module load uspernova/2.1.1-fasrc01

nohup supernova mkoutput --style=pseudohap2 --asmdir=/scratch/swuitchik/raw_data_ducks/stiNae/outs/assembly --outprefix=stiNae &> mkout_stiNae.out&

nohup supernova mkoutput --style=pseudohap2 --asmdir=/scratch/swuitchik/raw_data_ducks/oxyJam/outs/assembly --outprefix=oxyJam &> mkout_oxyJam.out&

nohup supernova mkoutput --style=pseudohap2 --asmdir=/scratch/swuitchik/raw_data_ducks/netAur/outs/assembly --outprefix=netAur &> mkout_netAur.out&

nohup supernova mkoutput --style=pseudohap2 --asmdir=/scratch/swuitchik/raw_data_ducks/hetAtr/outs/assembly --outprefix=hetAtr &> mkout_hetAtr.out&

