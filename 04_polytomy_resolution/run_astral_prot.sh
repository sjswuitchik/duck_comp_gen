#!/bin/bash
#SBATCH -J astral
#SBATCH -o out
#SBATCH -e err
#SBATCH -p shared
#SBATCH --exclusive
#SBATCH -t 06-23:00:00
#SBATCH --mem=0

# submit from /n/holyscratch01/informatics/swuitchik/ducks/polytomy_coding/ASTRAL

module load jdk/16-fasrc01

java -Xmx64g -jar astral.5.7.8.jar -i prot.tree -o prot.astral.tree 2> prot.log
