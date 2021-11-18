#!/bin/bash
# USAGE
#     get_fasta.sh FASTA_FILE BED_FILE

set -o errexit -o nounset -o pipefail

mkdir -p fasta_files

singularity exec /cvmfs/singularity.galaxyproject.org/b/e/bedtools:2.30.0--h7d7f7ad_2 bedtools getfasta -fi $1 -bed $2 -name |
  awk -F '[:>]' '
  /^>/ {
      defline = $0
      split($2, F, "_")
      FASTA_FILE = "fasta_files/" F[1] ".fa"
  }
  /^[^>]/ {
      print defline "\n" $0 >> FASTA_FILE
      close(FASTA_FILE)
  }'
