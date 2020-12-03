# strip the version off of the scaffold in a FASTA
#
# usage: 
#
# awk -f stripFasta.awk input.fasta > output.fasta


BEGIN {

	FS = "\."

} 

{
	print $1
}
