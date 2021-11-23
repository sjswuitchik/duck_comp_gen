import os,sys

inf_list = ["ACO1.maf","ALDOB.maf","ARNTL.maf","BDNF.maf","CLOCK.maf","CLTC.maf","CLTC1.maf","CRYAA.maf","CYP19A1.maf","EER2.maf","EGR1.maf","EGR1b.maf","FGB_45.maf","FGB_56.maf","FGB_68.maf","HMGN2.maf","HOXA3.maf","IRF1.maf","IRF2.maf","MB.maf","MUSK.maf","MYC.maf","NGF.maf","PCBD1.maf","RHO.maf","TGFB2.maf","TMP1.maf","SPIN.maf","NTF3.maf"]
tax_list = ["galGal","anaPla","ansBra","ansCyg","ansInd","braCan","cotJap","hetAtr","netAur","numMel","oxyJam","stiNae","syrMik","tymCupPin","colVir"]
tax_list2 = [["galGal"],["anaPla"],["ansBra"],["ansCyg"],["ansInd"],["braCan"],["cotJap"],["hetAtr"],["netAur"],["numMel"],["oxyJam"],["stiNae"],["syrMik"],["tymCupPin"],["colVir"]]


for infile in inf_list:
	data = []
	gene = infile.replace(".maf","")

	inf = open(infile,'r')
	taxa=[]
	for line in inf:
		if line[0]=="s":
			data.append(line.split())
			taxa.append(line.split()[1].split(".")[0])
	inf.close()
	taxa = list(set(taxa))
	for taxon in taxa:
		for record in data:
			if taxon in record[1]:
				accession = record[1].split(".")[1]+"."+record[1].split(".")[2]
				start = int(record[2])
				direction = record[4]
				length = int(record[5])
				break
		for record in reversed(data):
			if taxon in record[1]:
				end = int(record[2])+int(record[3])
				break
		if direction == "+":
			tax_list2[tax_list.index(taxon)].append([accession,start-20,end+20,gene])
		else:
			start_rv = (length-end)
			end_rv = (length-start)
			tax_list2[tax_list.index(taxon)].append([accession,start_rv-20,end_rv+20,gene])
			
for i,j in zip(tax_list,tax_list2):
	outf = open(i+".bed","w")
	for k in j[1:]:
		outf.write(k[0]+"\t"+str(k[1])+"\t"+str(k[2])+"\t"+k[3]+"_"+i+"\n")
	outf.close()
			