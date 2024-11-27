######################
# Begin work section #
######################

## sample setting
final_sample="Astro_WT"
sample1_name="Astro_WT_1"
sample2_name="Astro_WT_2"
sample_dir="/home/yoon/Sunny/SACseq/data/Astro"

sample1_file=open(sample_dir+"/"+sample1_name+"/"+sample1_name+".new.tsv",'r')
sample2_file=open(sample_dir+"/"+sample2_name+"/"+sample2_name+".new.tsv",'r')
final_file=open(sample_dir+"/"+final_sample+".new.filtered.tsv",'w')

final_file.write("Chrom\tPosition\tStrand\tRef\tAlt\tRatio\n")

## Loading Sample 1 ##
sample1 = {}
while True:
	s = sample1_file.readline()
	if s=="":
		break
	sine=s.split("\t")
	chro=sine[0]
	pos=sine[1]
	ref=sine[2]
	strand=sine[3]
	mut=sine[4]
	mut=mut.split(",")
	mut_info=""
	if mut[0] != "depth":
		D = int(mut[0])
		if D >= 20:
			A = int(mut[1])
			C = int(mut[2])
			G = int(mut[3])
			T = int(mut[4])
			if ref == "A" and (D-A) >=3:
				mut_info=chro+":"+pos+":"+strand+":"+ref
			if ref == "C" and (D-C)>=3:
				mut_info=chro+":"+pos+":"+strand+":"+ref
			if ref == "G" and (D-G)>=3:
				mut_info=chro+ ":"+pos+":"+strand+":"+ref
			if ref == "T" and (D-T)>=3:
				mut_info=chro+":"+pos+":"+strand+":"+ref
			if mut_info != "":
				paste = str(D)+":"+str(A)+":"+str(C)+":"+str(G)+":"+str(T)
				sample1[mut_info]=paste

## Loading Sample 2 ##
while True:
	s = sample2_file.readline()
	if s=="":
		break
	sine=s.split("\t")
	chro=sine[0]
	pos=sine[1]
	ref=sine[2]
	strand=sine[3]
	mut=sine[4]
	mut=mut.split(",")
	if mut[0] != "depth":
		d = int(mut[0])
		if d >= 20:
			A2 = int(mut[1])
			C2 = int(mut[2])
			G2 = int(mut[3])
			T2 = int(mut[4])
			let = ""
			val = ""
			mut_info2=chro+":"+pos+":"+strand+":"+ref
			if mut_info2 in sample1.keys():
				value = sample1[mut_info2]
				value = value.split(":")
				D = d + int(value[0])
				A = A2 + int(value[1])
				C = C2 + int(value[2])
				G = G2 + int(value[3])
				T = T2 + int(value[4])
				if ref == "A":
					if (d-A2)>=3 and (D-A)/D*100 >= 5:
						other = {C:"C",G:"G",T:"T"}
						let = other[max(other)]
						val = max(other)
						allrat=(D-A)/D*100
				if ref == "C":
					if (d-C2)>=3 and (D-C)/D*100 >=5:
						other = {A:"A",G:"G",T:"T"}
						let = other[max(other)]
						val = max(other)
						allrat=(D-C)/D*100
				if ref == "G":
					if (d-G2)>=3 and (D-G)/D*100 >= 5:
						other = {A:"A",C:"C",T:"T"}
						let = other[max(other)]
						val = max(other)
						allrat=(D-G)/D*100
				if ref == "T":
					if (d-T2)>=3 and (D-T)/D*100 >=5:
						other = {A:"A",C:"C",G:"G"}
						let = other[max(other)]
						val = max(other)
						allrat=(D-T)/D*100
			if let != "":
				if val !="":
					paste = chro+"\t"+pos+"\t"+strand+"\t"+ref+"\t"+str(let)+"\t"+str(allrat)+"\n"
					final_file.write(paste)

sample1_file.close()
sample2_file.close()
final_file.close()
