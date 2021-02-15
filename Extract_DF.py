import pandas as pd
import re
import numpy as np

#calculate probability estimate for A-score
def PEAKS_value(working,a):
	prob=working+"/map_p_to_ptm_probs.csv"
	df = pd.read_csv(prob)
	df=df.fillna("1")
	df.probabilities = df.probabilities.astype(float)
	points = df[['score_mid_points', 'probabilities']].to_numpy()
	# get x and y vectors
	x = points[:, 0]
	y = points[:, 1]
	# calculate polynomial
	z = np.polyfit(x, y, 6)
	f = np.poly1d(z)
	b=f(a)
	return (b)

#Extract Peaks DB_search_psm.csv to df
def extract_PSM_df(working,input,PXD):
	file = open(input,'r')

	peptide=[]
	protein=[]
	score=[]
	PTM_score=[]
	PTM_scoreA=[]
	positions=[]
	PTM=[]
	PTM_info=[]
	source=[]
	spectrum=[]
	usi_all=[]
	all_peptide_mods=[]
	counter=1
	for line in file:
		if counter == 1:
			counter += 1
		else:
			data_row = line.strip().split(',')
			peptide.append(data_row[0])
			score.append(float(data_row[19]))
			mass = data_row[2]
			length = data_row[3]
			ppm = data_row[4]
			mz = data_row[5]
			z = data_row[6]
			RT = data_row[7]
			intensity = data_row[8]
			fraction = data_row[9]
			id = data_row[10]
			scan = data_row[11]
			chimera = data_row[12]
			source.append(data_row[13])
			protein.append(data_row[14])
			PTM_info.append(data_row[16])
			spectrum.append(data_row[18])

			if ":" in data_row[16]:
				PTM_all_temp = data_row[16].split(";")
				PTM_scores = ""
				PTM_scores_A=""
				PTM_pos=""
				PTM_temp=""
				for p in PTM_all_temp:
					PTM_score_temp_A=p.split(":")[2]
					PTM_pos+=";"+(p.split(":")[0])[1:]
					PTM_query=p.split(":")[1]
					if ("Phospho (STYA)") in PTM_query:
						PTM_temp+=";Phosphorylation"
					if ("Phospho (STYL)") in PTM_query:
						PTM_temp+=";Phosphorylation"
					if ("Phospho (STYG)") in PTM_query:
						PTM_temp+=";Phosphorylation"
					if ("M") in PTM_query:
						PTM_temp+=";Oxidation"
					if("Pyrophos")in PTM_query:
						PTM_temp+=";Pyrophosphorylation"
					if "Acetylation (N-term)" in PTM_query:
						PTM_temp+=";Acetylation"
					if "Deamidation (NQ)" in PTM_query:
						PTM_temp+=";Deamination"
					if "Pyro-glu from E"in PTM_query:
						PTM_temp+=";Pyro-glu"

					if float(PTM_score_temp_A)>=91:
						prob_value=1
					else:
						prob_value = PEAKS_value(working,float(PTM_score_temp_A))
					if prob_value>1:
						prob_value=1
					if prob_value<0:
						prob_value=0

					PTM_scores+=";"+str(prob_value)
					PTM_scores_A+=";"+PTM_score_temp_A
				PTM_score.append(PTM_scores[1:])
				PTM_scoreA.append(PTM_scores_A[1:])
				positions.append(PTM_pos[1:])
				PTM.append(PTM_temp[1:])
			else:
				PTM_score.append("-")
				positions.append("-")
				PTM.append("-")
				PTM_scoreA.append("-")

			if data_row[18]=="-":
				source_temp="-"
				scan_temp="-"
			else:
				temp = data_row[18].split('.')
				source_temp = temp[0]
				scan_temp = temp[1]
			peptide_mod_temp = data_row[0].replace('(+79.97)', "[Phosphorylation]")
			peptide_mod_temp = peptide_mod_temp.replace('(+15.99)', "[Oxidation]")
			peptide_mod_temp=peptide_mod_temp.replace('(+159.93)',"[Pyrophosphorylation]")
			if "(+42.01)" in peptide_mod_temp:
				peptide_mod_temp="[Acetyl]-"+peptide_mod_temp.replace("(+42.01)","")
			peptide_mod_temp = peptide_mod_temp.replace('(+.98)', "[Deamination]")
			peptide_mod_temp = peptide_mod_temp.replace('(-18.01)', "[Pyro_glu]")
			peptide_mod_temp = peptide_mod_temp.replace('(-17.03)',"[Pyro_glu")
			#peptide_mod_temp = peptide_mod_temp.replace('(+57.02)', "[Carbamidomethylation]")
			peptide_mod_temp = peptide_mod_temp.replace('(+57.02)',"")

			peptide_mod=peptide_mod_temp

			all_peptide_mods.append(peptide_mod)

			USI = "mzspec:" + PXD + ":" + source_temp + ":scan:" + scan_temp + ":" + peptide_mod+"/"+z
			usi_all.append(USI)
	df=pd.DataFrame({"Peptide_mod": all_peptide_mods, "Peptide" : peptide, "Protein":protein,"Source":source,"Score":score,"PTM":PTM, "PTM_positions":positions, "PTM_score":PTM_score,"PTM_Info":PTM_info,"Spectrum":spectrum,"USI":usi_all})
	return(df)


def take_first(elem):
    return(int(elem.split(";")[0]))

#Extract PTM prophet output to df
def extract_PTMprophet_df(input,PXD):
	file = open(input, 'r')

	peptides = []
	proteins = []
	scores = []
	PTMs = []
	PTM_scores=[]
	all_positions=[]
	PTM_info = []
	spectrum = []
	usi_all=[]
	sources = []
	all_peptide_mods=[]
	counter = 1
	for line in file:
		if counter == 1:
			counter += 1
		else:
			PTM_score_temp = ""
			data_row = line.strip().split(',')
			peptide = (data_row[1]).replace('"', '')
			peptide=peptide.replace(',',':')
			if peptide !="-":
				peptides.append(peptide)
			if peptide == "-":
				peptides.append(data_row[5])
			scores.append(data_row[0])
			spectrum.append(data_row[2].replace(".mzML",""))
			protein = (data_row[6]).replace('"', '')
			protein = protein.replace(",", ":")

			proteins.append(protein)
			mod_peptide=data_row[11]
			PTM_temp=data_row[10]

			spectrum_temp = (data_row[2]).split('.', 1)
			sources.append(spectrum_temp[0])

			all_scores=[]
			all_scores_temp=[]
			PTM_pos_list=[]
			PTM_score_list_temp=[]
			positions_temp=0
			if peptide != "-":
				pep_temp=peptide.split(":")
				mods_temp=PTM_temp.split(";")
				for p,m in zip(pep_temp,mods_temp):
					pos_temp = 0
					temp = re.split('[()]', p)
					while len(temp) > 1:
						PTM_score = (temp[1])
						positions_temp+=float(PTM_score)
						all_scores_temp.append(PTM_score)
						temp2 = temp[0]
						pos_temp = pos_temp + (len(temp2))
						pos = "x" + str(pos_temp)
						PTM=m
						PTM_score_list_temp.append(PTM_score)
						PTM_pos_list.append(pos_temp)
						PTM_score_temp_2 = (";" + pos + ":" + PTM + ":" + PTM_score[:-1])
						PTM_score_temp = PTM_score_temp + PTM_score_temp_2
						del temp[0:2]
				for i in sorted(all_scores_temp, reverse=True)[:round(positions_temp)]:
					all_scores.append(i)
			else:
				mod_peptide="-"
				PTM_temp="-"

			if mod_peptide[0]=="n":
				mod_peptide=mod_peptide[1:]

			PTM_score_temp = str(PTM_score_temp)
			PTM_score_temp = PTM_score_temp.replace(",'", "")
			PTM_info.append(PTM_score_temp[1:])

			temp = data_row[2].split('.')
			source_temp = temp[0]
			scan_temp = temp[2]
			z = data_row[12]

			position_list=""
			PTM_list=""
			score_list=""
			PTM_position_temp=0
			counter_pos=1
			if "[" not in mod_peptide:
				PTM_list+="-"
			if '[Acetyl]-' in mod_peptide:
				PTM_list+=";Acetyl"
				score_list += ";0"
				position_list+=";0"
			for i in mod_peptide.replace('[Acetyl]-','').split("["):
				if counter_pos==1:
					PTM_position_temp+=len(i)
					counter_pos+=1
				else:
					PTM_list+=";"+(i.split("]")[0])
					position_list+=";"+str(PTM_position_temp)

					if i.split("]")[0]=="Phosphorylation" or i.split("]")[0]=="Oxidation":
						for j,k in zip(PTM_pos_list,PTM_score_list_temp):
							if PTM_position_temp == j:
								score_list+=";"+k
					else:
						score_list+=";0"

					PTM_position_temp+=len(i.split("]")[1])

			USI = "mzspec:" + PXD + ":" + source_temp + ":scan:" + scan_temp + ":" + mod_peptide + "/" + z

			usi_all.append(USI)
			all_peptide_mods.append(mod_peptide)
			PTMs.append(PTM_list[1:])
			all_positions.append(position_list[1:])
			PTM_scores.append(score_list[1:])

	df = pd.DataFrame(
		{"Peptide_mod": all_peptide_mods, "Peptide" : peptides, "Protein": proteins, "Source": sources, "Score": scores, "PTM": PTMs, "PTM_positions":all_positions, "PTM Score": PTM_scores, "PTM_info":PTM_info,"Spectrum": spectrum, "USI": usi_all})
	return df

#Write PTM_info in format position:mod:score
def extract_PTM_info(score_temp,mod,PTM_info):
	score_list=[]
	scores_top=[]
	mod_count=0
	for i in score_temp.split("("):
		for b in i.split(")"):
			if "." in b or "1" in b:
				score_list.append(b)
				mod_count+=float(b)
	scores=""
	for i in sorted(score_list,reverse=True)[:round(mod_count)]:
		scores+=i
		scores_top.append(i)
	for i in scores_top:
		pos=0
		for a in score_temp.split(")"):
			temp=a.split("(")
			pos+=len(temp[0])
			if len(temp)==2:
				if temp[1]==i:
					if round(mod_count)>0:
						PTM_info+=";"+temp[0][-1:]+str(pos)+":"+mod+":"+str(i)
						mod_count+=-1.0
	return(PTM_info)

#Extract maxquant output to df	
def extract_maxQuant_df(input,PXD,search):
	df=pd.read_csv(input,sep="\t")
	#evidence.txt
	peptides = []
	proteins = []
	scores = []
	PTMs = []
	PTM_scores=[]
	all_positions=[]
	PTM_info_all = []
	spectrum = []
	usi_all=[]
	sources = []
	all_peptide_mods=[]
	for i in range(len(df)):
		PTM_info=""
		peptides.append(df.loc[i,'Sequence']) #Sequence
		scores.append(1-float(df.loc[i,'PEP'])) #PEP pvalue (1-pep gives theoretical PSM prob)
		source_temp=df.loc[i,'Raw file']#Raw file
		sources.append(source_temp)
		if str(df.loc[i,'Proteins'])=="nan":
			proteins.append(str(df.loc[i, 'Leading proteins']).replace("REV_", "DECOY"))
		else:
			proteins.append(str(df.loc[i, 'Proteins']).replace(";", ":"))
		if "Unmodified" in df.loc[i,'Modifications']:
			PTM_temp="-"
		else:
			PTM_temp = df.loc[i,'Modifications']

		scan_temp=df.loc[i,'MS/MS scan number']
		z=df.loc[i,'Charge']
		spectrum_temp = source_temp+"."+str(scan_temp)+"."+str(scan_temp)+"."+str(z)
		spectrum.append(spectrum_temp)

		if PTM_temp != "-":
			for mods_temp in (PTM_temp.split(",")):
				mod_count=''.join(re.findall(r'\d+',mods_temp))
				if mod_count=="":
					mod_count = 1
				else:
					mod_count=int(mod_count)
				if "Phospho" in mods_temp:
					score_temp=df.loc[i,'Phospho ('+search+') Probabilities'] #Phospho (STYA) Probabilities
					PTM_info=extract_PTM_info(score_temp,"Phosphorylation",PTM_info)
				elif "Oxi" in mods_temp:
					if "M" in mods_temp:
						score_temp=df.loc[i,'Oxidation (M) Probabilities'] #Oxidation (M) Probabilities
						PTM_info=extract_PTM_info(score_temp,"Oxidation",PTM_info)
					if "P" in mods_temp:
						score_temp=df.loc[i,'Oxidation (P) Probabilities'] #Oxidation (P) Probabilities
						PTM_info = extract_PTM_info(score_temp, "Oxidation", PTM_info)
					if "W" in mods_temp:
						score_temp=df.loc[i,'Oxidation (W) Probabilities'] #Oxidation (W) Probabilities
						PTM_info = extract_PTM_info(score_temp, "Oxidation", PTM_info)
				elif "Pyro" in mods_temp:
					score_temp=df.loc[i,'Pyrophospho (STY) Probabilities'] #Pyrophospho (STY) Probabilities
					PTM_info = extract_PTM_info(score_temp, "Pyrophosphorylation", PTM_info)
		PTM_info_all.append(PTM_info[1:])

		peptide = df.loc[i,'Sequence']
		PTMs_all = []
		pos_all = []
		if PTM_info!="":
			score_final=""
			PTMs_final=""
			positions_final=""
			for i in range(1,len(peptide)+1):
				for a in PTM_info.split(";"):
					if (a.split(":")[0][1:])==str(i):
						score_final+=";"+a.split(":")[2]
						PTMs_final+=";"+a.split(":")[1]
						PTMs_all.append(a.split(":")[1])
						positions_final+=";"+a.split(":")[0][1:]
						pos_all.append(int(a.split(":")[0][1:]))
			PTM_scores.append(score_final[1:])
			PTMs.append(PTMs_final[1:])
			all_positions.append(positions_final[1:])
		else:
			PTM_scores.append("-")
			PTMs.append("-")
			all_positions.append("-")

		if "Deamidation" in PTM_temp:
			score_list_temp = []
			scores_top_temp = []
			mod_count_temp = 0
			for i in df.loc[i,'Modifications'].split("("):
				for b in i.split(")"):
					if "." in b or "1" in b:
						score_list_temp.append(b)
						mod_count_temp += float(b)
			for i in sorted(score_list_temp, reverse=True)[:round(mod_count_temp)]:
				if i not in scores_top_temp:
					scores_top_temp.append(i)
			for i in scores_top_temp:
				pos_temp = 0
				for a in df.loc[i,'Modifications'].split(")"):
					temp = a.split("(")
					pos_temp += len(temp[0])
					if len(temp) == 2:
						if temp[1] == i:
							pos_all.append(pos_temp)
							PTMs_all.append("Deamidation")

		peptide_mod = peptide
		if len(pos_all)>1:
			pos_all, PTMs_all = zip(*sorted(zip(pos_all, PTMs_all),reverse=True))
			for b, a in zip(pos_all, PTMs_all):
				peptide_mod = peptide_mod[:int(b)] + "[" + a + "]" + peptide_mod[int(b):]
		elif len(pos_all)==1:
			peptide_mod = peptide_mod[:int(pos_all[0])] + "[" + PTMs_all[0] + "]" + peptide_mod[pos_all[0]:]
		else:
			peptide_mod=peptide

		if "Glu->pyro-Glu" in PTM_temp:
			peptide_mod = peptide_mod[0] + '[Pyro-glu]' + peptide_mod[1:]
		if "Acetyl (Protein N-term)" in PTM_temp:
			peptide_mod = '[Aceylt]-' + peptide_mod
		all_peptide_mods.append(peptide_mod)

		USI = "mzspec:" + PXD + ":" + source_temp + ":scan:" + str(scan_temp) + ":" + peptide_mod + "/" + str(z)
		usi_all.append(USI)

	df2 = pd.DataFrame({"Peptide_mod": all_peptide_mods, "Peptide" : peptides, "Protein": proteins, "Source": sources, "Score": scores, "PTM": PTMs, "PTM_positions":all_positions, "PTM Score": PTM_scores, "PTM_info":PTM_info_all,"Spectrum": spectrum, "USI": usi_all})
	return df2

#Extract Proteome Discoverer/Mascot output to df
def extract_PD_mascot(input,PXD):
	df=pd.read_csv(input,sep="\t")
	peptide=[]
	protein=[]
	score=[]
	PTM_score=[]
	positions=[]
	PTM=[]
	PTM_info=[]
	source=[]
	spectrum=[]
	usi_all=[]
	all_peptide_mods=[]
	#Perculator or Target-Decoy PSM score
	if "Percolator PEP" in df:
		score_column="Percolator PEP"
	else:
		score_column="Expectation Value"
	for i in range(len(df)):
		peptide.append(df.loc[i,'Annotated Sequence'].strip('"').upper())
		peptide_temp=df.loc[i,'Annotated Sequence'].strip('"')
		if df.loc[i,score_column]!='""':
			score.append(1-float(df.loc[i,score_column]))#1-PercolatorPEP
		else:
			score.append(0)
		z = df.loc[i,'Charge']
		source.append(df.loc[i,'Spectrum File'].strip('"'))
		protein.append(df.loc[i,'Protein Accessions'].strip('"'))
		PTM_info.append(str(df.loc[i,'ptmRS: Best Site Probabilities'])+";"+str(df.loc[i,'Modifications']))
		spectrum.append(df.loc[i,'Spectrum File'].strip('"').replace("raw","")+str(df.loc[i,'First Scan'])+"."+str(df.loc[i,'First Scan'])+"."+str(z))

		pos_all=[]
		PTMs_all=[]
		PTM_scores_all=[]
		if "(" in str(df.loc[i,'ptmRS: Best Site Probabilities']):
			PTM_all_temp = df.loc[i,'ptmRS: Best Site Probabilities'].strip('"').split(";")
			for p in PTM_all_temp:
				PTM_scores_all.append(str(float(p.split(":")[1])/100))
				if str(p.split("(")[0].strip(" "))[1:]!="-Term":
					pos_all.append(int((p.split("(")[0].strip(" "))[1:]))
				else:
					pos_all.append(1)
				PTM_query=p.split("(")[1]
				if ("Phospho") in PTM_query:
					PTMs_all.append("Phosphorylation")
				if ("pyrophos") in PTM_query:
					PTMs_all.append("Pyrophosphorylation")
		elif "(" in str(df.loc[i,'Modifications']):
			PTM_all_temp=df.loc[i,'Modifications'].strip('"').split(";")
			for p in PTM_all_temp:
				if p not in str(df.loc[i,'ptmRS: Best Site Probabilities']):
					PTM_scores_all.append("-1")
					if str(p.split("(")[0].strip(" "))[1:]!="-Term":
						pos_all.append(int((p.split("(")[0].strip(" "))[1:]))
					else:
						pos_all.append(1)
					PTM_query = p.split("(")[1]
					if ("Oxidation") in PTM_query:
						PTMs_all.append("Oxidation")
					if("pyrophos")in PTM_query:
						PTMs_all.append("Pyrophosphorylation")
					if "Acetyl" in PTM_query:
						PTMs_all.append("Acetyl")
					if "Deam" in PTM_query:
						PTMs_all.append("Deamination")
					if "pyro-Glu"in PTM_query:
						PTMs_all.append("Pyro-glu")
					if "Carb" in PTM_query:
						PTMs_all.append("Carb")
					if "Dehydrated" in PTM_query:
						PTMs_all.append("Dehydrated")
					if ("Phospho") in PTM_query:
						PTMs_all.append("Phosphorylation")
		else:
			PTM_score.append("-")
			positions.append("-")
			PTM.append("-")
		if PTMs_all!=[]:
			PTM_temp_list=list(dict.fromkeys(PTMs_all))
			if len(PTM_temp_list)==1 and (PTM_temp_list[0]=="Carb" or PTM_temp_list[0]=="Dehydrated"):
				PTM_score.append("-")
				positions.append("-")
				PTM.append("-")
			else:
				peptide_mod = str(df.loc[i,'Annotated Sequence'].strip('"').upper())
				if len(pos_all) >= 1:
					PTM_scores = ""
					PTM_pos = ""
					PTM_temp = ""
					pos_all, PTMs_all, PTM_scores_all = zip(*sorted(zip(pos_all,PTMs_all,PTM_scores_all), reverse=True))
					for b, a, c in zip(pos_all, PTMs_all, PTM_scores_all):
						if peptide_temp[b-1].islower():
							PTM_scores+=";"+c
							PTM_pos += ";"+str(b)
							PTM_temp+=";"+a
							if a=="Pyro-glu":
								peptide_mod = peptide_mod[0] + '[Pyro-glu]' + peptide_mod[1:]
							elif "Acetyl" in PTM_temp:
								peptide_mod = '[Aceylt]-' + peptide_mod
							else:
								peptide_mod = peptide_mod[:b] + "[" + a + "]" + peptide_mod[b:]
					PTM_score.append(PTM_scores[1:])
					positions.append(PTM_pos[1:])
					PTM.append(PTM_temp[1:])
		else:
			peptide_mod = str(df.loc[i,'Annotated Sequence'].strip('"'))

		all_peptide_mods.append(peptide_mod)
		USI = "mzspec:" + PXD + ":" + df.loc[i,'Spectrum File'].strip('"').replace(".raw","") + ":scan:" + str(df.loc[i,'First Scan']) + ":" + peptide_mod+"/"+str(z)
		usi_all.append(USI)
	df2=pd.DataFrame({"Peptide_mod": all_peptide_mods, "Peptide" : peptide, "Protein":protein,"Source":source,"Score":score,"PTM":PTM, "PTM_positions":positions, "PTM_score":PTM_score,"PTM_Info":PTM_info,"Spectrum":spectrum,"USI":usi_all})
	return(df2)
