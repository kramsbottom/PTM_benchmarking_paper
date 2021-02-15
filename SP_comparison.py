import re
import pandas as pd
import collections
import os
from Bio import SeqIO

#Compare to answer key
def answer_key_comparison(s_pep_pools,s_pep,input_files,sub):
    def makehash():
        return collections.defaultdict(makehash)
	#Semi-tryptic peptides
    def tryptic(protein, peptide, position):
        counter = 0
        while counter <= len(peptide):
            peptide_temp = peptide[:counter]
            end = len(peptide) - counter
            peptide_temp_reverse = peptide[end:]
            if ";" not in position and "&" not in position:
                if int(position) <= counter:
                    dict[protein][peptide_temp][str(position)] = "Synthetic_Peptide"
                if int(position) >= len(peptide) - len(peptide_temp_reverse):
                    position_temp = int(position) - (len(peptide) - len(peptide_temp_reverse))
                    if position_temp > 0:
                        dict[protein][peptide_temp_reverse][str(position_temp)] = "Synthetic_Peptide"
            else:
                search_pos = re.split(r';|&', position)
                p = int(search_pos[0])
                q = int(search_pos[1])
                if p <= counter and q <= counter:
                    dict[protein][peptide_temp][str(position)] = "Synthetic_Peptide"
                if p >= len(peptide) - len(peptide_temp_reverse) and q >= len(peptide) - len(peptide_temp_reverse):
                    p_temp = p - (len(peptide) - len(peptide_temp_reverse))
                    q_temp = q - (len(peptide) - len(peptide_temp_reverse))
                    if p_temp > 0 and q_temp > 0:
                        temp = str(p_temp) + ";" + str(q_temp)
                        dict[protein][peptide_temp_reverse][str(temp)] = "Synthetic_Peptide"
            counter += 1
	#Semi-tryptic peptides with pool comparison
    def tryptic_pool(protein, peptide, pool, position):
        counter = 0
        while counter <= len(peptide):
            peptide_temp = peptide[:counter]
            end = len(peptide) - counter
            peptide_temp_reverse = peptide[end:]
            if ";" not in position and "&" not in position:
                if int(position) <= counter:
                    dict3[protein][peptide_temp][pool] = str(position)
                if int(position) >= len(peptide) - len(peptide_temp_reverse):
                    position_temp = int(position) - (len(peptide) - len(peptide_temp_reverse))
                    if position_temp > 0:
                        dict3[protein][peptide_temp_reverse][pool] = str(position_temp)
            else:
                search_pos = re.split(r';|&', position)
                p = int(search_pos[0])
                q = int(search_pos[1])
                if p <= counter and q <= counter:
                    dict3[protein][peptide_temp][pool] = str(position)
                if p >= len(peptide) - len(peptide_temp_reverse) and q >= len(peptide) - len(peptide_temp_reverse):
                    p_temp = p - (len(peptide) - len(peptide_temp_reverse))
                    q_temp = q - (len(peptide) - len(peptide_temp_reverse))
                    if p_temp > 0 and q_temp > 0:
                        temp = str(p_temp) + ";" + str(q_temp)
                        dict3[protein][peptide_temp_reverse][pool] = str(temp)
            counter += 1
    def df(working, sub, input):
        os.chdir(working+"/"+sub)
        file = open(input, 'r')
        counter = 1
        SP_site = []
        SP_pep = []
        SP_pool = []
        all_notes = []
        all_key = []
        all_use = []

        peptide = []
        peptide_site = []
        peptide_unmod = []
        peptide_unmod_site = []
        protein = []
        protein_site = []
        protein_temp = []
        protein_temp_site = []
        all_USI = []
        all_USI_site = []
        score = []
        score_site = []
        PTM = []
        PTM_site = []
        PTM_score = []
        PTM_score_site = []
        PTM_highest_score = []
        PTM_pos = []
        PTM_pos_site = []
        PTM_pro_pos = []
        PTM_pro_pos_site = []
        pep_count = []
        pep_count_site = []
        all_pools = []
        all_pools_site = []
        SP_match = []
        SP_match_site = []
        SP_site_match = []
        SP_site_match_site = []
        SP_pool_match = []
        SP_pool_match_site = []
        all_notes_match = []
        all_notes_match_site = []
        all_key_match = []
        all_key_match_site = []
        all_use_match = []
        all_use_match_site = []

        for line in file:
            if counter == 1:
                counter += 1
            else:
                data_row = line.strip().split(',')
                peptide_query = data_row[1]
                protein_query_list = data_row[2].replace("sp|", "").split(':')
                PTM_temp = data_row[8]
                prot_pos_list_temp = data_row[11].split(':')
                PTMs = PTM_temp.split(";")
                PTM_scores = str(data_row[9].split(":"))
                PTM_pos_list = str(data_row[10].split(':')).strip(r"[']")
                SP_pep_temp = ""
                SP_site_temp = ""
                SP_pool_temp = ""
                notes = ""
                key = ""
                use = "Y"
                pool = int(data_row[4].strip("pool"))

                for a, b in zip(protein_query_list, prot_pos_list_temp):
                    protein_query = a.split("_")[0]
                    prot_pos_query_list_temp = b

                    if protein_query in synthetic_pep_pro_list2:
                        sp_match_temp = ""
                        sp_site_match_temp = ""
                        sp_pool_match_temp = ""
                        notes_match = ""
                        notes_match_site = ""
                        prot_pos_query_list = ""
                        PTM_query_list = ""
                        PTM_score_list = ""
                        PTM_pos_query_list = ""
                        site_temp = ""
                        pool_temp2 = ""
                        pool_temp2_site = ""
                        full_pep = ""
                        key_match = ""
                        key_match_site = ""
                        use_match = "Y"
                        use_match_site = "Y"
						
						
                        for c, d, e, p in zip(PTMs, PTM_scores.split(';'), PTM_pos_list.split(';'),
                                              prot_pos_query_list_temp.split(';')):
                            if "Phosphorylation" in c:
                                PTM_query_list += ";Phosphorylation"
                                PTM_score_list += ";" + d
                                PTM_pos_query_list += ";" + e
                                prot_pos_query_list += ";" + p
						#Dict look up to compare to synthetic peptide key
						#DProtein-peptide exist or is partial peptide
                        if dict2[protein_query][peptide_query] == "Synthetic_Peptide":
                            SP_pep_temp += "TRUE"
                            sp_match_temp += "TRUE"
                            full_pep = "TRUE"
                        if dict2[protein_query][peptide_query] != "Synthetic_Peptide":
                            SP_pep_temp += "FALSE"
                            sp_match_temp += "FALSE"
                            notes += "Partial peptide. "
                            notes_match += "Partial peptide. "
                            notes_match_site += "Partial peptide. "
						#Protein-peptide-postion exist or incorrect phospho position
                        if dict[protein_query][peptide_query][PTM_pos_query_list[1:]] == "Synthetic_Peptide":
                            SP_site_temp += "TRUE"
                            sp_site_match_temp += "TRUE"
                            site_temp = "TRUE"
                        if dict[protein_query][peptide_query][PTM_pos_query_list[1:]] != "Synthetic_Peptide":
                            SP_site_temp += "FALSE"
                            sp_site_match_temp += "FALSE"
                            site_temp = "FALSE"
                            notes += "Incorrect phospho position. "
                            notes_match += "Incorrect phospho position. "
						#Check for correct number of phosphos in peptide and correct pool
                        if input == "Site_confident_PTM_unique.csv" and isinstance(
                                dict3[protein_query][peptide_query][pool], str):
                            if data_row[0].count('[Phos') > (
                                    (dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                notes += "Incorrect phospho count - too many."
                                notes_match += "Incorrect phospho count - too many. "
                            if data_row[0].count('[Phos') < (
                                    (dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                notes_match += "Incorrect phospho count - too few. "
                                notes += "Incorrect phospho count - too few. "
                            if PTM_pos_query_list[1:] in dict3[protein_query][peptide_query][pool].split(";"):
                                SP_pool_temp += "TRUE"
                                sp_pool_match_temp += "TRUE"
                                key = dict3[protein_query][peptide_query][pool]
                                key_match = dict3[protein_query][peptide_query][pool]
                            else:
                                SP_pool_temp += "FALSE"
                                sp_pool_match_temp += "FALSE"
                                pool_temp2 = "FALSE"
                                notes += "Incorrect phospho position. "
                                notes_match += "Incorrect phospho position. "
                                key = dict3[protein_query][peptide_query][pool]
                                key_match = dict3[protein_query][peptide_query][pool]
                        elif input == "Site_confident_PTM_unique.csv" and not isinstance(
                                dict3[protein_query][peptide_query][pool], str) and full_pep == "TRUE":
                            notes += "Incorrect pool - not present. "
                            notes_match += "Incorrect pool - not present. "

                        else:
                            if dict3[protein_query][peptide_query][pool] == PTM_pos_query_list[1:]:
                                SP_pool_temp += "TRUE"
                                sp_pool_match_temp += "TRUE"
                                key = dict3[protein_query][peptide_query][pool]
                                key_match = dict3[protein_query][peptide_query][pool]
                            if dict3[protein_query][peptide_query][pool] != PTM_pos_query_list[1:]:
                                SP_pool_temp += "FALSE"
                                sp_pool_match_temp += "FALSE"
                                pool_temp2 = "FALSE"
                                if isinstance(dict3[protein_query][peptide_query][pool], str):
                                    key = dict3[protein_query][peptide_query][pool]
                                    key_match = dict3[protein_query][peptide_query][pool]
                                if not isinstance(dict3[protein_query][peptide_query][pool],
                                                  str) and full_pep == "TRUE":
                                    notes += "Incorrect pool - not present. "
                                    notes_match += "Incorrect pool - not present. "
                                else:
                                    if isinstance(dict3[protein_query][peptide_query][pool], str) and (
                                            (PTM_pos_query_list.count(';')) < (
                                    dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                        notes += "Incorrect phospho count - too few. "
                                        notes_match += "Incorrect phospho count - too few. "
                                        for p in PTM_pos_query_list.split(";"):
                                            for q in dict3[protein_query][peptide_query][pool].split(";"):
                                                if p != "" and q != "" and int(p) == int(q):
                                                    notes += "Too few - match. "
                                                    notes_match += "Too few - match. "
                                    elif isinstance(dict3[protein_query][peptide_query][pool], str) and (
                                            (PTM_pos_query_list.count(';')) > (
                                    dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                        notes += "Incorrect phospho count - too many. "
                                        notes_match += "Incorrect phospho count - too many. "
                                        for p in PTM_pos_query_list.split(";"):
                                            for q in dict3[protein_query][peptide_query][pool].split(";"):
                                                if p != "" and q != "" and int(p) == int(q):
                                                    notes += "Too many - match. "
                                                    notes_match += "Too many - match. "
                                    else:
                                        notes += "Incorrect phospho position. "
                                        notes_match += "Incorrect phospho position. "
                        if site_temp == "TRUE" and pool_temp2 == "FALSE":
                            notes += "Incorrect pool. "
                            notes_match += "Incorrect pool. "
                        if "A[Phosphorylation]" in data_row[0]:
                            notes += "pAla"
                            notes_match += "pAla. "
						#filter out other modifications as well as those with incorrect phospho counts and partial peptides
                        if "Phosphorylation" not in data_row[
                            0] or "Partial" in notes_match or "Incorrect phospho count" in notes_match or "Pyro" in \
                                data_row[0] or "Oxi" in data_row[0]:
                            use = "N"
                            use_match = "N"
                        if data_row[
                            20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep2:scan:01946:GRRS[Phosphorylation]PSPGNSPSGR/3" or \
                                data_row[
                                    20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep1:scan:01937:GRRS[Phospho]PSPGNSPSGR/3" or \
                                data_row[
                                    20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep1:scan:1937:GRRS[Phosphorylation]PSPGNSPSGR/3" or \
                                data_row[
                                    20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep2:scan:1946:GRRS[Phosphorylation]PSPGNSPSGR/3":
                            use = "N"
                            use_match = "N"

                        if "Site" not in input and isinstance(dict3[protein_query][peptide_query][pool], str):
                            if data_row[0].count('[Phos') > (
                                    (dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                notes_match_site += "Incorrect phospho count - too many. "
                            if data_row[0].count('[Phos') < (
                                    (dict3[protein_query][peptide_query][pool]).count(";") + 1):
                                notes_match_site += "Incorrect phospho count - too few. "

                            for w, x, y, z in zip(prot_pos_query_list[1:].split(";"), PTM_pos_query_list[1:].split(";"),
                                                  PTM_score_list[1:].split(';'), PTM_query_list[1:].split(";")):
                                sp_pool_match_temp_site = ""
                                if str(x) in dict3[protein_query][peptide_query][pool]:
                                    sp_pool_match_temp_site += "TRUE"
                                    key_match_site = dict3[protein_query][peptide_query][pool]
                                else:
                                    if isinstance(dict3[protein_query][peptide_query][pool], str):
                                        sp_pool_match_temp_site += "FALSE"
                                        pool_temp2_site = "FALSE"
                                        key_match_site = dict3[protein_query][peptide_query][pool]
                                    if not isinstance(dict3[protein_query][peptide_query][pool],
                                                      str) and "Partial peptide" not in notes_match_site:
                                        notes_match_site += "Incorrect pool - not present. "

                                if site_temp == "TRUE" and pool_temp2_site == "FALSE":
                                    notes_match_site += "Incorrect pool. "
                                if "A[Phosphorylation]" in data_row[0]:
                                    notes_match_site += "pAla. "
                                if "Phosphorylation" not in data_row[
                                    0] or "Partial" in notes_match_site or "Incorrect phospho count" in notes_match_site or "Pyro" in \
                                        data_row[0] or "Oxi" in data_row[0]:
                                    use_match_site = "N"
                                if data_row[
                                    20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep2:scan:01946:GRRS[Phosphorylation]PSPGNSPSGR/3" or \
                                        data_row[
                                            20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep1:scan:01937:GRRS[Phosphorylation]PSPGNSPSGR/3" or \
                                        data_row[
                                            20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep1:scan:1937:GRRS[Phosphorylation]PSPGNSPSGR/3" or \
                                        data_row[
                                            20] == "mzspec:PXD007058:SF_200217_pPeptideLibrary_pool5_HCDOT_rep2:scan:1946:GRRS[Phosphorylation]PSPGNSPSGR/3":
                                    use_match_site = "N"
                                if z != "":
                                    peptide_unmod_site.append(peptide_query)
                                    peptide_site.append(data_row[0])
                                    protein_site.append(a)
                                    protein_temp_site.append(a.split("_")[0])
                                    all_USI_site.append(data_row[20])
                                    score_site.append(data_row[7])
                                    PTM_site.append(z)
                                    PTM_score_site.append(y.strip("']").strip("['"))
                                    PTM_pos_site.append(x)
                                    PTM_pro_pos_site.append(w)
                                    pep_count_site.append(data_row[17])
                                    all_pools_site.append(pool)
                                    SP_match_site.append(sp_match_temp)
                                    SP_site_match_site.append(sp_site_match_temp)
                                    SP_pool_match_site.append(sp_pool_match_temp_site)
                                    all_notes_match_site.append(notes_match_site)
                                    all_key_match_site.append(key_match_site)
                                    all_use_match_site.append(use_match_site)

                        if PTM_query_list != "":
                            peptide_unmod.append(peptide_query)
                            peptide.append(data_row[0])
                            protein.append(a)
                            protein_temp.append(a.split("_")[0])
                            all_USI.append(data_row[20])
                            score.append(data_row[7])
                            PTM.append(PTM_temp)
                            PTM_highest_score.append(max(data_row[9].split(';')))
                            PTM_score.append(re.sub("[][']", '', PTM_score_list)[1:])
                            PTM_pos.append(re.sub("[][']", '', PTM_pos_query_list)[1:])
                            PTM_pro_pos.append(prot_pos_query_list[1:])
                            pep_count.append(data_row[17])
                            all_pools.append(pool)
                            SP_match.append(sp_match_temp)
                            SP_site_match.append(sp_site_match_temp)
                            SP_pool_match.append(sp_pool_match_temp)
                            all_notes_match.append(notes_match)
                            all_key_match.append(key_match)
                            all_use_match.append(use_match)
                SP_site.append(SP_site_temp)
                SP_pep.append(SP_pep_temp)
                SP_pool.append(SP_pool_temp)
                all_notes.append(notes)
                all_key.append(key)
                all_use.append(use)
		#Syntehtic peptides with pool
        df = pd.read_csv(input)
        df['Synthetic_Pep_site_match'] = SP_site
        df['Synthetic_Pep_match'] = SP_pep
        df['Synthetic_Pep_Pool-site_match'] = SP_pool
        df['Notes'] = all_notes
        df2 = pd.DataFrame(
            {"Peptide_mod": peptide, "Peptide": peptide_unmod, "Protein": protein, "Protein_temp": protein_temp,
             "Pool": all_pools, "All_USI": all_USI, "Score": score, "PTM": PTM, "PTM_score": PTM_score,
             "PTM_best_score": PTM_highest_score, "PTM_positions": PTM_pos, "PTM_Protein_Positions": PTM_pro_pos,
             "PSM_count": pep_count, "Synthetic_peptide_match": SP_match, "Synthetic_peptide_site_match": SP_site_match,
             "Synthetic_peptide_pool-site_match": SP_pool_match, "Notes": all_notes_match, "Key": all_key_match,
             "Use_for_SP_analysis": all_use_match})
        df2['Temp'] = df2['Protein_temp'] + ";" + df2['Peptide_mod'] + ";" + df2['PTM_score'].astype(str) + ";" + df2[
            'Pool'].astype(str)
        df2 = df2.drop_duplicates(subset=('Temp'), keep='last', inplace=False)
        file2 = input.replace(".csv", "_Synthetic_Peptides_w-pool.csv")
		#Synthetic peptides - site based
        df2 = df2.drop(['Protein_temp', 'Temp'], axis=1)
        df2 = df2.reset_index(drop=True)
        df2.to_csv(file2, index=False)
        df3 = pd.DataFrame({"Peptide_mod": peptide_site, "Peptide": peptide_unmod_site, "Protein": protein_site,
                            "Protein_temp": protein_temp_site, "Pool": all_pools_site, "All_USI": all_USI_site,
                            "Score": score_site, "PTM": PTM_site, "PTM_score": PTM_score_site,
                            "PTM_best_score": PTM_score_site, "PTM_positions": PTM_pos_site,
                            "PTM_Protein_Positions": PTM_pro_pos_site, "PSM_count": pep_count_site,
                            "Synthetic_peptide_match": SP_match_site,
                            "Synthetic_peptide_site_match": SP_site_match_site,
                            "Synthetic_peptide_pool-site_match": SP_pool_match_site, "Notes": all_notes_match_site,
                            "Key": all_key_match_site, "Use_for_SP_analysis": all_use_match_site})
        df3['Temp'] = df3['Protein_temp'] + ";" + df3['Peptide_mod'] + ";" + df3['PTM_score'].astype(str) + ";" + df3[
            'Pool'].astype(str)
        df3 = df3.drop_duplicates(subset=('Temp'), keep='last', inplace=False)
        file3 = input.replace(".csv", "_Synthetic_Peptides_w-pool_Site-based.csv")
        df3 = df3.drop(['Protein_temp', 'Temp'], axis=1)
        df3 = df3.reset_index(drop=True)
        df3.to_csv(file3, index=False)
        return df2
	
	#make dictionary of synthetic peptides
    dict=makehash() #protein-seq-position
    dict2=makehash() #protein-seq
    dict3=makehash() #protein-seq-pool-pos
    synthetic_pep_pro_list=[]
	#Loop through fasta file of synthetic peptide key and create dictionaries for look up
    for seq_record in SeqIO.parse(s_pep, "fasta"):
        seq = str(seq_record.seq)
        pep_temp = re.sub(r'(?<=[RK])(?=[^P])', '\n', seq, 0, re.DOTALL)
        pep = (pep_temp.split())
        protein = seq_record.id

        synthetic_pep_pro_list.append(protein)
        dict2[protein][seq]="Synthetic_Peptide"
        description = seq_record.description
        pos_temp=description.split("|")[1].strip(" ")
        pos_temp=re.sub('[S|T|Y|V|x|E]', '', pos_temp)

        if "or" in pos_temp:
            positions=pos_temp.split("or")
            for p in positions:
                pos=p.strip(" ")
                dict[protein][seq][pos]="Synthetic_Peptide"
                tryptic(protein,seq,pos)
        else:
            pos_temp=pos_temp.replace("&",";")
            dict[protein][seq][pos_temp]="Synthetic_Peptide"
            tryptic(protein,seq,pos_temp)

        if len(pep[0]) != len(seq):
            search_pos=pos_temp.split("or")
            for p in search_pos:
                if p!="" and "&" not in p and ";" not in p:
                    counter = 0
                    while counter<len(pep):
                        missed=1
                        while missed<=3 and missed<=len(pep[counter:]):
                            missed_temp=pep[counter:counter+missed]
                            missed_pep="".join(str(x) for x in missed_temp)
                            if int(p)<=len(missed_pep) and int(p)>0:
                                dict[protein][missed_pep][str(p)]="Synthetic_Peptide"
                                tryptic(protein,missed_pep,str(p))
                            missed+=1
                        if int(p)<=len(pep[counter])and int(p)>0:
                            dict[protein][pep[counter]][str(p)]="Synthetic_Peptide"
                            tryptic(protein,pep[counter],str(p))
                        p=int(p)-len(pep[counter])
                        counter+=1
                else:
                    search_pos2=re.split(r';|&', p)
                    counter=0
                    if len(search_pos2)==2:
                        s = int(search_pos2[0])
                        q = int(search_pos2[1])
                        while counter<len(pep) and s>=1 and q>=1:
                            missed=1
                            while missed<=4 and missed<=len(pep[counter:]):
                                missed_temp=pep[counter:counter+missed]
                                missed_pep="".join(str(x) for x in missed_temp)
                                if s<=len(missed_pep) and q<=len(missed_pep):
                                    temp=str(s)+";"+str(q)
                                    dict[protein][missed_pep][temp]="Synthetic_Peptide"
                                    tryptic(protein,missed_pep,temp)
                                missed+=1
                            if s<=len(pep[counter]) and q<=len(pep[counter]):
                                temp=str(s)+";"+str(q)
                                dict[protein][pep[counter]][temp]="Synthetic_Peptide"
                                tryptic(protein,pep[counter],temp)
                            else:
                                s=s-len(pep[counter])
                                q=q-len(pep[counter])
                            counter+=1
    synthetic_pep_pro_list2=[]
    for seq_record in SeqIO.parse(s_pep_pools, "fasta"):
        seq = str(seq_record.seq)
        pep_temp = re.sub(r'(?<=[RK])(?=[^P])', '\n', seq, 0, re.DOTALL)
        pep = (pep_temp.split())
        protein = seq_record.id
        synthetic_pep_pro_list2.append(protein)
        description = seq_record.description
        pos_temp=description.split("|")[2].strip(" ")
        pos_temp=re.sub('[S|T|Y|V|x|E]', '', pos_temp)
        pool=int(description.split("|")[1].strip(" "))

        if "or" in pos_temp:
            positions=pos_temp.split("or")
            for p in positions:
                pos=p.strip(" ")
                dict3[protein][seq][pool]=pos
                tryptic_pool(protein,seq,pool,pos)
        else:
            pos_temp=pos_temp.replace("&",";")
            dict3[protein][seq][pool]=pos_temp
            tryptic_pool(protein,seq,pool,pos_temp)

        if len(pep[0]) != len(seq):
            search_pos=pos_temp.split("or")
            for p in search_pos:
                if p!="" and "&" not in p and ";" not in p:
                    counter = 0
                    while counter<len(pep):
                        missed=1
                        while missed<=3 and missed<=len(pep[counter:]):
                            missed_temp=pep[counter:counter+missed]
                            missed_pep="".join(str(x) for x in missed_temp)
                            if int(p)<=len(missed_pep) and int(p)>0:
                                dict3[protein][missed_pep][pool]=str(p)
                                tryptic_pool(protein,missed_pep,pool,str(p))
                            missed+=1
                        if int(p)<=len(pep[counter])and int(p)>0:
                            dict3[protein][pep[counter]][pool]=str(p)
                            tryptic_pool(protein,pep[counter],pool,str(p))
                        p=int(p)-len(pep[counter])
                        counter+=1
                else:
                    search_pos2=re.split(r';|&', p)
                    counter=0
                    if len(search_pos2)==2:
                        s = int(search_pos2[0])
                        q = int(search_pos2[1])
                        while counter<len(pep) and s>=1 and q>=1:
                            missed=1
                            while missed<=4 and missed<=len(pep[counter:]):
                                missed_temp=pep[counter:counter+missed]
                                missed_pep="".join(str(x) for x in missed_temp)
                                if s<=len(missed_pep) and q<=len(missed_pep):
                                    temp=str(s)+";"+str(q)
                                    dict3[protein][missed_pep][pool]=temp
                                    tryptic_pool(protein,missed_pep,pool,temp)
                                missed+=1
                            if s<=len(pep[counter]) and q<=len(pep[counter]):
                                temp=str(s)+";"+str(q)
                                dict3[protein][pep[counter]][pool]=temp
                                tryptic_pool(protein,pep[counter],pool,temp)
                            else:
                                s=s-len(pep[counter])
                                q=q-len(pep[counter])
                            counter+=1

    for i in input_files:
        file=i+"/"+sub+"/All_confident_PTM_no_collapse.csv"
        df(i,sub,file)