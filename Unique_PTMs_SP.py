import pandas as pd
import re
import numpy as np
import time

start_time = time.time()

#Collapse for best scoring for each peptide, protein site, mass shift on protein or no collapse
def unique(output, peptide_based,site_based,mass_based,non_collapse):
    file = open(output,'r')

    peptide=[]
    peptide_unmod=[]
    protein=[]
    spectrum=[]
    USI=[]
    source=[]
    score=[]
    PTMs=[]
    PTM_score=[]
    PTM_positions=[]
    PTM_protein_positions=[]
    PTM_info=[]
    FDR=[]
    q=[]
    mass_shift=[]
    peptide_split = []
    peptide_unmod_split = []
    protein_split = []
    spectrum_split = []
    USI_split = []
    source_split = []
    score_split = []
    PTMs_split = []
    PTM_score_split = []
    PTM_positions_split = []
    PTM_protein_positions_split = []
    PTM_info_split = []
    FDR_split = []
    q_split = []
    mass_shift_split=[]
    pool=[]
    pool_split = []
    for line in file:
        data_row = line.strip().split(',')
        peptide.append(data_row[0])
        peptide_unmod.append(data_row[1])
        protein_list=(data_row[2])
        protein.append(protein_list)
        spectrum.append(data_row[3])
        USI.append(data_row[4])
        source.append(data_row[5])
        score.append(data_row[6])
        PTM_list=(data_row[7])
        PTMs.append(PTM_list)
        PTM_score_list=(data_row[8])
        PTM_score.append(PTM_score_list)
        PTM_positions_list=(data_row[9])
        PTM_positions.append(PTM_positions_list)
        pro_pos=(data_row[10])
        PTM_protein_positions.append(pro_pos)
        PTM_info.append(data_row[11])
        FDR.append(data_row[12])
        q.append(data_row[13])

        if "_" in data_row[5]:
            pool.append(data_row[5].split("_")[3])
        else:
            pool.append("")

        mass_temp=0
        for p in PTM_list.split(";"):
            if "Phosphorylation" in p:
                mass_temp+=80
            if "Pyrophosphorylation" in p:
                mass_temp+=160
            if "Oxidation" in p:
                mass_temp+=16
        mass_shift.append(str(mass_temp))

        for a,b in zip(protein_list.split(":"),pro_pos.split(":")):
            query_protein=a
            query_pro_pos_list=b
            for w,x,y,z in zip(query_pro_pos_list.split(";"),PTM_list.split(";"),PTM_score_list.split(";"),PTM_positions_list.split(";")):
                protein_split.append(query_protein)
                PTMs_split.append(x)
                PTM_protein_positions_split.append(w)
                PTM_score_split.append(y)
                PTM_positions_split.append(z)

                peptide_split.append(data_row[0])
                peptide_unmod_split.append(data_row[1])
                spectrum_split.append(data_row[3])
                USI_split.append(data_row[4])
                source_split.append(data_row[5])
                score_split.append(data_row[6])
                PTM_list=data_row[7]
                PTM_info_split.append(data_row[11])
                FDR_split.append(data_row[12])
                q_split.append(data_row[13])

                if "_" in data_row[5]:
                    pool_split.append(data_row[5].split("_")[3])
                else:
                    pool_split.append("")

                mass_split_temp = 0
                for p in PTM_list.split(";"):
                    if "Phosphorylation" in p:
                        mass_split_temp += 80
                    if "Pyrophosphorylation" in p:
                        mass_split_temp += 160
                    if "Oxidation" in p:
                        mass_split_temp += 16
                mass_shift_split.append(str(mass_split_temp))
				
	#Collapse by mass shift
    df = pd.DataFrame({"Peptide_mod": peptide[1:], "Peptide": peptide_unmod[1:], "Protein": protein[1:], "Spectrum": spectrum[1:], "Pool": pool[1:], "USI": USI[1:], "Source": source[1:], "Score": score[1:], "PTM": PTMs[1:],
                       "PTM_score": PTM_score[1:], "PTM_positions": PTM_positions[1:],"PTM_Protein_Positions": PTM_protein_positions[1:], "PTM_info": PTM_info[1:], "Mass shift": mass_shift[1:], "FDR": FDR[1:], "Qvalue": q[1:]})
    df=df.sort_values(['Peptide','Pool','Mass shift','Score','PTM_score'],ascending=[True,True,True,True,True])
    df['Pep-pool'] = df['Peptide'] +";" + df['Pool']+";"+df['Mass shift']
    df['PSM_count'] = df.groupby(['Pep-pool'], sort=False).cumcount() + 1
    USI_list=[]
    for i in df['USI'].groupby(df['Pep-pool']).apply(';'.join):
        USI_list.append(i)
    protein_list=[]
    for i in df['Protein'].groupby(df['Pep-pool']).apply(','.join):
        protein_list_temp = []
        all_protein = ""
        for p in re.split(",|:",i):
            if p not in protein_list_temp:
                protein_list_temp.append(p)
                all_protein+=";"+p
        protein_list.append(all_protein[1:])
    df['Unique'] = np.where(df['PSM_count'] == 1, "TRUE", "FALSE")
    df=df.drop_duplicates(subset=('Pep-pool'),keep='last', inplace=False)
    df['All_Proteins']=protein_list
    df['All_USI']=USI_list
    df.to_csv(mass_based,index=False)
	
	#Collapse by peptide
    df2 = pd.DataFrame({"Peptide_mod": peptide[1:], "Peptide": peptide_unmod[1:], "Protein": protein[1:], "Spectrum": spectrum[1:],"Pool": pool[1:],"USI": USI[1:], "Source": source[1:], "Score": score[1:], "PTM": PTMs[1:],
                        "PTM_score": PTM_score[1:],"PTM_positions": PTM_positions[1:],"PTM_Protein_Positions": PTM_protein_positions[1:], "PTM_info": PTM_info[1:], "Mass shift": mass_shift[1:], "FDR": FDR[1:], "Qvalue": q[1:]})
    df2=df2.sort_values(['Pool','Peptide','Score','PTM_score'],ascending=[True,True,True,True])
    df2['Pro-pep'] =+df2['Pool'] + "-" + df2['Peptide']
    df2['Peptide_count'] = df2.groupby(['Pro-pep'], sort=False).cumcount()+1

    USI_list2 = []
    for i in df2['USI'].groupby(df2['Pro-pep']).apply(';'.join):
        USI_list2.append(i)
    protein_list2 = []
    for i in df2['Protein'].groupby(df2['Pro-pep']).apply(','.join):
        protein_list_temp2 = []
        all_protein2 = ""
        for p in re.split(",|:",i):
            if p not in protein_list_temp2:
                protein_list_temp2.append(p)
                all_protein2 += ";" + p
        protein_list2.append(all_protein2[1:])

    df2['Unique'] = np.where(df2['Peptide_count'] == 1, "TRUE", "FALSE")
    df2=df2.drop_duplicates(subset=('Pro-pep'),keep='last', inplace=False)
    df2['All_Proteins'] = protein_list2
    df2['All_USI']=USI_list2
    df2.to_csv(peptide_based,index=False)
	
	#Collapse by protein position
    df3 = pd.DataFrame({"Peptide_mod": peptide_split[1:], "Peptide":peptide_unmod_split[1:], "Protein": protein_split[1:],"Spectrum":spectrum_split[1:], "Pool":pool_split[1:],"USI":USI_split[1:], "Source":source_split[1:], "Score": score_split[1:], "PTM": PTMs_split[1:],
                        "PTM_score":PTM_score_split[1:],"PTM_positions":PTM_positions_split[1:],"PTM_Protein_Position":PTM_protein_positions_split[1:],"PTM_info": PTM_info_split[1:], "Mass shift":mass_shift_split[1:],"FDR": FDR_split[1:], "Qvalue": q_split[1:]})
    df3=df3.sort_values(['Protein','Pool','PTM_Protein_Position','Score','PTM_score'],ascending=[True,True,True,True,True])
    df3['Pro-pos'] = df3['Protein'] + "-" +df3['Pool'] + "-" + df3['PTM_Protein_Position']
    df3['Protein_count'] = df3.groupby(['Pro-pos'], sort=False).cumcount()+1
    USI_list3 = []
    for i in df3['USI'].groupby(df3['Pro-pos']).apply(';'.join):
        USI_list3.append(i)

    df3['Unique'] = np.where(df3['Protein_count'] == 1, "TRUE", "FALSE")
    #print("DF3 drop--- %s seconds ---" % (time.time() - start_time))
    df3=df3.drop_duplicates(subset=('Pro-pos'),keep='last', inplace=False)
    df3['All_Proteins']=df3['Protein']
    df3['All_USI']=USI_list3
    df3.to_csv(site_based,index=False)
	
	#No collapse
    df4 = pd.DataFrame({"Peptide_mod": peptide[1:], "Peptide": peptide_unmod[1:], "Protein": protein[1:], "Spectrum": spectrum[1:],"Pool": pool[1:],"USI": USI[1:], "Source": source[1:], "Score": score[1:], "PTM": PTMs[1:],"PTM_score": PTM_score[1:],
                        "PTM_positions": PTM_positions[1:],"PTM_Protein_Positions": PTM_protein_positions[1:], "PTM_info": PTM_info[1:], "Mass shift": mass_shift[1:],"FDR": FDR[1:], "Qvalue": q[1:]})
    df4 = df4.sort_values(['Pool','Peptide', 'Score', 'PTM_score'], ascending=[True, True, True,True])
    df4["x"] = "-"
    df4["y"] = "-"
    df4['z'] = "-"
    df4['All_Proteins']=df4['Protein']
    df4['All_USI']=df4['USI']
    df4.to_csv(non_collapse,index=False)


