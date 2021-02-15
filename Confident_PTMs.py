import pandas as pd
import Extract_peptides_DB
import re
import os
import collections
import time

start_time = time.time()

def makehash():
    return collections.defaultdict(makehash)

#Filter for FDR (and A score) cutoff
def confident(input,fasta_input,FDR,A,final_output):
    # set params
    FDR_cutoff = FDR
    Ascore_cutoff = A

    peptides=[]
    peptide_mod_all=[]
    proteins=[]
    peptide_unmod=[]
    spectrum=[]
    USI=[]
    source=[]
    score=[]
    PTMs=[]
    PTM_scores=[]
    filtered_PTM_list=[]
    filtered_PTM_score_list=[]
    filtered_PTM_position_list=[]
    filtered_PTM_protein_position_list=[]
    PTM_info=[]
    FDR=[]
    q_all=[]

    # import or create dictionary for protein position
    dict_name=fasta_input.replace(".fasta","_Peptides_DB.csv")
    if not os.path.isfile(dict_name):
        print("Writing peptide dict")
        dict = Extract_peptides_DB.extract_peptide_positions(fasta_input)
        print("Done")
    else:
        print("Peptide dict already exists")
        dict = makehash()
        file = open(dict_name, "r")
        for line in file:
            data_row = line.strip().split(' ')
            dict[data_row[0]][data_row[1]] = data_row[2]

    #print("DICT done --- %s seconds ---" % (time.time() - start_time))
    # extract values
    file = open(input, 'r')
    counter = 1
    for line in file:
        #ignore header
        if counter==1:
            counter+=1
        else:
            data_row = line.strip().split(',')
            q = float(data_row[16])
            #filter for FDR
            if q <= FDR_cutoff:
                peptide = (data_row[1])
                peptide_mod_all.append(data_row[0])
                peptides.append(peptide)
                proteins.append(data_row[2].replace("sp|",""))
                source.append(data_row[3])
                score.append(data_row[4])
                PTM = (data_row[5])
                PTMs.append(PTM)
                PTM_position = (data_row[6])
                PTM_score = (data_row[7])
                PTM_scores.append(PTM_score)
                PTM_info.append(data_row[8])
                spectrum.append(data_row[9])
                USI.append(data_row[10])
                FDR.append(data_row[15])
                q_all.append(q)

                positions = ""
                filtered_PTMs = ""
                filtered_scores = ""
                filtered_positions = ""

                #remove mods from peptide
                peptide_temp = re.split(':',peptide)
                peptide_edit = re.sub('[()+0-9. "]', '', peptide_temp[0])
                peptide_unmod.append(peptide_edit)

                if PTM != "-" and PTM!="":
                    for a, b, c in zip(PTM.split(";"), PTM_score.split(";"), PTM_position.split(";")):
                        if Ascore_cutoff <= float(b):
                            filtered_PTMs += ";" + a
                            filtered_scores += ";" + b
                            filtered_positions += ";" + c

                protein_list_temp = (data_row[2])
                protein_list = protein_list_temp.split(':')
                if PTM != "-" and PTM!="":
                    for p in protein_list:
                        for a, b, c in zip(PTM.split(";"), PTM_score.split(";"), PTM_position.split(";")):
                            if int(c)<=len(peptide_edit):
                                PTM_position_query = c
                                if Ascore_cutoff <= float(b):
                                    protein = p
                                    # add pos to pep start pos = modification pos on whole protein
                                    if protein in dict:
                                        query=str(dict[protein][peptide_edit])
                                        if query.isnumeric():
                                            PTM_pro_pos = ";"+str(int(query)+int(PTM_position_query))
                                        else:
                                            PTM_pro_pos=";0"
                                    # search protein accession with added prefix
                                    else:
                                        if protein != "":
                                            protein_pre = "sp|" + protein
                                            query = str(dict[protein_pre][peptide_edit])
                                            if query.isnumeric():
                                                PTM_pro_pos = ";" + str(int(query) + int(PTM_position_query))
                                            else:
                                                PTM_pro_pos = ";0"
                                    if positions==""or positions[-1:]==":":
                                        positions+=PTM_pro_pos[1:]
                                    else:
                                        positions+=PTM_pro_pos
                                else:
                                    filtered_PTMs+=""
                                    filtered_scores += ""
                                    filtered_positions += ""
                                    positions+=""
                        positions+=":"
                else:
                    filtered_PTMs += ""
                    filtered_scores += ""
                    filtered_positions += ""
                    positions += ""

                filtered_PTM_list.append(filtered_PTMs[1:])
                filtered_PTM_score_list.append(filtered_scores[1:])
                filtered_PTM_position_list.append(filtered_positions[1:])
                filtered_PTM_protein_position_list.append(str(positions[:-1]))

    df = pd.DataFrame({"Peptide_mod": peptide_mod_all, "Peptide":peptide_unmod,"Protein": proteins, "Spectrum": spectrum, "USI": USI, "Source": source, "Score": score,
                               "PTMs":filtered_PTM_list, "PTM_scores_all":filtered_PTM_score_list, "PTM_positions":filtered_PTM_position_list,
                                "PTM_Protein_Positions":filtered_PTM_protein_position_list, "PTM_info": PTM_info, "FDR": FDR, "Qvalue": q_all})
    df.to_csv(final_output,index=False)