import os
import pandas as pd
import csv
import re
import time

#Spectrum comparison between search files
def per_spectrum(wd,file_list,file_name_list,output):
    counter=1
    os.chdir(wd)
    search_spectrum = []
    dict = {}
    header = "Spectrum,"
    for f, n in zip(file_list, file_name_list):
        dict[f]={}
        reader = csv.reader(open(f, 'r'))
        counter=1
        for row in reader:
            counter+=1
            if row[0] == "Peptide_mod":
                list = [0, 2, 4, 6, 7, 8, 9, 10, 11]
                for i in list:
                    header += row[i] + "_" + n + ","
            else:
                spectrum = row[3].replace(".mzML", "").replace(".0", ".")
                dict[f][spectrum] = row
                if spectrum not in search_spectrum:
                    search_spectrum.append(spectrum)
    output_file = open(output, 'w')
    output_file.write(header + "\n")
    output_file.close()

    output_file = open(output, 'a')
    for s in search_spectrum:
        line = s + ","
        for f in file_list:
            if s in dict[f]:
                temp = dict[f][s]
                temp_all = temp[0] + "," + temp[2] + "," + temp[4] + "," + temp[6] + "," + temp[7] + "," + temp[8] + "," + temp[
                    9] + "," + temp[10] + "," + temp[11] + ","
                line += temp_all
            else:
                line += ',' + ',' + ',' + ',' + ',' + ',' + ',' + ',' + ','
        output_file.write(line + "\n")
        counter += 1
    output_file.close()

#Create csv of matching spectrum 
def matches(input, file_list, final_output):
    df = pd.read_csv(input, sep=",")
    for i in file_list:
        df['PTMs_'+i] = df['PTMs_'+i].astype(str)
        df['Peptide_mod_' + i] = df['Peptide_mod_' + i].astype(str)
    matches = []
    pep_matches = []
    mass_matches = []
    for a in range(len(df)):
        peptides=[]
        pep_unmod=[]
        pep_mass=[]
        mass=0
        for i in file_list:
            peptides.append(df.loc[a,'Peptide_mod_'+i])
            pep_temp= re.sub("[\[].*?[\]]", "",df.loc[a,'Peptide_mod_'+i])
            pep_unmod.append(pep_temp)
            if df.loc[a,'PTMs_'+i]=="nan":
                PTM_list=""
            else:
                PTM_list=df.loc[a,'PTMs_'+i]
            for p in PTM_list.split(";"):
                if "Phosphorylation" in p:
                    mass += 80
                if "Pyrophosphorylation" in p:
                    mass += 160
                if "Oxidation" in p:
                    mass += 16
            pep_mass.append(mass)
        if len(set(peptides)) == 1:
            matches.append("TRUE")
        else:
            matches.append("FALSE")
        if len(set(pep_unmod)) == 1:
            pep_matches.append("TRUE")
        else:
            pep_matches.append("FALSE")
        if len(set(pep_mass)) == 1:
            mass_matches.append("TRUE")
        else:
            mass_matches.append("FALSE")
    df['Match'] = matches
    df['Pep_Match'] = pep_matches
    df['Mass_shift_Match'] = mass_matches
    df.to_csv(final_output, index=False)


def comparison(working,sub,file_list,file_name_list):
    start_time = time.time()
    os.chdir(working)    
    if not os.path.exists("Spectrum_Comparisons"):
        os.mkdir("Spectrum_Comparisons")
    file_list_output=""
    for i in file_name_list:
        file_list_output+=i+"_"
    spectrum_file="Spectrum_Comparisons/"+file_list_output[:-1]
    if not os.path.exists(spectrum_file):
        os.mkdir(spectrum_file)
    output=working+"/"+spectrum_file+"/"+sub+"_spectrum_comparison_"+file_list_output[:-1]+".txt"
    per_spectrum(working,file_list,file_name_list,output)
    #Seperate matches
    final=working+"/"+spectrum_file+"/"+sub+"_spectrum_comparison_"+file_list_output[:-1]+".csv"
    matches(output,file_name_list,final)