import xml.sax
import csv
import pandas as pd
import os
import re
import collections
import time

start_time = time.time()

#convert xml output file from PTM prophet to csv
def convert(input):
    class TPPHandler(xml.sax.ContentHandler):
        def __init__(self):
            self.CurrentData = ""

        # Call when an element starts
        def startElement(self, tag, attributes):
            if tag == "spectrum_query":
                spectrum = attributes["spectrum"]
                output.write("\n" + spectrum + ",")
                rt = attributes["retention_time_sec"]
                output.write(rt + ",")
                z=attributes["assumed_charge"]
                output.write(z +",")
            if tag == "search_hit":
                no_proteins = int(attributes["num_tot_proteins"])
                output.write(str(no_proteins) + ",")
                protein = attributes["protein"]
                output.write(protein + ",")
                peptide = attributes["peptide"]
                output.write(peptide + ",")
                matched_ions = attributes["num_matched_ions"]
                total_ions = attributes["tot_num_ions"]
                ions = matched_ions + "/" + total_ions
                output.write(ions + ",")
                mass = attributes["calc_neutral_pep_mass"]
                output.write(mass + ",")
            if tag == "modification_info":
                mod_info = attributes["modified_peptide"]
                output.write("," + mod_info + ",")

            if tag == "search_score":
                if attributes["name"] == "expect":
                    expect = attributes["value"]
                    output.write(expect + ",")
            if tag == "alternative_protein":
                protein_alt = attributes["protein"]
                output.write(protein_alt + ";")
            if tag == "peptideprophet_result":
                prob = attributes["probability"]
                output.write(prob + ",")
            if tag == "ptmprophet_result":
                ptm_pep = attributes["ptm_peptide"]
                output.write(ptm_pep + ",")
                ptm = attributes["ptm"]
                output.write(ptm + ",")
                mod_mass = float(attributes["ptm"].split(":")[1])
                output.write(str(mod_mass) + ",")


    output = open("XML2CSV_temp.csv", "w")
    output.close()
    output = open("XML2CSV_temp.csv", "a")

    # create an XMLReader
    parser = xml.sax.make_parser()
    Handler = TPPHandler()
    parser.setContentHandler(Handler)
    parser.parse(input)
    output.close()
    print("xml parsing done--- %s seconds ---" )

    prob_all = []
    ptm_pep_all = []
    peptide_mod_all = []
    spectrum_all = []
    expect_all = []
    ions_all = []
    peptide_all = []
    protein_all = []
    no_proteins_all = []
    mass_all = []
    rt_all = []
    z_all=[]
    PTM_all = []

    file = open("XML2CSV_temp.csv", 'r')
    #print("file open--- %s seconds ---" % (time.time() - start_time))
    for line in file:
        if "," in line:
            data_row = line.strip().split(',')
            spectrum=data_row[0]
            rt=data_row[1]
            z=data_row[2]
            no_proteins=int(data_row[3])
            if no_proteins==1:
                proteins=data_row[4]
            else:
                proteins=data_row[4]
                for i in data_row[8].split(";")[:no_proteins-1]:
                    proteins+=":"+i
            peptide=data_row[5]
            ions=str(data_row[6])
            mass=data_row[7]
            regex=re.compile(r'[a-zA-Z]')
            if regex.search(data_row[9]):
                mod_peptide = data_row[9]
                expect = data_row[10]
                prob = data_row[11]
            else:
                mod_peptide=peptide
                expect=0
                prob=data_row[9]
            ptm_pep_query=12
            ptm_pep=""
            while len(data_row)>ptm_pep_query:
                if "(" in data_row[ptm_pep_query]:
                    ptm_pep+=":"+data_row[ptm_pep_query]
                ptm_pep_query+=3
            if ptm_pep=="":
                ptm_pep+="--"

            PTM=""
            if "S[247]" in mod_peptide or "T[261]" in mod_peptide or "Y[323]" in mod_peptide:
                PTM+=";Pyrophosphorylation"
            if "S[167]" in mod_peptide or "T[181]" in mod_peptide or "Y[243]" in mod_peptide or "A[151]" in mod_peptide or "G[137]" in mod_peptide or "H[217]" in mod_peptide or "D[195]" in mod_peptide or \
                    "E[209]" in mod_peptide or "K[208]" in mod_peptide or "R[236]" in mod_peptide or "P[177]" in mod_peptide:
                PTM+=";Phosphorylation"
            if "M[147]" in mod_peptide or "W[202]" in mod_peptide or "P[113]" in mod_peptide:
               PTM+=";Oxidation"
           #mod_peptide=mod_peptide.replace("[","(").replace("]",")")
            mod_peptide=mod_peptide.replace("S[247]","S[Pyrophosphorylation]")
            mod_peptide = mod_peptide.replace("T[261]", "T[Pyrophosphorylation]")
            mod_peptide=mod_peptide.replace("Y[323]","Y[Pyrophosphorylation]")
            mod_peptide=mod_peptide.replace("S[167]","S[Phosphorylation]")
            mod_peptide = mod_peptide.replace("T[181]", "T[Phosphorylation]")
            mod_peptide=mod_peptide.replace("Y[243]","Y[Phosphorylation]")
            mod_peptide=mod_peptide.replace("A[151]","A[Phosphorylation]")
            mod_peptide = mod_peptide.replace("G[137]", "G[Phosphorylation]")
            mod_peptide = mod_peptide.replace("L[193]", "L[Phosphorylation]")
            mod_peptide = mod_peptide.replace("D[195]", "D[Phosphorylation]")
            mod_peptide = mod_peptide.replace("E[209]", "E[Phosphorylation]")
            mod_peptide = mod_peptide.replace("P[177]","P[Phosphorylation]")

            mod_peptide=mod_peptide.replace("M[147]","M[Oxidation]")
            mod_peptide=mod_peptide.replace("W[202]","W[Oxidation]")
            mod_peptide=mod_peptide.replace("P[113]","P[Oxidation]")
            mod_peptide=mod_peptide.replace("N[115]","N[Deamination]")
            mod_peptide = mod_peptide.replace("Q[129]", "Q[Deamination]")
            mod_peptide = mod_peptide.replace("[111]", "[Pyro_glu]")
            mod_peptide = mod_peptide.replace("[143]", "[Pyro_glu]")
            mod_peptide = mod_peptide.replace("n[43]", "[Acetyl]-")

            prob_all.append(prob)
            ptm_pep_all.append(ptm_pep[1:])
            peptide_mod_all.append(mod_peptide)
            spectrum_all.append(spectrum)
            expect_all.append(expect)
            ions_all.append(ions)
            peptide_all.append(peptide)
            protein_all.append(proteins)
            no_proteins_all.append(no_proteins)
            mass_all.append(mass)
            rt_all.append(rt)
            PTM_all.append(PTM[1:])
            z_all.append(z)
    #print("make df--- %s seconds ---" % (time.time() - start_time))
    df = pd.DataFrame({"probability": prob_all, "ptm_peptide": ptm_pep_all, "spectrum": spectrum_all, "expect": expect_all,
         "ions": ions_all, "peptide": peptide_all, "protein": protein_all, "num_prots": no_proteins_all,
         "calc_neutral_pep_mass": mass_all, "RT": rt_all, "PTM": PTM_all, "mod_peptide": peptide_mod_all, "charge":z_all})
    df.to_csv("interact.ptm.pep.csv", index=False)