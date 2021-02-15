import os
from csv import writer
from csv import reader
from pyteomics import mgf, auxiliary
import pandas as pd
import collections

#model PSM probabilities from "map_p_to_probs.csv"
def value(x1,x2,y1,y2,query):
    m=(y2-y1)/(x2-x1)
    c=y1-(x1*m)
    y=(m*query)+c
    return(y)

def makehash():
    return collections.defaultdict(makehash)

#generates scan number from Peaks output
def scan(search_files_location, input, output, working):

    df = pd.read_csv(input)
    df.rename(columns={'Source File': 'Source'}, inplace=True)
    Search_files = df.Source.unique()

    os.chdir(search_files_location)
    #dict={'file':{'mass':{'time':'spectrum'}}}
    dict=makehash()
    for file in Search_files:
        with mgf.read(file) as reader:
            for entry in reader:
                spectrum_temp=(entry['params']['title'])
                spectrum_target = spectrum_temp.split( )[0]
                mass_target=round(((entry['params']['pepmass'])[0]))
                time_target=round(float(entry['params']['rtinseconds']))
                dict[file][mass_target][time_target]=spectrum_target

    os.chdir(working)
    counter=1
    file1 = open(input, 'r')
    spectrum_list=[]
    prob_value_list=[]
    prob = pd.read_csv("map_p_to_probs.csv")
    prob=prob.dropna()
    score_list=prob['score_mid_points'].to_list()
    prob_list = prob['probabilities'].tolist()
    for line in file1:
        if counter == 1:
            counter += 1
        else:
            data_row = line.strip().split(',')
            target_file = data_row[13]
            mass = round(float(data_row[5]))
            time = round(float(data_row[7])*60)
            spectrum = (dict[target_file][mass][time])
            if type(spectrum)!= str:
                spectrum=(dict[target_file][mass][time+1])
            if type(spectrum) != str:
                spectrum = (dict[target_file][mass][time-1])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass-1][time])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass-1][time-1])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass-1][time+1])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass+1][time])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass+1][time-1])
            if type(spectrum)!=str:
                spectrum=(dict[target_file][mass+1][time+1])
            if type(spectrum)!=str:
                print("Error"+str(mass)+str(time))
            spectrum_list.append(spectrum)

            lgP=float(data_row[1])
            if lgP<=score_list[0]:
                prob_value=0
            elif lgP>=score_list[-1]:
                prob_value=1
            else:
                x1 = 0
                y1 = 0
                counter2 = 0
                for i in score_list:
                    if i >= lgP:
                        x2 = i
                        break
                    else:
                        counter2 += 1
                if counter2 >= 1:
                    x1 = score_list[counter2 - 1]
                    y1 = prob_list[counter2 - 1]
                    y2 = prob_list[counter2]

                prob_value = value(x1, x2, y1, y2, lgP)
            prob_value_list.append(prob_value)
   # add spectrum and PSM probability as columns to csv
    df=pd.read_csv(input)
    df['Spectrum'] = spectrum_list
    df['Prob_value']=prob_value_list
    df.to_csv(output,index=False)








