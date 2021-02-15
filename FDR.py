import Extract_DF
import numpy as np
import matplotlib.pylab as plt
import time

start_time = time.time()

#calculate FDR
def calculateFDR(results_file, software, output,PXD,working,search):
    #extract results to df
    if software == "PEAKS":
        df=Extract_DF.extract_PSM_df(working,results_file,PXD)
    if software == "TPP":
        df=Extract_DF.extract_PTMprophet_df(results_file,PXD)
    if software =="MaxQuant":
        df=Extract_DF.extract_maxQuant_df(results_file,PXD,search)
    if software=="PD":
        df=Extract_DF.extract_PD_mascot(results_file,PXD)
    df['Score']=df['Score'].astype(float)
    df=df.sort_values(by='Score',ascending=False)
    df=df.reset_index(drop=True)
    #set decoy and target
    df['Protein_count'] = df['Protein'].str.count(":")+1
    df['Decoy_count'] = df['Protein'].str.count("DECOY")
    df['Decoy'] = np.where(df['Protein_count']==df['Decoy_count'],1,0)
    df=df.drop(columns=['Protein_count', 'Decoy_count'])
    #target_count = row_count-decoy_count
    df['decoy_count']=df['Decoy'].cumsum()
    df['row_count']=(df.index)+1
    df['target_count']=df['row_count']-df['decoy_count']
    #FDR = decoy_count/target_count
    df['FDR']= df['decoy_count']/df['target_count']
    df['FDR'] = df['FDR'].astype(float)
    #q_val = min FDR at this position or lower
    df['q_value']=df['FDR']
    df['q_value']=df.iloc[::-1]['FDR'].cummin()
    #return as csv
    df.to_csv(output,index=False)

    #print("Create plots --- %s seconds ---" % (time.time() - start_time))
    #FDR plots
    ax1=df.plot.scatter(x='FDR',y='row_count')
    ax2=df.plot.line(x='q_value',y='row_count')
    plt.ylabel("PSM count")
    plt.savefig('FDR.jpg',dpi=300)

    ax1=df.plot.scatter(x='FDR',y='row_count')
    ax2=df.plot.line(x='q_value',y='row_count',xlim=(0,0.2),ylim=(0,10000))
    plt.ylabel("PSM count")
    plt.savefig('FDR_filter.jpg',dpi=300)

    df.plot.line(x='Score', y='FDR',ylim=(0,1),xlim=(1,0))
    plt.xlabel("Peptide Probabilty")
    plt.ylabel("Global FDR")
    plt.savefig('FDR_score.jpg',dpi=300)
    plt.close('all')
    #print("--- %s seconds ---" % (time.time() - start_time))
