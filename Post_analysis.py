import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import pAla

def explode(df, lst_cols, fill_value='', preserve_index=False):
    # make sure `lst_cols` is list-alike
    if (lst_cols is not None
        and len(lst_cols) > 0
        and not isinstance(lst_cols, (list, tuple, np.ndarray, pd.Series))):
        lst_cols = [lst_cols]
    # all columns except `lst_cols`
    idx_cols = df.columns.difference(lst_cols)
    # calculate lengths of lists
    lens = df[lst_cols[0]].str.len()
    # preserve original index values
    idx = np.repeat(df.index.values, lens)
    # create "exploded" DF
    res = (pd.DataFrame({
                col:np.repeat(df[col].values, lens)
                for col in idx_cols},
                index=idx)
             .assign(**{col:np.concatenate(df.loc[lens>0, col].values)
                            for col in lst_cols}))
    # append those rows that have empty lists
    if (lens == 0).any():
        # at least one list in cells is empty
        res = (res.append(df.loc[lens==0, idx_cols], sort=False)
                  .fillna(fill_value))
    # revert the original index order
    res = res.sort_index()
    # reset index if requested
    if not preserve_index:
        res = res.reset_index(drop=True)
    return res

#Explode rows for each protein PTM position - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on each protein
def site_based_all(input):
    df = pd.read_csv(input)
    df.dropna(subset=['PTM'], inplace=True)
    df['Protein'] = df['Protein'].str.split(":")
    df['PTM_Protein_Positions'] = df['PTM_Protein_Positions'].str.split(':')
    df = explode(df, ['Protein', 'PTM_Protein_Positions'], fill_value='')
    df['PTM'] = df['PTM'].str.split(';')
    df['PTM_score'] = df['PTM_score'].str.split(';')
    df['PTM_positions'] = df['PTM_positions'].str.split(';')
    df['PTM_Protein_Positions'] = df['PTM_Protein_Positions'].str.split(';')
    df = explode(df, ['PTM', 'PTM_score', 'PTM_positions', 'PTM_Protein_Positions'], fill_value='')

    df = df[df.PTM == "Phosphorylation"]
    df=df[~df.Protein.str.contains("DECOY")]
    df=df[~df.Protein.str.contains("CONTAM")]
    df = df.reset_index(drop=True)

    output=input.replace(".csv","_Site-based_all_proteins.csv")
    df.to_csv(output)

#Explode rows for each PSM postion - for each mod in peptide, a seperate row will show PSM score and corresponding PTM position/score for each site on the peptide
def site_based(input):
    df = pd.read_csv(input)
    df.dropna(subset=['PTM'], inplace=True)
    df=df.fillna('')
    df['PTM'] = df['PTM'].str.split(';')
    df['PTM_score'] = df['PTM_score'].str.split(';')
    df['PTM_positions'] = df['PTM_positions'].str.split(';')
    df = explode(df, ['PTM', 'PTM_score', 'PTM_positions'], fill_value='')

    df = df[df.PTM == "Phosphorylation"]
    df = df[~df.Protein.str.contains("DECOY")]
    df = df[~df.Protein.str.contains("CONTAM")]
    df = df.reset_index(drop=True)

    output = input.replace(".csv", "_Site-based.csv")
    df.to_csv(output)

#Add column for if spectrum matches throughout comparisons
def spectrum_comparison(input,spectrum_match_list):
    df=pd.read_csv(input).fillna('N/A')
    match_counts=[]
    for i in range(len(df)):
        match_count=0
        spectrum=df.loc[i,'Spectrum']
        spectrum=spectrum.replace(".0",".")
        if spectrum in spectrum_match_list:
            match_count=1
        match_counts.append(match_count)
    df['Spectrum_matches']=match_counts
    output=input.replace(".csv","_spectrum_match.csv")
    df.to_csv(output,index=False)

#model FLR from combined probability: 1-(PSM prob*PTM prob)/Count
def model_FLR(file):
    df=pd.read_csv(file)
    df = df[df['PTM_positions'].notna()]
    df = df[df.PTM == "Phosphorylation"]
    df = df[~df.Protein.str.contains("DECOY",na=False)]
    df = df[~df.Protein.str.contains("CONTAM",na=False)]
    df = df.reset_index(drop=True)
    df['PTM_final_prob'] = df['Score'] * df['PTM_score']
    df = df.sort_values(by=(['PTM_final_prob']), ascending=[False])
    df = df.reset_index(drop=True)
    df['Count'] = (df.index) + 1
    df['final_temp'] = 1 - df['PTM_final_prob']
    df['PTM_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count']
    df = df.drop(['final_temp'], axis=1)
    df = df.reset_index(drop=True)
    output=file.replace(".csv","_FLR.csv")
    df.to_csv(output,index=False)

#Plot FLR comparisons - pX_FLR and model FLR for each software
def plot_FLR_comparisons(filter,file_list,software_list,working):
    folder=""
    legend=[]
    for i,j in zip(software_list,file_list):
        folder+="_"+i
        if "Site_confident" not in j:
            legend.append(i + " Model_FLR")
        if "pL" in j or "pLeu" in j:
            legend.append(i + " pLeu_FLR")
        elif "pG" in j or "pGly" in j:
            legend.append(i + " pGly_FLR")
        elif "pD" in j or "pAsp" in j:
            legend.append(i + " pAsp_FLR")
        elif "pE" in j or "pGlu" in j:
            legend.append(i + " pGlu_FLR")
        elif "pP" in j or "pPro" in j:
            legend.append(i + " pPro_FLR")
        else:
            legend.append(i + " pAla_FLR")
    output=working+"/Comparisons/"+folder[1:]
    if not os.path.exists(output):
        os.mkdir(output)
    #Combined prob FLR plot
    x_list=['PTM_final_prob','Count']
    for x in x_list:
        lines = ["blue","green","red","cyan","magenta","orange"]
        counter = 0
        max = 0
        for a, b in zip(file_list,lines):
            if "pGly" in a:
                p="pGly_FLR"
            elif "pLeu" in a:
                p="pLeu_FLR"
            elif "pAsp" in a:
                p="pAsp_FLR"
            elif "pGlu" in a:
                p="pGlu_FLR"
            elif "pPro" in a:
                p="pPro_FLR"
            else:
                p="pAla_FLR"
            df=pd.read_csv(a)
            if df[x].max()>max:
                max=df[x].max()
            if x == "PTM_final_prob":
                lim = (1, 0)
            else:
                lim=(0,max)
            if "Site_confident" not in file_list[0]:
                if counter==0:
                    ax=df.plot.line(x=x, y="PTM_final_prob_FLR", linestyle=":", color=b, xlim=lim, figsize=(14,7))
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax, color=b, xlim=lim)
                else:
                    df.plot.line(x=x, y="PTM_final_prob_FLR", linestyle=":",color=b, ax=ax, xlim=lim)
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax,color=b, xlim=lim)
                counter += 1
            else:
                if counter==0:
                    ax=df.plot.line(x=x, y=p, linestyle="-", color=b, xlim=lim, figsize=(14,7))
                else:
                    df.plot.line(x=x, y=p, linestyle="-", ax=ax,color=b, xlim=lim)
                counter += 1
        ax.legend(legend)
        ax.set_ylabel("FLR")
        if x=="PTM_final_prob":
            ax.set_xlabel("Combined probability")
        if "Site_confident_PTM_unique" in file_list[0]:
            plt.savefig(output+"/"+filter+"_"+x+"_All_software_Site_Confident_"+p+".jpg",dpi=300)
        else:
            plt.savefig(output + "/" + filter + "_" + x + "_All_software_"+p+".jpg", dpi=300)