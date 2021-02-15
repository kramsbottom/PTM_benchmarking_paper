#function for calculation pAla
#take any input, create pAla column
import pandas as pd
import re

#calculate decoy FLR from decoy input 
def calulate_decoy_FLR(input,decoy):
    if decoy=="pAla":
        x="A"
    if decoy=="pLeu":
        x="L"
    if decoy=="pGly":
        x="G"
    if decoy=="pAsp":
        x="D"
    if decoy=="pGlu":
        x="E"
    if decoy=="pPro":
        x="P"
    df = pd.read_csv(input,dtype={'PTM_positions': str})
    S_count=df['Peptide'].str.count('S').sum()
    T_count=df['Peptide'].str.count('T').sum()
    Y_count=df['Peptide'].str.count('Y').sum()
    A_count=df['Peptide'].str.count(x).sum()
    STY_count=S_count+T_count+Y_count
    STY_A_ratio=STY_count/A_count

    pA_count = 0
    df.fillna('-')
    for i in range(len(df)):
        for a in (str(df.loc[i,'PTM_positions']).split(";")):
            if a!="-" and a!="nan":
                peptide=re.sub('[()+0-9. "]', '',df.loc[i,'Peptide'])
                if peptide[int(float(a))-1]==x:
                    pA_count+=1
        df.loc[i, 'p'+x+'_count'] = pA_count
        #decoy pX_FLR = STY:X ratio * pX_count * 2 / Count
        df.loc[i,decoy+'_FLR']=(STY_A_ratio*df.loc[i,'p'+x+'_count']*2)/(i+1)
    output=input.replace(".csv","_"+decoy+".csv")
    df.to_csv(output)