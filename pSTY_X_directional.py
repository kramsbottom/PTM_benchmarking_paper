import pandas as pd
import re
import numpy as np
import matplotlib.pyplot as plt
from files_direction import *

df=pd.read_csv(file)
#Calculate model FLR
df['PTM_final_prob'] = df['Score'] * df['PTM_score']
df = df.sort_values(by=(['PTM_final_prob']), ascending=[False])
df = df.reset_index(drop=True)
df['Count'] = (df.index) + 1
df['final_temp'] = 1 - df['PTM_final_prob']
df['PTM_final_prob_FLR'] = df['final_temp'].cumsum() / df['Count']
#5% FLR cutoff
df=df.loc[df['PTM_final_prob_FLR']<=0.05]

A_dist_all = []
G_dist_all = []
L_dist_all = []
D_dist_all = []
E_dist_all = []
P_dist_all = []
STY_dist_all = []
S_dist_all = []
peptides = df['Peptide']
pos = df['PTM_positions']
regex = re.compile(r'[S|T|Y]')
for p,q in zip(peptides,pos):
    A_dist = 10000
    G_dist = 10000
    L_dist = 10000
    D_dist = 10000
    E_dist = 10000
    P_dist = 10000
    STY_dist = 10000
    S_dist = 10000
    direction_A = 0
    direction_G = 0
    direction_L = 0
    direction_D = 0
    direction_E = 0
    direction_P = 0
    direction_S = 0
    direction_STY = 0
    if p[q-1]=="S" or p[q-1]=="T" or p[q-1]=="Y":
        peps=[]
        peps_pos=[]
        peps_neg=[]
        peps.append(p[q:])
        peps.append(p[:q - 1][::-1])
        peps_pos.append(p[q:])
        peps_neg.append(p[:q-1][::-1])
        for i in peps:
            if "A" in i:
                dist_temp = len(i.split('A')[0]) + 1
                if dist_temp < A_dist:
                    A_dist = dist_temp
                    if i in peps_pos:
                        direction_A=1
                    if i in peps_neg:
                        direction_A=-1
            if "G" in i:
                dist_temp = len(i.split('G')[0]) + 1
                if dist_temp < G_dist:
                    G_dist = dist_temp
                    if i in peps_pos:
                        direction_G=1
                    if i in peps_neg:
                        direction_G=-1
            if "L" in i:
                dist_temp = len(i.split('L')[0]) + 1
                if dist_temp < L_dist:
                    L_dist = dist_temp
                    if i in peps_pos:
                        direction_L=1
                    if i in peps_neg:
                        direction_L=-1
            if "D" in i:
                dist_temp = len(i.split("D")[0]) + 1
                if dist_temp < D_dist:
                    D_dist = dist_temp
                    if i in peps_pos:
                        direction_D=1
                    if i in peps_neg:
                        direction_D=-1
            if "E" in i:
                dist_temp = len(i.split("E")[0]) + 1
                if dist_temp < E_dist:
                    E_dist = dist_temp
                    if i in peps_pos:
                        direction_E=1
                    if i in peps_neg:
                        direction_E=-1
            if "P" in i:
                dist_temp = len(i.split("P")[0]) + 1
                if dist_temp < P_dist:
                    P_dist = dist_temp
                    if i in peps_pos:
                        direction_P=1
                    if i in peps_neg:
                        direction_P=-1
            if regex.search(i):
                dist_temp = len(re.split("[S|T|Y]", i)[0]) + 1
                if dist_temp < STY_dist:
                    STY_dist = dist_temp
                if i in peps_pos:
                    direction_STY=1
                if i in peps_neg:
                    direction_STY=-1

    if A_dist!=10000:
        A_dist_all.append(direction_A*A_dist)
    else:
        A_dist_all.append(np.NaN)
    if G_dist!=10000:
        G_dist_all.append(direction_G*G_dist)
    else:
        G_dist_all.append(np.NaN)
    if L_dist!=10000:
        L_dist_all.append(direction_L*L_dist)
    else:
        L_dist_all.append(np.NaN)
    if D_dist!=10000:
        D_dist_all.append(direction_D*D_dist)
    else:
        D_dist_all.append(np.NaN)
    if E_dist!=10000:
        E_dist_all.append(direction_E*E_dist)
    else:
        E_dist_all.append(np.NaN)
    if P_dist!=10000:
        P_dist_all.append(direction_P*P_dist)
    else:
        P_dist_all.append(np.NaN)
    if STY_dist!=10000:
        STY_dist_all.append(direction_STY*STY_dist)
    else:
        STY_dist_all.append(np.NaN)
    if S_dist!=10000:
        S_dist_all.append(direction_S*S_dist)
    else:
        S_dist_all.append(np.NaN)

df2=pd.DataFrame({"Peptides":peptides,"Positions":pos,"A_dist":A_dist_all,"G_dist":G_dist_all,"L_dist":L_dist_all,"D_dist":D_dist_all,"E_dist":E_dist_all,"P_dist":P_dist_all,"STY_dist":STY_dist_all,"S_dist":S_dist_all})
A_dist_list = df2['A_dist'].dropna().values
STY_dist_list = df2['STY_dist'].dropna().values
S_dist_list = df2['S_dist'].dropna().values
L_dist_list = df2['L_dist'].dropna().values
G_dist_list = df2['G_dist'].dropna().values
D_dist_list = df2['D_dist'].dropna().values
E_dist_list = df2['E_dist'].dropna().values
P_dist_list = df2['P_dist'].dropna().values

fig,axes=plt.subplots(nrows=3,ncols=2,sharey=True,figsize=(14,10))
bins=[-5,-4,-3,-2,-1,0,1,2,3,4,5,6]
plt.setp(axes, xticks=[-4.5,-3.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5], xticklabels=['≤-5','-4','-3','-2','-1','0','1','2','3','4','≥5'])
ax0, ax1, ax2, ax3, ax4, ax5 = axes.flatten()

ax0.hist([np.clip(A_dist_list,-5,5)],bins=bins,density=True,width=1,label="A")
ax0.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax0.set_title("A")
ax0.legend(loc='upper right')
ax1.hist([np.clip(L_dist_list,-5,5)],bins=bins,density=True,width=1,label="L")
ax1.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax1.set_title("L")
ax1.legend(loc='upper right')
ax2.hist([np.clip(G_dist_list,-5,5)],bins=bins,density=True,width=1,label="G")
ax2.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax2.set_title("G")
ax2.legend(loc='upper right')
ax3.hist([np.clip(D_dist_list,-5,5)],bins=bins,density=True,width=1,label="D")
ax3.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax3.set_title("D")
ax3.legend(loc='upper right')
ax4.hist([np.clip(E_dist_list,-5,5)],bins=bins,density=True,width=1,label="E")
ax4.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax4.set_title("E")
ax4.legend(loc='upper right')
ax5.hist([np.clip(P_dist_list,-5,5)],bins=bins,density=True,width=1,label="P")
ax5.hist([np.clip(STY_dist_list,-5,5)],bins=bins,density=True,width=1,label="STY",alpha=0.5)
ax5.set_title("P")
ax5.legend(loc='upper right')

fig.tight_layout()
fig.text(0.5, 0, 'Min dist pSTY-X', ha='center')
fig.text(0, 0.5, 'Normalised frequency', va='center', rotation='vertical')
plt.savefig(output + "/STY-pX_comparison_dist_directional_overlap_usingSTY_5%modelFLR.jpg", dpi=300)