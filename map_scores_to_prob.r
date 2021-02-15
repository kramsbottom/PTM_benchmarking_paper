library(ggplot2)
library(stringr)
library(tidyverse)

code_folder = "C:/Users/krams/Dropbox/PTMExchange/Scripts(26.11.19)/Modular/"
setwd("C:/Users/krams/Dropbox/PTMExchange/PXD008355/All_tryptic/PXD008355_all_PEAKS_44")
df = read.csv("DB search psm.csv")

decoy_df = subset(df,substr(df$Accession,0,5)=="DECOY")
target_df = subset(df,substr(df$Accession,0,5)!="DECOY")

bins = seq(0, 200, length.out = 50)
decoy_hist = hist(decoy_df$X.10lgP,  breaks = bins)
target_hist = hist(target_df$X.10lgP, breaks = bins)

final_probs = (target_hist$counts - decoy_hist$counts) / target_hist$counts

map_p_to_probs = data.frame(score_mid_points = target_hist$mids, probabilities = final_probs )
write.csv(map_p_to_probs,"map_p_to_probs.csv")

df_ascore=data.frame(AScore=df$AScore)
df_ascore$AScore <- as.character(df_ascore$AScore)
df_ascore=df_ascore %>%
  mutate(AScore=strsplit(AScore,";")) %>%
  unnest(AScore)
df_ascore= df_ascore %>%
  separate(AScore,c("Pos","PTM","Score"),":")
df_ascore$Residue=gsub('[[:digit:]]+', '', df_ascore$Pos)
df_ascore$Score<- as.numeric(df_ascore$Score)
df_ascore=filter(df_ascore, Score<=200)
df_sty = subset(df_ascore,df_ascore$Residue == "S" | df_ascore$Residue == "T" | df_ascore$Residue == "Y")
df_pAlas = subset(df_ascore,df_ascore$Residue == "A") 

#A_to_STY_ratio = 0.343

ptm_decoy_hist = hist(df_pAlas$Score,  breaks = bins)
ptm_target_hist = hist(df_sty $Score, breaks = bins)


magic_ratio = ptm_decoy_hist$counts[1] /ptm_target_hist$counts[1]
magic_ratio

final_probs_ptms = (ptm_target_hist$counts - ptm_decoy_hist$counts /magic_ratio ) / ptm_target_hist$counts
map_p_to_ptm_probs = data.frame(score_mid_points = ptm_target_hist$mids, probabilities = final_probs_ptms )
write.csv(map_p_to_ptm_probs,"map_p_to_ptm_probs.csv")
