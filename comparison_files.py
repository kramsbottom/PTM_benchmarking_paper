#Working directory
p="D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058"
#Working directory of each search to be compared
TPP_wd = p + "\TPP_tryptic" #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
Peaks_wd = p + "\PXD007058_SyntheticPeptides_PEAKS_43"
MQ_wd = p + "\MaxQuant-Ananth-Trypticsearch-nodecoys"
PD_wd = p + "\PD\PD_TD_STYA/NL_FALSE"
#Location of .mgf search files
search_files_location = "Z:\hlkramsb\PXD007058\Peaks\Input files"
#location of search database
database = "Z:\hlkramsb\Search_databases\pools_crap_targetdecoy.fasta"
#Synthetic peptide fasta with pools
s_pep_pools = "D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058\pool_ALL_name.fasta"
#Synthetic peptide fasta without pools
s_pep = "D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058\pool_ALL_name_no_pool.fasta"
#PXD ID
PXD = "PXD007058"
#FDR cutoff (0.0-1.0)
FDR_cutoff = 0.01
#List all working directories
w = [TPP_wd, Peaks_wd, MQ_wd, PD_wd]
#list all software names
s = ["TPP", "PEAKS", "MaxQuant", "PD"]

# #Working directory
# p="D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058"
# #Working directory of each search to be compared
# TPP_wd = p + "\TPP_tryptic" #if decoy amino acid not specified in folder name (ie pL/pLeu), will automatically use pAla method
# wd = p + "\TPP_pGly"
# wd2 = p + "\TPP_pLeu"
# wd3 = p + "\TPP_pD"
# wd4 = p + "\TPP_pE"
# wd5 = p + "\TPP_pP"
# #Location of .mgf search files
# search_files_location = "Z:\hlkramsb\PXD007058\Peaks\Input files"
# # #location of search database
# database = "Z:\hlkramsb\Search_databases\pools_crap_targetdecoy.fasta"
# # #Synthetic peptide fasta with pools
# s_pep_pools = "D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058\pool_ALL_name.fasta"
# # #Synthetic peptide fasta without pools
# s_pep = "D:\Dropbox\PTMExchange\PXD007058\PA_Phospho_PXD007058\pool_ALL_name_no_pool.fasta"
# # #PXD ID
# PXD = "PXD007058"#FDR cutoff (0.0-1.0)
# FDR_cutoff = 0.01
# #List all working directories
# w = [TPP_wd, wd, wd2, wd3, wd4, wd5]
# #list all software names
# s = ["TPP","TPP","TPP","TPP","TPP","TPP"]