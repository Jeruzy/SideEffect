# Vizualize heatmaps and clusters of side-effect similarity and chemical similarity



library(tidyr)
library(ggplot2)
library(ChemmineR)


# Import datasets from 
cmp_sim_mat <- readr::read_tsv("Data/cmp_similarity_matrix.tsv")

cmp_ae_sim <- readr::read_tsv("Data/drug_se_binary.tsv")

cid_smiles <- readr::read_tsv("Data/SMILES_Pubchem.tsv", col_names = TRUE)

column_to_rownames(cmp_ae_sim, )

ae_smiles <- dplyr::left_join(cmp_ae_sim, cid_smiles)

ae_smiles <- ae_smiles %>% 
  select(SMILES, pubchem_id, everything())

#dist_ae <- dist.(cmp_ae_sim[,2:ncol(cmp_ae_sim)])

ChemmineR::cmp.cluster()
