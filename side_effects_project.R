# Description: This script serves the data collection and cleaning for the bachelorthesis .
# of Jerome Blessing. (jerome.blessining[@]bayer.com). 
# 
# It works with data sets from SIDER 4.1 side effect ressource and calls the PubChem-API
# to obtain InChi representation of the drug molecules.
# After an update, it know works with the SDF representation of the molecules.


library("tidyr")
library("dplyr")
library("ChemmineR")
library("ChemmineOB")
library("RColorBrewer")

setwd("C:/Users/Samsung/source/repos/R/Side_effects")
source("C:/Users/Samsung/source/repos/R/Side_effects/functions.R")



drug_names <- readr::read_tsv("Data/drug_names.tsv", col_names = c("stitch_id_flat","drug_name"))

meddra_freq <- readr::read_tsv("Data/meddra_freq.tsv", col_names =c('stitch_id_flat',
                                                                    'stitch_id_sterio',
                                                                    'umls_cui_from_label',
                                                                    'placebo',
                                                                    'frequency',
                                                                    'lower',
                                                                    'upper',
                                                                    'meddra_type',
                                                                    'umls_cui_from_meddra',
                                                                    'side_effect_name'))


tibble::add_column(drug_names, pubchem_id = "")

#Converts the STITCH ID to a pubchem ID
for (id in drug_names[1]) {
  drug_names$pubchem_id <- stringr::str_replace(id, "(CID[10]0*)","")
}


drug_names$InChi <- "NA"

# Run only for first time
if(FALSE)
{
  for (id in 1:nrow(drug_names) ) {
    
    # Prgoramm executes so slow, making this code obsolete.
    # This block is necesseary as PubChem might block IP-Adresses if #API-Calls > 5/sec.
    # As a precautionary measure, the script will wait 1.1 sec after every 4th call
    # if(id %% 5 == 0){
    #   Sys.sleep(1.1)
    # }
    
    drug_names$InChi[id] <- get_inchi(drug_names$pubchem_id[id])
    
  }
}

# Run only for first time.
# write.table(drug_names, "Data/drugs_with_inchi.tsv", sep = "\t")

drug_names_inchi <- readr::read_tsv("Data/drugs_with_inchi.tsv")

meddra_freq_cleaned <- meddra_freq %>%
  dplyr::select(stitch_id_flat,frequency,meddra_type,side_effect_name,lower)

meddra_freq_cleaned <- dplyr::left_join(meddra_freq_cleaned, drug_names_inchi, by="stitch_id_flat")
meddra_freq_cleaned <- meddra_freq_cleaned %>% select(-row_num)

meddra_freq_cleaned <- meddra_freq_cleaned %>% 
  dplyr::filter(meddra_type != "LLT") %>% 
  dplyr::filter(lower >= 0.1)

meddra_freq_cleaned <- meddra_freq_cleaned %>% 
  dplyr::filter(!frequency %in% c("rare", "very rare", "uncommon", "postmarketing","infrequent"))
  

##################################################

ae_df <- meddra_freq_cleaned %>% 
  select(pubchem_id, side_effect_name,lower)

ae_df_spread <- ae_df %>% 
  mutate(id = row_number()) %>% 
  spread(side_effect_name,lower,fill = 0)

ae_df_spread <- ae_df_spread %>% 
  select(-id)

ae_df_spread <- ae_df_spread %>% 
  group_by(pubchem_id) %>% 
  summarise_all(funs(mean))

# Replaces all values unequal 0 with 1 (expression evaluats to TRUE(1) or FALSE(0))
ae_df_spread[,-1] <- (ae_df_spread[,-1] != 0)*1 

readr::write_tsv(ae_df_spread, "Data/drug_se_binary.tsv")
save(ae_df_spread, file = "drug_se_binary.rda", compress =TRUE)
##################################################



readr::write_tsv(meddra_freq_cleaned, "Data/drug_inchi_se_freq.tsv")

cid_inchi <- drug_names_inchi %>% 
  select(pubchem_id, InChi)

cid_inchi <- left_join(cmp_node, cid_inchi, by=c("Source"="pubchem_id"))

cid_inchi$InChi <- stringr::str_replace_all(cid_inchi$InChi,"[\r\n]","")


drug_names_inchi$sdf <- "NA"
drug_names_inchi&smiles <- "NA"
for(i in drug_names_inchi$InChi)
{
  drug_names_inchi$smiles <- ChemmineOB::convertFormat("inchi","sdf",i,options = data.frame(names="gen2d"))
}

sdfset <- ChemmineR::read.SDFset("Data/compounds.sdf")
valid <- ChemmineR::validSDF(sdfset)
sdfset <- sdfset[valid]

cids_valid <- drug_names$pubchem_id[valid]
unique_cids <- unique(meddra_freq_cleaned$pubchem_id)
sdfset@ID <- cids_valid 

apset <- ChemmineR::sdf2ap(sdfset)
sdfset_valid <- sdfset[sdfset@ID %in% unique_cids,]
apset_valid <- sdf2ap(sdfset_valid)

#Save all files to compressed .rda format. Makes it quicker to load in the future
save(sdfset, file = "sdfset.rda", compress = TRUE)
save(sdfset_valid, file = "sdfset_valid.rda", compress = TRUE)
save(apset, file = "apset.rda", compress = TRUE)
save(apset_valid, file ="apset_valid.rda", compress = TRUE)
load("apset_valid.rda")

cmp_sim_df <- data.frame(matrix(ncol = 700, nrow = 700))
colnames(cmp_sim_df) <- unique_cids
rownames(cmp_sim_df) <- unique_cids

# Calculate similarity: duration of 45 min currently. do not execute
# for(i in 1:nrow(cmp_sim_df))
# {
#   for(j in 1:ncol(cmp_sim_df))
#   {
#     cmp_sim_df[i,j] <- cmp.similarity(apset_valid[rownames(cmp_sim_df[i,])], apset_valid[colnames(cmp_sim_df[j])])
#   }
# }

sim_matrix <- as.matrix(cmp_sim_df)

cmp_long <- tibble::rownames_to_column(cmp_sim_df)

cmp_long <- cmp_long %>% 
  tidyr::gather(key = "Destination", value = "similarity",-rowname)

cmp_long$Type <- "undirected"
cmp_long <- cmp_long %>% dplyr::rename(Source = rowname)

cmp_long[is.na(cmp_long)] <- 0

cmp_node <- dplyr::left_join(cmp_long,drug_names, by=c("Source" = "pubchem_id")) %>% 
  dplyr::select(Source,drug_name) %>% 
  dplyr::distinct()

readr::write_tsv(cmp_long, "Data/PID_Similarity.tsv")
readr::write_tsv(cmp_node, "Data/PID_Node.tsv")
readr::write_tsv(cid_inchi, "Data/cid_drugname_inchi.tsv")

###############################################################
# create and read similarity matrix
# Matrix has to be spread again from gather savefile, 
# after sim_matrix was not saved  and takes a long time
# to be computed.


cmp_matrix <- readr::read_tsv("Data/PID_Similarity.tsv")

cmp_matrix <- cmp_matrix %>% 
  select(-Type) %>% 
  spread(Target, Similarity)

readr::write_tsv(cmp_matrix, "Data/cmp_similarity_matrix.tsv")
###############################################################
