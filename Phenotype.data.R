setwd("/")

library(dplyr)
library(tibble)

# load muttable for list of tx421 pts
load('camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/SNV/Mut_table/20200129/patientMutTable.RData')

# create df of patients

phenotype.df <- data.frame(SampleID = rep(NA, length(unique(muttable_df$SampleID))))
phenotype.df$SampleID <- unique(muttable_df$SampleID)
pattern = "([A-Z]{3})([0-9]{3})"
phenotype.df$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), phenotype.df$SampleID)


# things to add 

# clinical data - physical activity, smoking, prev cancer, fhx cancer, etoh, ethnicity
load("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/working/clinicalData/tracerx_clin.nooutcomes.20200506.ERBB.RData")
tracerx_clin2 <- tibble::rownames_to_column(tracerx_clin2, "Hospital_ID_short")
clinical_phenotypes_of_interest <- tracerx_clin2[,c("Hospital_ID_short", "FreqPhyActivity", "DEMSmkPerDay","PrevCancer", "FamHistCaMember",
                                                    "DEMEthnicity", "DaysSinceBirth", "REGSex", "PathHist.R_TRACERx_Lesion1Form")]
# make blank cells NA
empty_as_na <- function(x){
  if("factor" %in% class(x)) x <- as.character(x) 
  ifelse(as.character(x)!="", x, NA)
}
clinical_phenotypes_of_interest <- clinical_phenotypes_of_interest %>% mutate_each(funs(empty_as_na)) 

phenotype.df <- merge(phenotype.df, clinical_phenotypes_of_interest, by = "Hospital_ID_short", all.y = FALSE)

# mutation signatures

signature_weights <- read.csv("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/SNV/Mut_Signatures/20200505/signature_weights_tumour.txt", sep = "\t", stringsAsFactors = F, header = T)
signature_counts <- read.csv("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/SNV/Mut_Signatures/20200505/signature_counts_tumour.txt", sep = "\t", stringsAsFactors = F, header = T)
colnames(signature_weights) <- paste0(colnames(signature_weights), "weights")
colnames(signature_counts) <- paste0(colnames(signature_counts), "counts")
# removing second cluster of two tumour cases for now
pattern = "([A-Z]{3})([0-9]{3})(_Cluster1)"
signature_counts$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), signature_counts$Tumourcounts)
signature_weights$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), signature_weights$Tumourweights)
phenotype.df <- merge(phenotype.df, signature_counts[,-1], by = "Hospital_ID_short", all.y = FALSE)
phenotype.df <- merge(phenotype.df, signature_weights[,-1], by = "Hospital_ID_short", all.y = FALSE)


# 3. oncogene associations (EGFR) - muttable - could also consider scna - ask Tom W
pattern = "([A-Z]{3})([0-9]{3})"
muttable_df$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), muttable_df$SampleID)
EGFR_pts <- subset(muttable_df, Hugo_Symbol == "EGFR" & driverCategory %in% c("1A", "2A"))
phenotype.df$EGFR_mut[phenotype.df$Hospital_ID_short %in% EGFR_pts$Hospital_ID_short] <- "Yes"
phenotype.df$EGFR_mut[!phenotype.df$Hospital_ID_short %in% EGFR_pts$Hospital_ID_short] <- "No"

# 4. tmb
load('camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/SNV/Mut_table/20200129/RegionMutTable.RData')
TMB.df.patient <- muttable_region_df %>%
  filter(var_count != 0) %>%
  group_by(RegionID) %>%
  summarise(TMB = n())
pattern = "([A-Z]{3})([0-9]{3})(.+)"
TMB.df.patient$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), TMB.df.patient$RegionID)
for (i in 1:nrow(TMB.df.patient)){
  TMB.df.patient$av.TMB[i] <- mean(TMB.df.patient$TMB[TMB.df.patient$Hospital_ID_short == TMB.df.patient$Hospital_ID_short[i]])
}
TMB.df.patient.unique <- TMB.df.patient[,c("Hospital_ID_short", "av.TMB")]
TMB.df.patient.unique <- TMB.df.patient.unique[!duplicated(TMB.df.patient.unique$av.TMB),]
phenotype.df <- merge(phenotype.df, TMB.df.patient.unique, by = "Hospital_ID_short", all = TRUE)

# 5. wgd frequency
wgd <- read.table("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/WGD/20191127/wgd_df.txt", sep = "\t", stringsAsFactors = F, header = T)
pattern = "([A-Z]{3})([0-9]{3})"
wgd$Hospital_ID_short <- gsub(pattern = pattern, paste0("\\1", "0", "\\2"), wgd$patient_id)
# removing second cluster of two tumour cases for now
wgd <- wgd[- grep("Cluster2", wgd$tumour_id),]
wgd_to_merge <- wgd[,c("gd_status", "Hospital_ID_short")]
phenotype.df <- merge(phenotype.df, wgd_to_merge, by = "Hospital_ID_short", all = TRUE)
  
# 6. wgii scores - don't know where these are


# 8. hla loh
# 17/05/20 JRMB - havent worked with this data before and not sure how best to use it, so leaving it out of df for now
LOHHLA <- read.csv("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/release/LOHHLA/421_lohhla.csv", stringsAsFactors = F, header = T)

# 9. danaher scores 
load("camp/project/proj-tracerx-lung/tctProjects/lungTx/Tx421/working/Expression/DanaherImmuneScores/til_score_tx400_patient_df.Rda")
til_score_tx400_patient_df$Hospital_ID_short <- til_score_tx400_patient_df$TRACERxID
phenotype.df <- merge(phenotype.df, til_score_tx400_patient_df, by = "Hospital_ID_short", all.x = TRUE)

save(phenotype.df, file = "camp/lab/swantonc/working/blackj/Germline_project/phenotype.df.Rdata")

