library(ukbtools)

new_col <- function(ukb, pheno_pos, all_pheno_eids)
{
  pheno_pos = data.frame(V1=pheno_pos, V2=1)
  ukb = data.frame(eid=ukb[,"eid"])
  merge_tmp = merge(x = ukb, y = pheno_pos, by.x = "eid", by.y = "V1", all.x = T)
  merge_tmp[is.na(merge_tmp$V2),"V2"] = "0"
  merge_tmp[merge_tmp$V2 == "0" & merge_tmp$eid %in% all_pheno_eids, "V2"] = "-9"
  
  #col_ret = merge_tmp[,"V2"]
  return(merge_tmp)
}

setwd("~/jl11/Documents/PhD/UK biobank/")
# once
#ukb_data = ukb_df("ukb10043")
#saveRDS(ukb_data, file = "ukb_data.Rdata")

ukb_data = readRDS("ukb_data.Rdata")
all_pheno_eids = read.table("~/Documents/Postdoc/uk_biobank/all_pos_pheno.txt", header = F, stringsAsFactors = F)
all_pheno_eids = all_pheno_eids$V1

# test of ID matches to genetic data
merge_test = merge(genotyped, ukb_data, by.x = "X1", by.y = "eid")[,c("X1","X5","sex_0_0")]
merge_test # looks ok

# self reported sepsis eids
sr_sepsis_eids <- c(NA)
for(i in 15:101)
{
  sr_sepsis_eids = c(sr_sepsis_eids,ukb_data[!is.na(ukb_data[,i]) & ukb_data[,i] == "1657",1])
}
sr_sepsis_eids = sr_sepsis_eids[-1]

# self reported meningitis eids
sr_men_eids <- c(NA)
for(i in 15:101)
{
  sr_men_eids = c(sr_men_eids,ukb_data[!is.na(ukb_data[,i]) & ukb_data[,i] == "1247",1])
}
sr_men_eids = sr_men_eids[-1]

all_sr = unique(c(sr_men_eids, sr_sepsis_eids))

# Phenotypes
icd_codes <- read_delim("~/Documents/Postdoc/uk_biobank/icd_codes.txt", "\t", 
                        escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
colnames(icd_codes) = c("eid", "icd10")

all_icd = icd_codes[icd_codes$icd10 != "A872" & icd_codes$icd10 != "A879" & icd_codes$icd10 != "B021",]$eid
all_men = icd_codes[icd_codes$icd10 == "G01" | icd_codes$icd10 == "G001" | icd_codes$icd10 == "G002" | icd_codes$icd10 == "G003" | icd_codes$icd10 == "G008" | icd_codes$icd10 == "G009",]$eid
# too few
#pneu_men = icd_codes[icd_codes$icd10 == "G001",]$eid
pneu_sep = icd_codes[icd_codes$icd10 == "A403" | icd_codes$icd10 == "A409" | icd_codes$icd10 == "A408" | icd_codes$icd10 == "A40",]$eid

pheno = new_col(ukb_data, all_sr, all_pheno_eids)
colnames(pheno)[2] = "ALL_SR"
pheno = merge(pheno, new_col(ukb_data, sr_men_eids, all_pheno_eids))
colnames(pheno)[3] = "MEN_SR"
pheno = merge(pheno, new_col(ukb_data, sr_sepsis_eids, all_pheno_eids))
colnames(pheno)[4] = "SEP_SR"
pheno = merge(pheno, new_col(ukb_data, all_icd, all_pheno_eids))
colnames(pheno)[5] = "ALL_ICD"
pheno = merge(pheno, new_col(ukb_data, all_men, all_pheno_eids))
colnames(pheno)[6] = "MEN_ICD"
pheno = merge(pheno, new_col(ukb_data, pneu_sep, all_pheno_eids))
colnames(pheno)[7] = "SEP_ICD"
pheno = data.frame(FID=pheno$eid, IID=pheno$eid, pheno[,-1])

write.table(pheno, file = "phenotypes.txt", quote = F, col.names = T, row.names = F, sep = " ")

# Covariates
dobs = as.Date(paste(ukb_data[,"year_of_birth_0_0"], ukb_data[,"month_of_birth_0_0"], "1", sep = "-"), format = "%Y-%B-%d")
ages = (ukb_data[,"date_of_attending_assessment_centre_0_0"] - dobs)/365
covariates = data.frame(FID=ukb_data[,"eid"], IID=ukb_data[,"eid"], AGE=ages, CENTRE=ukb_data[,"uk_biobank_assessment_centre_0_0"])
write.table(covariates, file="ukbb.covar", quote = F, row.names = F)
