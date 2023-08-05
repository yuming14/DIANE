library(readr)
library(dplyr)
library(stringr)

MS = read_csv("/n/data1/hsph/biostat/celehs/lab/SHARE/UPMC/MS/data_table/data_processed/UPMC_MS_2011_to_2021_Codified_monthly_counts_2023-05-23.csv")
MS = MS[grep(":",MS$code),]
AD = read_csv("/n/data1/hsph/biostat/celehs/lab/SHARE/UPMC/AD/data_table/data_processed/UPMC_AD_2011_to_2021_Codified_NLP_monthly_counts_2023-05-15.csv")
AD = AD[grep(":",AD$code),]
dict = data.frame(code = unique(c(unique(AD$code),
                                  unique(MS$code))))
dict$code = gsub("CCS-PCS","CCS",dict$code)
dict$group = sapply(dict$code, function(x) strsplit(x,":")[[1]][1])
table(dict$group)
# CCS      CPT4     HCPCS   ICD10CM  ICD10PCS    ICD9CM LOCAL|LAB  LOCAL|PX     LOINC 
# 247       195       138      1752        19       821      2775      3329      4267 
# NDC   PheCode    RXNORM 
# 4149      1851      2496 
sum(dict$group%in%c("CCS","LOINC","PheCode","RXNORM"))
# 8861
save(dict, file = "~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")

intersect(unique(AD$patient_num),unique(MS$patient_num))
# 9332036023 9332192647 9332806774 9332827350 9332846860 9332861771 9332868047
month_count = rbind(AD[,c('patient_num','num_months','code','count')],
                    MS[,c('patient_num','num_months','code','count')])
rm(MS, AD)
save(month_count, file = "~/Ising_model/result/2306/UPMC_MS_AD_2011_to_2021_Codified_monthly_occur.Rdata")
patient_list = unique(month_count$patient_num)

n = length(patient_list)
df = month_count %>%
  filter(patient_num %in% patient_list)

pe_code = lapply(patient_list, function(pe){
  return(unique(df$code[which(df$patient_num==pe)]))
})
save(pe_code, file = "~/Ising_model/result/2306/UPMC_MS_AD_Codified_per_patient.Rdata")
pe_code = unlist(pe_code)
pe_code = sort(table(pe_code), decreasing = TRUE)
sum(pe_code>n*0.01)
# 3275
sum(pe_code>n*0.001)
# 8351

code = names(pe_code)[which(pe_code>n*0.01)]
idx = grep("PheCode",code)
phecode = code[idx]
code = code[-idx]
phecode = sapply(phecode, function(x) strsplit(x,":")[[1]][2])
phecode_int = floor(as.numeric(phecode))
phecode_select = lapply(unique(phecode_int), function(x){
  idx = which(phecode_int==x)
  phe = phecode[idx]
  nphe = nchar(phe)
  x = str_pad(x, 3, "left", "0")
  if(length(phe)==1 | length(unique(nphe))==1){
    return(phe)
  }else{
    id_select = which(nphe==nchar(x)+2)
    if(length(id_select)>0){
      return(phe[id_select])
    }else{
      return(character())
    }
  }
})
summary(unlist(lapply(phecode_select, length)))
phecode_select = unlist(phecode_select)
code = c(code, names(phecode_select))
tail(code)
code = gsub('CCS-PCS',"CCS",code)
dict$select = dict$code%in%code
table(dict[,c('select','group')])
#       group
#select   CCS CPT4 HCPCS ICD10CM ICD10PCS ICD9CM LOCAL|LAB LOCAL|PX LOINC  NDC PheCode RXNORM
#  FALSE   94  174   133    1632       19    747      2603     3127  3495 3983    1146   1957
#. TRUE   153   21     5     120        0     74       172      202   772  166     705    539
save(dict, file = "~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")

##type 1 code to obtain samples ####
load("~/Ising_model/result/2306/UPMC_MS_AD_2011_to_2021_Codified_monthly_occur.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
code = dict$code[which(dict$select)]
code = gsub("CCS","CCS-PCS",code)
get_code_per_patient = function(df, windowsize = 6){
  mon = seq(min(df$num_months), max(df$num_months), windowsize)
  codelist = list()
  if(length(mon)==1){
    codelist[[1]] = unique(df$code)
  }else{
    for(i in 2:length(mon)){
      idx = which(df$num_months<mon[i])
      codelist[[i-1]] = unique(df$code[idx])
      df = df[-idx,]
    }
  }
  return(codelist[which(lapply(codelist,length)>0)])
}

get_code_vector = function(codelist, code){
  code_vector = lapply(codelist, function(x){
    tmp = code %in% x
    return(tmp*2-1)
  })
  code_vector = do.call("rbind", code_vector)
  colnames(code_vector) = code
  return(code_vector)
}

get_codelist = function(df, patient, code, windowsize = 6, file = NULL){
  code_matrix = matrix(nrow = 0, ncol = length(code))
  for(i in 1:length(patient)){
    pe = patient[i]
    idx = which(df$patient_num==pe)
    tmp = get_code_per_patient(df[idx,], windowsize)
    code_matrix = rbind(code_matrix, get_code_vector(tmp, code))
    df = df[-idx,]
    if(i%%200==0 & !is.null(file)){
      save(df, i, patient, code_matrix,
           file = file)
    }
  }
  colnames(code_matrix) = code
  return(code_matrix)
}

length(unique(month_count$patient_num))
# 17428

sample_matrix = get_codelist(df = month_count, 
                             patient = unique(month_count$patient_num),
                             code = code, windowsize = 6,
                             file = "~/Ising_model/result/2306/MS_AD_sample_matrix.Rdata")
colnames(sample_matrix) = gsub("CCS-PCS","CCS",colnames(sample_matrix))
save(sample_matrix, file = "~/Ising_model/result/2306/MS_AD_sample_matrix.Rdata")
rm(list = ls())

##type 2 code to obtain samples ####
load("~/Ising_model/result/2306/UPMC_MS_AD_2011_to_2021_Codified_monthly_occur.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
code.list = dict$code[which(dict$select)]
code.list = gsub("CCS","CCS-PCS",code.list)

month_count$month.idx = ceiling((month_count$num_months +2) /6)
month_count = month_count %>%
  filter(code %in% gsub("CCS","CCS-PCS",dict$code[which(dict$select)]))
month_count = month_count %>% group_by(patient_num, month.idx, code) %>% 
  summarise(total.count = sum(count))
month_count = tidyr::spread(month_count, key = code, value = total.count, fill = 0)
dim(month_count)
# 236875   2928
save(month_count, file = "~/Ising_model/result/2306/MS_AD_sample_count_matrix.Rdata")

sample_matrix = month_count[,-c(1:2)]
sample_matrix = as.matrix(sample_matrix)
sample_matrix[sample_matrix>0] = 1
sample_matrix[sample_matrix==0] = -1
rs = rowSums(sample_matrix)
min(rs)
# -2926
colnames(sample_matrix) = gsub("CCS-PCS","CCS",colnames(sample_matrix))
save(sample_matrix, file = "~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")

## scaled AUC ####

load("~/Ising_model/result/2306/MS_AD_codified_pairs.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
source("~/Ising_model/code/2305/get_AUC_MS.R")

load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_15.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5946533 20075
# related RxNorm from UMLS 0.8216099  1107
# RxNorm Hierachy (rm)     0.7525460   703
# RxNorm Hierachy          0.8118664   660
# PheCode Hierachy         0.9353133   510
# LOINC Hierachy (LP1)     0.6207692   477
AUC1 = rbind(AUC[[1]],AUC[[2]])
AUC1 = AUC1[match(c("RxNorm Hierachy","PheCode Hierachy","classifies","may_prevent","ssc","ddx","may_treat"), rownames(AUC1)),]
AUC1$group = rownames(AUC1)
AUC1$type = factor(c(rep("similar pairs",2), rep("related pairs",5)),
                   levels = c("similar pairs","related pairs"))
AUC1$group = factor(AUC1$group, levels = c("PheCode Hierachy","RxNorm Hierachy",
                                           "classifies","may_prevent","ssc",
                                           "ddx","may_treat"),
                    labels = c("PheCode Hierachy","RxNorm Hierachy",
                               "classifies","may prevent","ssc",
                               "ddx","may treat"))
colnames(AUC1) = c("AUC","num","group","type")
library(ggplot2)
ggplot(aes(x = group, y = AUC, fill = type), data = AUC1) +
  geom_bar(stat = "identity") + 
  labs(x = NULL) + 
  coord_cartesian(ylim=c(0.6,1)) + 
  theme_bw() + 
  theme(text = element_text(size = 22)) +
  theme(legend.position="top") +
  theme(axis.text.x = element_text(angle = 30, hjust = 0.5, vjust = 0.6))
