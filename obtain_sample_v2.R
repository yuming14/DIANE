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

## obtain Theta ####
load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
source("~/Ising_model/code/2306/Ising_function_0530.R")
est_Theta = est_Theta_DC_nonconvex(x = sample_matrix, 
                                   Theta0 = NULL, m = 8, eta = 0.005, d = 200, 
                                   maxstep = c(20,20), epsilon = c(1e-2,1e-2),
                                   diag_add = 0.01, PSD = TRUE,
                                   Theta_star = NULL, info = FALSE,
                                   lambda = 0.001, maxstep_convex = 20,
                                   ini_convex = 1, timing = TRUE,
                                   stopsign = "Theta", 
                                   file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01.Rdata")
save(est_Theta, file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01.Rdata")

load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
source("~/Ising_model/code/2306/Ising_function_0530.R")
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01.Rdata")
est_Theta = est_Theta_DC_nonconvex(x = sample_matrix, 
                                   Theta0 = Theta_new, m = 8, eta = 0.005, d = 200, 
                                   maxstep = c(20,20), epsilon = c(1e-2,1e-2),
                                   diag_add = 0.01, PSD = TRUE,
                                   Theta_star = NULL, info = FALSE,
                                   lambda = 0.001, maxstep_convex = 20,
                                   ini_convex = 1, timing = TRUE,
                                   stopsign = "Theta", 
                                   file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01_1.Rdata")
save(est_Theta, file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01_1.Rdata")


load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_2.Rdata")
eta = 0.001
est_Theta = est_Theta_DC_nonconvex(x = sample_matrix, 
                                   Theta0 = Theta_new, m = 8, eta = 0.001, d = 250, 
                                   maxstep = c(50,75), epsilon = c(1e-10,1e-10),
                                   diag_add = 0, PSD = TRUE,
                                   Theta_star = NULL, info = FALSE,
                                   lambda = 0.001, maxstep_convex = 100,
                                   ini_convex = 1, timing = TRUE,
                                   stopsign = "Theta", 
                                   file = paste("~/Ising_model/result/2306/MS_AD_Theta_est_eta_",eta,"_0620_3.Rdata",sep=""))
save(est_Theta, file = paste("~/Ising_model/result/2306/MS_AD_Theta_est_eta_",eta,"_0620_3.Rdata",sep=""))

## Weijing's code ####

load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
source("~/Ising_model/code/2305/ref/function_Ising.R")
Theta0 = LR_KLAIM(X = sample_matrix, M=NULL, B=NULL, gamma=NULL, 
                  lambda_LR=0.01, 
                  lambda_transfer=0, tau_B=0.01, tau_gamma=0, 
                  max.iter = 100, tol = 1e-4,
                  file = "~/Ising_model/result/2306/MS_AD_Theta_est_tau_0.01_lambda_0.01.Rdata")

## nonconvex after weijing's code ####
load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_0.01_lambda_0.01.Rdata")
load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
source("~/Ising_model/code/2306/Ising_function_0530.R")
rownames(B) = colnames(B) = colnames(sample_matrix)
rm(gamma, i)
est_Theta = est_Theta_DC_nonconvex(x = sample_matrix, 
                                   Theta0 = B, m = 8, eta = 0.001, d = 200, 
                                   maxstep = c(20,20), epsilon = c(1e-2,1e-2),
                                   diag_add = 0.01, PSD = TRUE,
                                   Theta_star = NULL, info = FALSE,
                                   lambda = 0.001, maxstep_convex = 50,
                                   ini_convex = 1, timing = TRUE,
                                   stopsign = "Theta", 
                                   file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_ini_weijing.Rdata")
save(est_Theta, file = "~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_ini_weijing.Rdata")



## AUC ####
library(dplyr)
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.005_diag_0.01.Rdata")
load("/n/data1/hsph/biostat/celehs/lab/zig728/Ising_Model/PairsWithNull_MS.Rdata")
pairs0 = pairs
load("~/Ising_model/result/2306/pairs_AD_codified.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
pairs1 = pairs
pairs = rbind(pairs0, pairs1) %>%
  select(code1, code2, source, type, group) %>%
  unique()
source("~/Ising_model/code/2305/get_AUC_MS.R")
get_AUC_embed(U, pairs, WithNull = FALSE,
              dict = dict, normalize = FALSE)

load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_0.01_lambda_0.01.Rdata")
load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_pairs.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
rownames(B) = colnames(B) = colnames(sample_matrix)
source("~/Ising_model/code/2305/get_AUC_MS.R")
fit = svd(B)
idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
rownames(U) = rownames(B)
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.5118201 20075
# related RxNorm from UMLS 0.5597491  1107
# RxNorm Hierachy (rm)     0.6368925   703
# RxNorm Hierachy          0.6402984   660
# PheCode Hierachy         0.4757632   510
# LOINC Hierachy (LP1)     0.4473188   477
AUC = lapply(AUC, function(x){
  x = as.data.frame(x)
  x = x[order(x$auc, decreasing = TRUE),]
  return(x)
})

## scaled AUC ####

load("~/Ising_model/result/2306/MS_AD_codified_pairs.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
source("~/Ising_model/code/2305/get_AUC_MS.R")

load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0618.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.5219411 20075
# related RxNorm from UMLS 0.6760632  1107
# RxNorm Hierachy (rm)     0.6397779   703
# RxNorm Hierachy          0.6886777   660
# PheCode Hierachy         0.6342407   510
# LOINC Hierachy (LP1)     0.6046614   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0618_1.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5131946 20075
# related RxNorm from UMLS 0.6463492  1107
# RxNorm Hierachy (rm)     0.6494378   703
# RxNorm Hierachy          0.6818526   660
# PheCode Hierachy         0.5992080   510
# LOINC Hierachy (LP1)     0.5661212   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_2.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5564326 20075
# related RxNorm from UMLS 0.7488676  1107
# RxNorm Hierachy (rm)     0.6758881   703
# RxNorm Hierachy          0.7459320   660
# PheCode Hierachy         0.7828374   510
# LOINC Hierachy (LP1)     0.6104321   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_3.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5604831 20075
# related RxNorm from UMLS 0.7558389  1107
# RxNorm Hierachy (rm)     0.6842328   703
# RxNorm Hierachy          0.7537672   660
# PheCode Hierachy         0.8028643   510
# LOINC Hierachy (LP1)     0.6184882   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_4.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5657589 20075
# related RxNorm from UMLS 0.7636915  1107
# RxNorm Hierachy (rm)     0.6949287   703
# RxNorm Hierachy          0.7641483   660
# PheCode Hierachy         0.8293349   510
# LOINC Hierachy (LP1)     0.6273090   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_5.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5698206 20075
# related RxNorm from UMLS 0.7712504  1107
# RxNorm Hierachy (rm)     0.7041960   703
# RxNorm Hierachy          0.7720340   660
# PheCode Hierachy         0.8495348   510
# LOINC Hierachy (LP1)     0.6295945   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_6.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5732370 20075
# related RxNorm from UMLS 0.7784118  1107
# RxNorm Hierachy (rm)     0.7114986   703
# RxNorm Hierachy          0.7782415   660
# PheCode Hierachy         0.8659208   510
# LOINC Hierachy (LP1)     0.6312206   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_7.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5753222 20075
# related RxNorm from UMLS 0.7825091  1107
# RxNorm Hierachy (rm)     0.7153633   703
# RxNorm Hierachy          0.7817929   660
# PheCode Hierachy         0.8755440   510
# LOINC Hierachy (LP1)     0.6317305   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_8.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5783124 20075
# related RxNorm from UMLS 0.7883796  1107
# RxNorm Hierachy (rm)     0.7206526   703
# RxNorm Hierachy          0.7857117   660
# PheCode Hierachy         0.8881238   510
# LOINC Hierachy (LP1)     0.6323106   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_9.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5806180 20075
# related RxNorm from UMLS 0.7928669  1107
# RxNorm Hierachy (rm)     0.7246934   703
# RxNorm Hierachy          0.7887397   660
# PheCode Hierachy         0.8969127   510
# LOINC Hierachy (LP1)     0.6321612   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_10.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5826052 20075
# related RxNorm from UMLS 0.7969732  1107
# RxNorm Hierachy (rm)     0.7281211   703
# RxNorm Hierachy          0.7917493   660
# PheCode Hierachy         0.9038754   510
# LOINC Hierachy (LP1)     0.6324776   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_11.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5872054 20075
# related RxNorm from UMLS 0.8066488  1107
# RxNorm Hierachy (rm)     0.7367571   703
# RxNorm Hierachy          0.7981152   660
# PheCode Hierachy         0.9180123   510
# LOINC Hierachy (LP1)     0.6306757   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_12.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5902225 20075
# related RxNorm from UMLS 0.8127160  1107
# RxNorm Hierachy (rm)     0.7414333   703
# RxNorm Hierachy          0.8021625   660
# PheCode Hierachy         0.9255402   510
# LOINC Hierachy (LP1)     0.6280298   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_13.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5926357 20075
# related RxNorm from UMLS 0.8169863  1107
# RxNorm Hierachy (rm)     0.7461499   703
# RxNorm Hierachy          0.8055808   660
# PheCode Hierachy         0.9301192   510
# LOINC Hierachy (LP1)     0.6256521   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_14.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
# auc   num
# RxNorm High Level        0.5940344 20075
# related RxNorm from UMLS 0.8195861  1107
# RxNorm Hierachy (rm)     0.7503182   703
# RxNorm Hierachy          0.8095455   660
# PheCode Hierachy         0.9336409   510
# LOINC Hierachy (LP1)     0.6236040   477
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
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.0001.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.5220296 20075
# related RxNorm from UMLS 0.6763668  1107
# RxNorm Hierachy (rm)     0.6398831   703
# RxNorm Hierachy          0.6889486   660
# PheCode Hierachy         0.6343676   510
# LOINC Hierachy (LP1)     0.6047581   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_1e-04_0618.Rdata")
AUC = get_AUC_embed(est_Theta$U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.4838735 20075
# related RxNorm from UMLS 0.4058561  1107
# RxNorm Hierachy (rm)     0.5967860   703
# RxNorm Hierachy          0.5953076   660
# PheCode Hierachy         0.4601230   510
# LOINC Hierachy (LP1)     0.4125276   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_1e-04_0618_1.Rdata")
AUC = get_AUC_embed(est_Theta$U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.4800792 20075
# related RxNorm from UMLS 0.4661524  1107
# RxNorm Hierachy (rm)     0.4738785   703
# RxNorm Hierachy          0.4541644   660
# PheCode Hierachy         0.4583506   510
# LOINC Hierachy (LP1)     0.4752098   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_1e-05_0618.Rdata")
AUC = get_AUC_embed(est_Theta$U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.4811954 20075
# related RxNorm from UMLS 0.3808677  1107
# RxNorm Hierachy (rm)     0.4761427   703
# RxNorm Hierachy          0.4636226   660
# PheCode Hierachy         0.4561053   510
# LOINC Hierachy (LP1)     0.3757719   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_1e-05_0618_1.Rdata")
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
#                                auc   num
# RxNorm High Level        0.4819114 20075
# related RxNorm from UMLS 0.3787763  1107
# RxNorm Hierachy (rm)     0.4701027   703
# RxNorm Hierachy          0.4509183   660
# PheCode Hierachy         0.4554018   510
# LOINC Hierachy (LP1)     0.3602793   477

namelist = rownames(U)

load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_0.01_lambda_0_0618.Rdata")
fit = svd(B)
idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
rownames(U) = namelist
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
#                                auc   num
# RxNorm High Level        0.5314393 20075
# related RxNorm from UMLS 0.7033291  1107
# RxNorm Hierachy (rm)     0.6724948   703
# RxNorm Hierachy          0.7179293   660
# PheCode Hierachy         0.6795425   510
# LOINC Hierachy (LP1)     0.5041028   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_0.001_lambda_0_0618.Rdata")
fit = svd(B)
idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
rownames(U) = namelist
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
#                                auc   num
# RxNorm High Level        0.5302954 20075
# related RxNorm from UMLS 0.7302132  1107
# RxNorm Hierachy (rm)     0.6810621   703
# RxNorm Hierachy          0.7353306   660
# PheCode Hierachy         0.7595079   510
# LOINC Hierachy (LP1)     0.5527471   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_1e-04_lambda_0_0618.Rdata")
fit = svd(B)
idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
rownames(U) = namelist
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
#                                auc   num
# RxNorm High Level        0.4940409 20075
# related RxNorm from UMLS 0.5094402  1107
# RxNorm Hierachy (rm)     0.6198673   703
# RxNorm Hierachy          0.6339486   660
# PheCode Hierachy         0.5285160   510
# LOINC Hierachy (LP1)     0.4801014   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_tau_1e-05_lambda_0_0618.Rdata")
fit = svd(B)
idx = which(sign(fit$u[1,])==sign(fit$v[1,]))
U = fit$u[,idx] %*% diag(sqrt(fit$d[idx]))
rownames(U) = namelist
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
#                                auc   num
# RxNorm High Level        0.4835499 20075
# related RxNorm from UMLS 0.3984858  1107
# RxNorm Hierachy (rm)     0.5942061   703
# RxNorm Hierachy          0.5970363   660
# PheCode Hierachy         0.4658824   510
# LOINC Hierachy (LP1)     0.4174457   477

load("~/Ising_model/result/2306/MS_AD_Theta_est_convex_eta_0.1_0618.Rdata")
rownames(est_Theta$U) = namelist
AUC = get_AUC_embed(est_Theta$U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
# auc   num
# RxNorm High Level        0.5203555 20075
# related RxNorm from UMLS 0.6552961  1107
# RxNorm Hierachy (rm)     0.6155250   703
# RxNorm Hierachy          0.6594949   660
# PheCode Hierachy         0.6823260   510
# LOINC Hierachy (LP1)     0.5678529   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_convex_eta_0.01_0618.Rdata")
rownames(est_Theta$U) = namelist
AUC = get_AUC_embed(est_Theta$U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
# auc   num
# RxNorm High Level        0.5346549 20075
# related RxNorm from UMLS 0.7201352  1107
# RxNorm Hierachy (rm)     0.6431328   703
# RxNorm Hierachy          0.7039325   660
# PheCode Hierachy         0.7191888   510
# LOINC Hierachy (LP1)     0.5797283   477
load("~/Ising_model/result/2306/MS_AD_Theta_est_convex_eta_0.1_0618_1.Rdata")
rownames(U) = namelist
AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC
# auc   num
# RxNorm High Level        0.5368381 20075
# related RxNorm from UMLS 0.6924213  1107
# RxNorm Hierachy (rm)     0.6248348   703
# RxNorm Hierachy          0.6997314   660
# PheCode Hierachy         0.6605113   510
# LOINC Hierachy (LP1)     0.5752981   477

## filtered pairs AUC ####

load("~/Ising_model/result/2306/MS_AD_codified_pairs.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
source("~/Ising_model/code/2305/get_AUC_MS.R")

load("~/Ising_model/result/2306/MS_AD_sample_full_matrix.Rdata")
namelist = colnames(sample_matrix)
rm(sample_matrix)

library(readxl)
MS_Drugs <- read_excel("/n/data1/hsph/biostat/celehs/lab/zig728/Ising_Model/MS Drugs Mechanism Citations Efficacy 20230521.xlsx")
related_code = unique(paste("RXNORM:",MS_Drugs$`RxNorm Ingredient id`,sep=""))
rm(MS_Drugs)
related_code = c(related_code, "PheCode:335", "PheCode:290.11")

add_pairs = data.frame(code1 = related_code[1:36],
                       code2 = "PheCode:335",
                       source = "MS_treat",
                       type = "related",
                       group = "codi-codi",
                       nullcode1 = "",
                       nullcode2 = "")
pairs = rbind(pairs, add_pairs)

more_related_code = pairs %>%
  filter(source %in% c("related RxNorm from UMLS","RxNorm Hierachy","PheCode Hierachy")) %>%
  filter(code1 %in% related_code | code2 %in% related_code) %>%
  select(code1, code2)
more_related_code = unique(c(more_related_code$code1,more_related_code$code2))
more_related_code = unique(c(more_related_code,related_code))

subpairs1 = pairs %>%
  filter(code1 %in% related_code | code2 %in% related_code)

subpairs2 = pairs %>%
  filter(code1 %in% more_related_code | code2 %in% more_related_code)

load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_2.Rdata")
rownames(U) = namelist

AUC = get_AUC_embed(U, pairs, WithNull = FALSE,
                    dict = dict, normalize = FALSE)
AUC

AUC1 = get_AUC_embed(U, subpairs1, WithNull = FALSE,
                     dict = dict, normalize = FALSE)
AUC1

AUC1 = get_AUC_embed(U/apply(U,1,norm,'2'), 
                     get_merge_pairs(subpairs1,
                                     merged = list(RxNorm_Hierachy = c("RxNorm Hierachy (rm)","RxNorm Hierachy"))), 
                     WithNull = FALSE,
                     dict = dict, normalize = FALSE)
AUC1
AUC1[[1]] = AUC1[[1]][2,,drop = FALSE]
AUC1[[2]] = AUC1[[2]][c(3,6,10),,drop = FALSE]

AUC2 = get_AUC_embed(U, subpairs2, WithNull = FALSE,
                     dict = dict, normalize = FALSE)
AUC2

load("/n/data1/hsph/biostat/celehs/lab/zig728/Ising_Model/PairsWithNull_MS.Rdata")
pairsadd = data.frame(code1 = related_code[1:36],
                      code2 = "PheCode:335",
                      source = "MS related",
                      type = "related",
                      group = "codi-codi",
                      nullcode1 = "",
                      nullcode2 = "")
newpairs = rbind(pairs, add_pairs)

np = get_merge_pairs(newpairs, merged = list(all = setdiff(unique(pairs$source),"MS related"),
                                             ddx = c('ddx','possibly_equivalent_to','same_as'),
                                             ms_related = c('MS related')))
get_AUC_embed(U, np %>% filter(code1%in%related_code | code2%in%related_code), WithNull = FALSE,
              dict = dict, normalize = FALSE)
