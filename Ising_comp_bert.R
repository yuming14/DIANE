## AUC compare ####

load("~/Ising_model/result/2306/MS_AD_codified_pairs.Rdata")
# load("~/Ising_model/result/2306/MS_AD_codified_dict.Rdata")
load("~/Ising_model/result/2306/MS_AD_codified_selected_dict.Rdata")
source("~/Ising_model/code/2305/get_AUC_MS.R")
load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_15.Rdata")
dict = dict[which(!is.na(dict$desc)),]
pairs_new = pairs %>%
  filter(code1!=code2)
pairs = get_nullpairs(pairs_new, dict)
AUC = get_AUC_embed(U, pairs, WithNull = TRUE,
                    dict = dict, normalize = FALSE)

library(readr)
AUClist = lapply(c("bert","biobert","pubmedbert","sapbert"), function(x){
  bert <- read_csv(paste("/n/data1/hsph/biostat/celehs/lab/zig728/PLM/data/processed/embeddings/",x,
                         "_MEAN_AD_MS_code_desc_0803.csv",sep=""),
                   col_names = FALSE)
  namelist <- bert$X1
  bert <- as.matrix(bert[,-1])
  rownames(bert) <- namelist
  return(get_AUC_embed(bert, pairs, WithNull = TRUE,
                       dict = dict, normalize = TRUE))
})

AUCcomp <- lapply(AUClist, function(x) do.call("rbind",x))
AUCcomp <- do.call("cbind", AUCcomp)
AUCcomp <- AUCcomp[,c(2,1,3,5,7)]
AUCcomp <- cbind(do.call("rbind",AUC)[,1], AUCcomp)
colnames(AUCcomp) = c("Ising","num","bert","biobert","pubmedbert","sapbert")

library(knitr)
idx = match(c("codi-codi(similar).RxNorm Hierachy",
              "codi-codi(similar).PheCode Hierachy",
              "codi-codi(related).ddx",
              "codi-codi(related).may_treat",
              "codi-codi(related).classifies",
              "codi-codi(related).may_prevent"),
            rownames(AUCcomp))
kable(AUCcomp[idx,c(2,1,3,4,5,6)], "latex", 3)

## risk prediction ####
library(readr)
library(MASS)
library(tidymodels)
library(tidyr)
library(pROC)
library(survival)
library(survminer)
load("~/Ising_model/result/2306/Patient_code_counts_AD_MS.Rdata")
pi = list(AD = AD_p_count$patient_num,
          MS = MS_p_count$patient_num)
AD_p_count <- as.matrix(AD_p_count[,-1])
rownames(AD_p_count) <- pi$AD
MS_p_count <- as.matrix(MS_p_count[,-1])
rownames(MS_p_count) <- pi$MS

idx = intersect(rownames(AD_p_count), rownames(MS_p_count))
ADMS_p_count <- MS_p_count[match(idx,rownames(MS_p_count)),]
AD_p_count <- AD_p_count[which(!rownames(AD_p_count)%in%idx),]
MS_p_count <- MS_p_count[which(!rownames(MS_p_count)%in%idx),]

codelist <- read_csv("/n/data1/hsph/biostat/celehs/lab/zig728/PLM/data/processed/embeddings/bert_MEAN_AD_MS_code_desc_0803.csv",
                     col_names = FALSE)
codelist <- codelist$X1
codelist <- setdiff(codelist, c("PheCode:290.1","PheCode:335"))

train_index <- list(AD = sample(rownames(AD_p_count), nrow(AD_p_count)/2),
                    MS = sample(rownames(MS_p_count), nrow(MS_p_count)/2))
test_index <- list(AD = setdiff(rownames(AD_p_count), train_index$AD),
                   MS = setdiff(rownames(MS_p_count), train_index$MS))

get_patient_embedding = function(embed, p_count, codelist){
  embed <- as.matrix(embed[match(codelist, rownames(embed)),])
  p_count <- p_count[,match(codelist, colnames(p_count))]
  p_count[which(is.na(p_count))] <- 0
  p_count <- p_count/rowSums(p_count)
  return(p_count%*%embed)
}

get_classify_tree = function(AD, MS, train_index, test_index){
  AD_train <- as.data.frame(AD[match(train_index$AD,rownames(AD)),])
  AD_test <- as.data.frame(AD[match(test_index$AD,rownames(AD)),])
  rm(AD)
  MS_train <- as.data.frame(MS[match(train_index$MS,rownames(MS)),])
  MS_test <- as.data.frame(MS[match(test_index$MS,rownames(MS)),])
  rm(MS)
  tree_spec <- decision_tree() %>%
    set_engine("rpart") %>%
    set_mode("regression")
  colnames(AD_train) <- colnames(MS_train) <- 
    colnames(AD_test) <- colnames(MS_test) <- paste("X",1:ncol(AD_train),sep="")
  AD_train$disease <- 1
  AD_test$disease <- 1
  MS_train$disease <- 0
  MS_test$disease <- 0
  train_data <- as.data.frame(rbind(AD_train, MS_train))
  test_data <- as.data.frame(rbind(AD_test, MS_test))
  tree_fit <- tree_spec %>%
    fit(disease ~ ., data = train_data)
  predictions <- tree_fit %>%
    predict(test_data) %>%
    pull(.pred)
  ans = c(auc = roc(test_data$disease, predictions, direction = "<")$auc,
          brier = mean((test_data$disease - predictions)^2))
  cat(ans)
  return(ans)
}

accu_comp = lapply(c("Ising","bert","biobert","pubmedbert","sapbert"), function(x){
  if(x=="Ising"){
    load("~/Ising_model/result/2306/MS_AD_Theta_est_eta_0.001_0620_15.Rdata")
    bert <- as.matrix(U)
    rm(Theta_new, U, V, delta, t)
  }else{
    bert <- read_csv(paste("/n/data1/hsph/biostat/celehs/lab/zig728/PLM/data/processed/embeddings/",x,
                           "_MEAN_AD_MS_code_desc_0803.csv",sep=""),
                     col_names = FALSE)
    namelist <- bert$X1
    bert <- as.matrix(bert[,-1])
    rownames(bert) <- namelist
  }
  AD <- get_patient_embedding(embed = bert, p_count = AD_p_count, 
                              codelist = codelist)
  MS <- get_patient_embedding(embed = bert, p_count = MS_p_count, 
                              codelist = codelist)
  return(get_classify_tree(AD, MS, train_index, test_index))
})
