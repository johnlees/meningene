## NOTES
# See paper Gene 2002 Iannelli et. al.
# Distinguishing 1 and 5 is difficult
# 3 and 6 too, but not as bad
# Look at alignments in Fig 1 visually
# See also trees drawn from alignment

require(stringr)
require(randomForest)
require(e1071)
require(ggbiplot)
require(adegenet)
require(glmnet)
require(kknn)
require(EMCluster)
require(dplyr)
require(caret)

# Functions

imputeData <- function(x)
{
  for (i in 1:ncol(x))
  {
    # NA Fields which go to max
    if(grepl("indels",colnames(x)[i], fixed=T) || 
       grepl("mismatches",colnames(x)[i], fixed=T) || 
       grepl("truncated",colnames(x)[i], fixed=T) ||
       grepl("gaps",colnames(x)[i], fixed=T))
    {
      if (length(which(is.na(x[,i]))) > 0)
      {
        x[which(is.na(x[,i])),i] <- max(x[-which(is.na(x[,i])),i])  
      }
    }
    # NA Fields which go to min
    else if(grepl("coverage",colnames(x)[i], fixed=T) || 
            grepl("bitscore",colnames(x)[i], fixed=T) || 
            grepl("evalue",colnames(x)[i], fixed=T) ||
            grepl("length",colnames(x)[i], fixed=T) ||
            grepl("p.id",colnames(x)[i], fixed=T))
    {
      if (length(which(is.na(x[,i]))) > 0)
      {
        x[which(is.na(x[,i])),i] <- min(x[-which(is.na(x[,i])),i])    
      }
    }
    
    # Change cap of log values to maximum
    if(grepl("p.id",colnames(x)[i], fixed=T) || 
       grepl("coverage",colnames(x)[i], fixed=T) || 
       grepl("evalue",colnames(x)[i], fixed=T))
    {
      if (length(which(x[,i]==709)) > 0)
      {
        x[which(x[,i]==709),i] <- max(x[-which(x[,i]==709),i])
      }
    }
  }
  
  return(x)
}

# Useful variables
cbpA_alleles = c("pspC-1.1","pspC-2.2","pspC-2.3","pspC-2.4","pspC-2.5","pspC-3.10","pspC-3.11","pspC-3.12","pspC-3.13","pspC-3.1","pspC-3.5","pspC-3.6","pspC-3.7","pspC-3.8","pspC-3.9","pspC-4.1","pspC-5.1","pspC-5.2","pspC-6.10","pspC-6.11","pspC-6.12","pspC-6.13","pspC-6.14","pspC-6.1","pspC-6.2","pspC-6.3","pspC-6.4","pspC-6.5","pspC-6.6","pspC-6.7","pspC-6.8","pspC-6.9");
pspC_alleles = c("pspC-10.1","pspC-11.1","pspC-11.2","pspC-11.3", "pspC-7.1","pspC-7.2","pspC-7.3","pspC-7.4","pspC-8.1","pspC-8.2","pspC-8.3","pspC-8.4","pspC-9.1","pspC-9.2","pspC-9.3","pspC-9.4");

# Set up data

setwd("~/Documents/PhD/Interaction/pspC/")

design_matrix <- dplyr::as_data_frame(read.delim("design_matrix.txt", 
                                                 header=T, stringsAsFactors = F))

cbpA_cols <- which(unlist(lapply(colnames(design_matrix),
                                 function (x) str_match(x, "^pspC\\.(\\d+)\\.\\d+\\.")[,2]) < 7))
pspC_cols <- which(unlist(lapply(colnames(design_matrix),
                                 function (x) str_match(x, "^pspC\\.(\\d+)\\.\\d+\\.")[,2]) >= 7))

# For datasets without spades
#cbpA_cols <- which(unlist(lapply(colnames(design_matrix),
#                                 function (x) str_match(x, "^pspC\\.(\\d+)\\.\\d+\\.")[,2] < 7 & is.na(str_match(x, "spades")))))
#pspC_cols <- which(unlist(lapply(colnames(design_matrix),
#                                 function (x) str_match(x, "^pspC\\.(\\d+)\\.\\d+\\.")[,2] >= 7 & is.na(str_match(x, "spades")))))

training <- filter(design_matrix, !is.na(pspC.allele) & !is.na(cbpA.allele))
rownames(training) <- unlist(training[,"sample.name"])
training <- imputeData(training)
cbpA_training <- training[,c(1:3,cbpA_cols)]
pspC_training <- training[,c(1:3,pspC_cols)]

test <- filter(design_matrix, is.na(pspC.allele) | is.na(cbpA.allele))
test <- imputeData(test)

# All with labels
cbpA_all <- rbind(cbpA_training, dplyr::select(test[,c(1:3,cbpA_cols)]))
pspC_all <- rbind(pspC_training, dplyr::select(test[,c(1:3,pspC_cols)]))

# Test data only predictors, no labels
cbpA_test <- dplyr::select(test[,c(1:3,cbpA_cols)], -sample.name, -pspC.allele, -cbpA.allele)
rownames(cbpA_test) <- unlist(test[,"sample.name"])

pspC_test <- dplyr::select(test[,c(1:3,pspC_cols)], -sample.name, -pspC.allele, -cbpA.allele)
rownames(pspC_test) <- unlist(test[,"sample.name"])

# Manually labelled data to judge performance
manually_typed <- dplyr::as_data_frame(read.csv("~/Documents/PhD/Interaction/pspC/manually_typed.csv", stringsAsFactors=FALSE))

####
### cbpA
####

# Cut into predictors (same for both models) and labels (different for each model)
cbpA_train <- dplyr::select(cbpA_training, -sample.name, -pspC.allele)
rownames(cbpA_train) <- unlist(cbpA_training[,"sample.name"])

# PCA/DAPC
cbpA.pca <- prcomp(dplyr::select(cbpA_train,-cbpA.allele), center = T, scale. = T)
ggbiplot(cbpA.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(cbpA_train$cbpA.allele), ellipse = TRUE, 
         circle = TRUE,var.axes = F)

cbpA.dapc <- dapc(x = dplyr::select(cbpA_train,-cbpA.allele), grp=as.factor(cbpA_train$cbpA.allele), scale = T)
scatter(cbpA.dapc) # Good separation 
contrib <- loadingplot(cbpA.dapc$var.contr, axis=2,thres=0.007, lab.jitter=1)
predict.dapc(cbpA.dapc, cbpA_test)$assign

dapc_predict <- predict.dapc(cbpA.dapc, cbpA_test)$assign
names(dapc_predict) <- test$sample.name
confusionMatrix(dapc_predict[manually_typed$lane], as.factor(manually_typed$cbpA))

# Random forests
cbpA_forest <- randomForest(as.factor(cbpA.allele) ~ ., data=cbpA_train, importance=T)
#cbpA_predict <- predict(cbpA_forest) # Just checking it works - predicting from training data is otherwise uninformative
cbpA_predict <- predict(cbpA_forest, newdata = test)
names(cbpA_predict) <- test$sample.name
confusionMatrix(cbpA_predict[manually_typed$lane], as.factor(manually_typed$cbpA))

# SVM
cbpA_svm <- svm(as.factor(cbpA.allele) ~ ., data=cbpA_train, kernel="linear",cost=10)
cbpA_predict <- predict(cbpA_svm, newdata = cbpA_test)
confusionMatrix(cbpA_predict[manually_typed$lane], as.factor(manually_typed$cbpA))

svm.tune <- tune(svm, as.factor(cbpA.allele) ~ ., data=cbpA_train, tunecontrol = tune.control(cross=5),
                 ranges=list(cost=c(0.1,1,10,100,1000), gamma=c(0.5,1,2,3,4)))

# Nearest neighbours
cbpA_kknn<-kknn(as.factor(cbpA.allele) ~ ., cbpA_train, test)
cbpA_predict <- fitted(cbpA_kknn)
names(cbpA_predict) <- test$sample.name
confusionMatrix(cbpA_predict[manually_typed$lane], as.factor(manually_typed$cbpA))

# Gaussian mixture (semi-supervised)
semi_labels <- unlist(lapply(cbpA_all$cbpA.allele, function(x) if(!is.na(x)) {x+1} else {0}))
# Take 90% of variance from PCA
semi_pca <- prcomp(dplyr::select(cbpA_all, -sample.name, -pspC.allele, -cbpA.allele), center = T, scale. = T)
components <- which(cumsum(semi_pca$sdev^2/sum(semi_pca$sdev^2))<0.8)
projected <- t(semi_pca$rotation[,components]) %*% t(dplyr::select(cbpA_all, -sample.name, -pspC.allele, -cbpA.allele))

# slow when more dimensions are included
#init.EM(t(projected), nclass=7, lab=semi_labels)
cbpA.EM <- init.EM(t(projected)[,1:15], nclass=7, lab=semi_labels)
plotem(cbpA.EM, t(projected)[,1:2])

####
### pspC
####
pspC_train <- dplyr::select(pspC_training, -sample.name, -cbpA.allele)
rownames(cbpA_train) <- unlist(training[,"sample.name"])

pspC_svm <- svm(as.factor(pspC.allele) ~ ., data=pspC_train, kernel="linear",cost=10)
pspC_predict <- predict(pspC_svm, newdata = pspC_test)

semi_labels <- unlist(lapply(pspC_all$pspC.allele, function(x) if(is.na(x)) {0} else if (x==0) {1} else {x-6}))
components <- which(cumsum(pspC.pca$sdev^2/sum(pspC.pca$sdev^2))<0.8)
projected <- t(pspC.pca$rotation[,components]) %*% t(dplyr::select(pspC_all, -sample.name, -pspC.allele, -cbpA.allele))
#pspC.EM <- init.EM(t(projected)[,1:2], nclass=5, lab=semi_labels)

pspC.pca <- prcomp(dplyr::select(pspC_train,-pspC.allele), center = T, scale. = T)
ggbiplot(pspC.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(pspC_train$pspC.allele), ellipse = TRUE, 
         circle = TRUE,var.axes = F)

pspC.dapc <- dapc(x = dplyr::select(pspC_train,-pspC.allele), 
                  grp=as.factor(pspC_train$pspC.allele), scale = T)
predict.dapc(pspC.dapc, newdata=pspC_test)$assign

write.table(data.frame(cbpA_predict, pspC_predict), file="pspC_predictions.txt", 
            quote = F, sep="\t", row.names=T, col.names=F)



