require(stringr)
require(e1071)
require(ggbiplot)
require(adegenet)
require(dplyr)

# Functions

imputeData <- function(x)
{
  # Remove non-informative rows (all same value)
  x <- x[,!apply(apply(x, 2, duplicated), 2, sum) == nrow(x) - 1]
  
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

# Set up data

setwd("~/Documents/PhD/Interaction/zmpA/")

design_matrix <- dplyr::as_data_frame(read.delim("design_matrix.txt", 
                                                 header=T, stringsAsFactors = F))

training <- filter(design_matrix, !is.na(zmpA.allele))
rownames(training) <- unlist(training[,"sample.name"])
training <- imputeData(training)
# After imputation, fields with only one informative value and the rest NA will be imputed the same
# Need to be removed before PCA
training <- training[,!apply(apply(training, 2, duplicated), 2, sum) == nrow(training) - 1]

test <- filter(design_matrix, is.na(zmpA.allele))
zmp_allele <- test$zmpA.allele

test <- imputeData(test) # Lose the zmpA.allele col as all NA
test <- test[,!apply(apply(test, 2, duplicated), 2, sum) == nrow(test) - 1]
test <- dplyr::bind_cols(data.frame(sample.name=test$sample.name, stringsAsFactors = F),
                          data.frame(zmpA.allele=zmp_allele), test[,c(-1,-2)])

# Some columns are removed due to all being the same
keep_columns <- intersect(colnames(test), colnames(training))
training <- training[,intersect(colnames(test), colnames(training))]
test <- test[,intersect(colnames(test), colnames(training))]
rownames(test) <- unlist(test[,"sample.name"])

####
### zmpA
####

# Cut into predictors (same for both models) and labels (different for each model)
zmpA_train <- dplyr::select(training, -sample.name)
rownames(zmpA_train) <- unlist(training[,"sample.name"])

# PCA
zmpA.pca <- prcomp(dplyr::select(zmpA_train,-zmpA.allele), center = T, scale. = T)
ggbiplot(zmpA.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(zmpA_train$zmpA.allele), ellipse = TRUE, 
         circle = TRUE,var.axes = F)

# SVM
zmpA_svm <- svm(as.factor(zmpA.allele) ~ ., data=zmpA_train, kernel="linear",cost=10)
zmpA_predict <- predict(zmpA_svm, newdata = test)

write.table(zmpA_predict, file="zmpA_predictions.txt", 
            quote = F, sep="\t", row.names=T, col.names=F)



