require(stringr)
require(e1071)
require(ggbiplot)
require(adegenet)
require(dplyr)

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

# Set up data

setwd("~/Documents/PhD/Interaction/pspA/")

design_matrix <- dplyr::as_data_frame(read.delim("design_matrix.txt", 
                                                 header=T, stringsAsFactors = F))
# For datasets without spades
#keep_cols = which(unlist(lapply(colnames(design_matrix), function (x) is.na(str_match(x, "spades")))))

training <- filter(design_matrix, !is.na(pspA.allele))
rownames(training) <- unlist(training[,"sample.name"])
training <- imputeData(training)

test <- filter(design_matrix, is.na(pspA.allele))
test <- imputeData(test)
rownames(test) <- unlist(test[,"sample.name"])

####
### pspA
####

# Cut into predictors (same for both models) and labels (different for each model)
pspA_train <- dplyr::select(training, -sample.name)
rownames(pspA_train) <- unlist(training[,"sample.name"])

# PCA
pspA.pca <- prcomp(dplyr::select(pspA_train,-pspA.allele), center = T, scale. = T)
ggbiplot(pspA.pca, obs.scale = 1, var.scale = 1, 
         groups = as.factor(pspA_train$pspA.allele), ellipse = TRUE, 
         circle = TRUE,var.axes = F)

# SVM
pspA_svm <- svm(as.factor(pspA.allele) ~ ., data=pspA_train, kernel="linear",cost=10)
pspA_predict <- predict(pspA_svm, newdata = test)

write.table(pspA_predict, file="pspA_predictions.txt", 
            quote = F, sep="\t", row.names=T, col.names=F)

# Features (for 1 vs 2)
sort((t(pspA_svm$coefs) %*% pspA_svm$SV)[1,]^2)
