require(fmsb)
require(glmnet)

setwd("~/Documents/PhD/dutch_carriage/charge/")

sero_charge <- read.delim("~/Documents/PhD/dutch_carriage/charge/sero_charge.covar", header=FALSE, stringsAsFactors=FALSE)
lane_serotypes <- read.delim("~/Documents/PhD/dutch_carriage/charge/lane_serotypes.txt", header=FALSE, stringsAsFactors=FALSE)
for_R_invasive <- read.delim("~/Documents/PhD/dutch_carriage/charge/for_R_invasive.pheno", header=FALSE, stringsAsFactors=FALSE)

charge_data <- merge(sero_charge, for_R_invasive, by.x = "V1", by.y="V1")
charge_data <- merge(charge_data, lane_serotypes, by.x = "V1", by.y="V1", all.x = T, all.y = T)

rownames(charge_data) <- charge_data$V1
charge_data <- charge_data[,c(3,6,4)]
colnames(charge_data) <- c("charge", "serotype", "invasive")

charge_glm <- glm(factor(invasive) ~ charge, data = charge_data, family = binomial())
summary(charge_glm)
#charge       0.12392    0.01377   9.001   <2e-16 ***
NagelkerkeR2(charge_glm)
#$N
#[1] 1893
#
#$R2
#[1] 0.06184394

charge_glm <- glm(factor(invasive) ~ charge, data = charge_data[charge_data$charge!=-12.710,], family = binomial())
NagelkerkeR2(charge_glm)
#$N
#[1] 1458
#
#$R2
#[1] 0.07964819

lasso_X <- as.matrix(model.matrix(factor(invasive) ~ serotype, charge_data)[,-1])
lasso_y = charge_data[!is.na(charge_data$serotype),"invasive"]

lasso.sero <- glmnet(lasso_X, lasso_y, alpha = 1,
                     family = 'binomial')
plot(lasso.sero,xvar="lambda",label=T)
cv.lasso.sero <- cv.glmnet(lasso_X, lasso_y, family = 'binomial',
                           alpha = 1, nfolds = length(lasso_y))
plot(cv.lasso.sero)
coef(cv.lasso.sero, s="lambda.1se")

# use selected predictors in regression
selected <- lasso_X[,which(coef(cv.lasso.sero, s="lambda.1se")[-1] != 0)]
glm.sero <- glm(lasso_y ~ selected, family = binomial())
summary(glm.sero)
NagelkerkeR2(glm.sero)

#$N
#[1] 1735
#
#$R2
#[1] 0.4493011

# without selection
#NagelkerkeR2(glm(lasso_y ~ lasso_X, family = binomial()))
#$N
#[1] 1735
#
#$R2
#[1] 0.459293
