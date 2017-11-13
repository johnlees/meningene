require(snpStats)
require(Subtest)

library(foreach)
library(doParallel)

cores=detectCores()
cl <- makeCluster(cores[1])
registerDoParallel(cl)

# Wish to test
# 1) Sp vs. Nm meningitis, using NL data + controls
# 2) Sp meningitis vs. bacteremia, using DK data + controls

subtest_p = function(snp_genotypes, ya, yd, both_covariates=NULL, ldak_weights,
                     start_pars=c(0.8, 0.1, 2, 2, 3, 0.5))
{
  rand_samples = 1500
  sample_size = 400
  
  # Generate z scores
  z_mat <- matrix()
  if(is.null(both_covariates))
  {
    z_mat <- z_scores(snp_genotypes, 
                          ya, 
                          yd,
                          signed = TRUE, control = TRUE)
  }
  else
  {
    z_mat <- z_scores(snp_genotypes, 
                      ya, 
                      yd, 
                      both_covariates, 
                      both_covariates, 
                      signed = TRUE, control = TRUE)
  }
  good_idx <- which(!is.nan(z_mat[,1]) & !is.nan(z_mat[,2]))
  
  # Fit 3 Gaussian mixture to alt and null models
  print("Fitting model")
  alt_fit <- fit.3g(z_mat[good_idx,], pars = start_pars, weights = ldak_weights[good_idx],
                      C = 1, fit_null = FALSE, maxit = 10000, tol = 1e-04, sgm = 0.8,
                      one_way = FALSE, syscov = 0, accel = TRUE, verbose = TRUE,
                      file = NULL, n_save = 20, incl_z = TRUE, em = TRUE,
                      control = list(factr = 10))
  
  print("Fitting null model")
  null_fit <- fit.3g(z_mat[good_idx,], pars = start_pars, weights = ldak_weights[good_idx],
                    C = 1, fit_null = TRUE, maxit = 10000, tol = 1e-04, sgm = 0.8,
                    one_way = FALSE, syscov = 0, accel = TRUE, verbose = TRUE,
                    file = NULL, n_save = 20, incl_z = TRUE, em = TRUE,
                    control = list(factr = 10))
  
  # Calculate PLR statistic from fit likelihoods
  plr <- alt_fit$logl - null_fit$logl - (alt_fit$logl_a - null_fit$logl_a)
  print(paste0("plr: ", plr))

  # Use subsampling to estimate null distribution of PLR
  # Generate starting parameter vector  - s1 and pi1 are fixed
  print("Generate starting parameters")
  ff = fit.em(z_mat[good_idx,1], weights=ldak_weights[good_idx], pi0_init=0.999)
  s1 = ff$sigma
  pi1 = 1 - ff$pi0; 
  pi0 = 1-(2*pi1)
  parsx=c(pi0,pi1,1,s1,1,0)
  print(parsx)
  
  # Create subsamples of the case group
  cores=detectCores()
  sims = matrix(data = NA, nrow = rand_samples, ncol = 18)
  case_idxs = which(ya == 1)
  sims <- foreach(it=1:rand_samples, .combine = rbind,
                  .packages=c('snpStats','Subtest')) %dopar%
  {
    if (it %% cores[1] == 0)
    {
      print(paste0("Subsample ", it))
    }
    sample_idxs <- sample(case_idxs, sample_size, replace = F)
    # New zd scores
    if(is.null(both_covariates))
    {
      sample_z <- cbind(zd_scores(snp_genotypes[sample_idxs,good_idx], 
                                  yd[sample_idxs], 
                                  signed = T, control = T), 
                        z_mat[good_idx,1])
    }
    else
    {
      sample_z <- cbind(zd_scores(snp_genotypes[sample_idxs,good_idx], 
                                yd[sample_idxs], 
                                both_covariates[sample_idxs], 
                                signed = T, control = T), 
                      ya[sample_idxs])
    }
    
    # fit the full model C1
    C1 = fit.cond(sample_z, pars=parsx, fit_null=FALSE, 
                  weights = ldak_weights[good_idx])
    # fit the null model C0
    C0 = fit.cond(sample_z, pars=parsx, fit_null=TRUE,
                  weights = ldak_weights[good_idx])
    
    # summary vector
    vec = c(C0$logl, C1$logl, C1$logl_a-C0$logl_a, 0, 0, 1, C0$pars, C1$pars)
    vec
  }
  saveRDS(sims, file = "sims.Rdata")
  
  # Estimate the null PLR distribution from these simulations:
  print("Calculating null PLR distribution")
  S = rand_analysis(sims)
  # Compute the p-value:
  return(p_value(plr,S))
}

# Test 1 - NL

# QC'd matrix
snp_mat <- read.plink("qc/nl_strict_qc.bed", 
                      "qc/nl_strict_qc.bim", 
                      "qc/nl_strict_qc.fam")

# Read from LDAK
ldak_weights <- read.delim("sections/weights.all", header = T, sep = " ", stringsAsFactors = F)

# Read phenotypes, subtypes and covariates
Ya_nl <- read.delim("Ya.nl.txt", header = F, sep=" ", stringsAsFactors = F)
Yd_unfav <- read.delim("Yd.unfav.nl.txt", header = T, sep = " ", stringsAsFactors = F)
Yd_pneu <- read.delim("Yd.pneu_mening.nl.txt", header = T, sep = " ", stringsAsFactors = F)

Cov_nl <- read.delim("Ca.nl.txt", header = T, sep=" ", stringsAsFactors = F)

merge1 <- merge(Ya_nl, Yd_unfav, by.x="V1", by.y="FID", all.x = T)
merge2 <- merge(merge1, Yd_pneu, by.x="V1", by.y="FID", all.x = T)
merge_all <- merge(merge2, Cov_nl, by.x="V1", by.y="IMMUNOCOMPROMISED", all.x = T)

# Get p-values
unfav_p <- subtest_p(snp_mat$genotypes, merge_all$V2, merge_all$UNFAV_CASE_ONLY,
                     merge_all$IMMUNOCOMPROMISED, ldak_weights$Weight)
print("nl unfav p-value")
print(unfav_p)

pneu_mening_p <- subtest_p(snp_mat$genotypes, merge_all$V2, merge_all$PNEU_MENING,
                     merge_all$IMMUNOCOMPROMISED, ldak_weights$Weight)
print("nl pneu p-value")
print(pneu_mening_p)

# Test 2 - DK

# QC'd matrix
snp_mat <- read.plink("qc/dk_strict_qc.bed", 
                      "qc/dk_strict_qc.bim", 
                      "qc/dk_strict_qc.fam")

# Read from LDAK
ldak_weights <- read.delim("sections/weights.all", header = T, sep = " ", stringsAsFactors = F)

# Read phenotypes, subtypes and covariates
Ya_dk <- read.delim("Ya.dk.txt", header = F, sep=" ", stringsAsFactors = F)
colnames(Ya_dk) = c("IID", "Ya")
Yd_dk <- read.delim("Yd.dk.txt", header = T, sep = " ", stringsAsFactors = F)
colnames(Yd_dk) = c("IID", "Yd")

merge_all <- merge(Ya_dk, Yd_dk, by.x="IID", by.y="IID", all.x = T)

#
dk_p <- subtest_p(snp_genotypes=snp_mat$genotypes, ya=merge_all$Ya, yd=merge_all$Yd,
                  ldak_weights=ldak_weights$Weight, start_pars=c(8.68695e-01,1.30997e-01,5.34020e+00,1.00551e+00,5.09695e+00,2.54296e+01))
print("dk p-value")
print(dk_p)
