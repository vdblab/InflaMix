rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/000_functions_constants.R")

library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)
library(survcomp)
library(pROC)
library(labelled)
library(grid)
library(gridExtra)
library(cowplot)

# Covariates for regression models
lbcl_covar="cluster+age+primary_ref+costim+bin_preld_ldh"

# Derivation cohort only
df_all_chrt <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type=="Lymphoma Outcomes") %>%
  filter(cohort=="MSK Development")

# For regression analysis filter any patients with missing covariate data
dev_df <- reformat_table_tp(df_all_chrt, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(bin_preld_ldh)) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  droplevels() %>% 
  mutate(record_id=record_id)

# Estimates for Cox PH and Logistic regression
cr_OR_mean <- tidy(glm(as.formula(paste0("everCR_100~", lbcl_covar)), weight=tau, data=dev_df, family=binomial())) %>%
  filter(term=="clusterInflammatory") %>% select(estimate, pvalue=p.value) %>% as.double()
crs_OR_mean  <- tidy(glm(as.formula(paste0("crs24~", lbcl_covar)),weights = tau,data=dev_df,family=binomial())) %>%
  filter(term=="clusterInflammatory") %>% select(estimate, pvalue=p.value) %>% as.double()
icans_OR_mean  <- tidy(glm(as.formula(paste0("icans24~", lbcl_covar)),weights = tau,data=dev_df,family=binomial())) %>%
  filter(term=="clusterInflammatory") %>% select(estimate, pvalue=p.value) %>% as.double()
pfs_HR_mean  <- tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~", lbcl_covar)), data = dev_df, weights=dev_df$tau))%>%
  filter(term=="clusterInflammatory") %>% select(estimate, pvalue=p.value) %>% as.double()
os_HR_mean  <- tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~", lbcl_covar)), data = dev_df, weights=dev_df$tau))%>%
  filter(term=="clusterInflammatory") %>% select(estimate, pvalue=p.value) %>% as.double()

reformat.dat <- dev_df %>% group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup() %>%
  select(!c(cluster, tau))

ref_inflam_vec <- dev_mu[,which.max(dev_mu["d0_il6",])]
unit_ref_inflam_vec <- ref_inflam_vec / sqrt(sum(ref_inflam_vec^2))

ref_LEASTinflam_vec <- dev_mu[,which.min(dev_mu["d0_il6",])]
unit_ref_LEASTinflam_vec <- ref_LEASTinflam_vec / sqrt(sum(ref_LEASTinflam_vec^2))

# Bootstrap inferences across 100 iterations
i <- 1
k <- 1
n = 100
cr_OR_est <- 1:n
crs_OR_est <- 1:n
icans_OR_est <- 1:n
pfs_HR_est <- 1:n
os_HR_est <- 1:n

while(i <= n){
  set.seed(k)
  # Resampling with replacement
  resamples1 <- reformat.dat[sample(1:dim(reformat.dat)[1], size=dim(reformat.dat)[1], replace = TRUE),] %>% as.matrix() %>% as.data.frame()
  resamples <- resamples1 %>%
    mutate_at(vars(starts_with(c("age", "d0_",
                                 "tt_pfs_m", "ev_pfs", "tt_os_m", "ev_os"))), as.numeric)
  var_label(resamples) <- NULL
  resamp_d0clusters <- Mclust(resamples %>% select(starts_with("d0_")) %>% select(contains(cluster_labs)) %>% mutate_all(~as.numeric(.)), modelNames = "VVI", G = length(dev_pro))
  if(is.null(resamp_d0clusters)){
    k <- k+1
    next
  }

  # When we generate new clusters for each bootstrap, the cluster labels are not consistently carried over.
  # Identify which cluster is inflammatory by finding the cluster with a unit vector of mean laboratory values that is most parallel -
  # by calculating the dot product between the two - with the reference unit vector of our original model.
  dotprods <-  as.data.frame(resamp_d0clusters$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
    rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
    dplyr::group_by(cluster) %>% dplyr::summarize(dotprod = unit_ref_inflam_vec %*% meanvalues)
  LEASTdotprods <-  as.data.frame(resamp_d0clusters$parameters$mean) %>% transmute_all(function(x) { x / sqrt(sum(x^2)) } ) %>%
    rownames_to_column("labs") %>% pivot_longer(!labs, names_to = "cluster", values_to= "meanvalues") %>%
    group_by(cluster) %>% summarize(dotprod = unit_ref_LEASTinflam_vec %*% meanvalues)
  inflamm_index <- which.max(as.double(unlist(dotprods[,2])))
  LEASTinflamm_index <- which.max(as.double(unlist(LEASTdotprods[,2])))

  nclust <- resamp_d0clusters$G
  reformat.dat_temp <- do.call(rbind, replicate(nclust, resamples, simplify=FALSE))
  reformat.dat_temp$tau <- as.vector(resamp_d0clusters$z)
  nrows <- nrow(resamples)
  clustclass <- rep(1,nrows)
  for(j in 2:nclust){
    clustclass <- c(clustclass, rep(j,nrows))
  }
  reformat.dat_temp$class <- as.factor(clustclass)
  reformat.dat_temp <- reformat.dat_temp %>% dplyr::rename(cluster=class)

  reformat.dat_temp <- reformat.dat_temp %>%
    mutate(ev_os = as.numeric(ev_os)) %>%
    mutate(ev_pfs = as.numeric(ev_pfs)) %>%
    mutate(cluster=factor(ifelse(cluster==inflamm_index, "Inflammatory", ifelse(cluster==LEASTinflamm_index, "Non-Inflammatory", ifelse(length(dev_pro)>3, paste0("Neutral Cluster", cluster), "Neutral Cluster"))))) %>%
    mutate(cluster=factor(cluster, levels=rev(levels(cluster)))) %>%
    mutate(cluster=relevel(cluster, "Non-Inflammatory")) %>%
    mutate(everCR_100=factor(everCR_100, levels=c("CR", "Not in CR"))) %>%
    mutate(crs24=factor(crs24, levels=c("CRS 0-1", "CRS >1"))) %>%
    mutate(icans24=factor(icans24, levels=c("ICANS 0-1", "ICANS >1"))) %>%
    filter(tau!=0)

  # Inference estimates for bootstrapped model i
  cr_OR_est[i] <- tidy(glm(as.formula(paste0("everCR_100~", lbcl_covar)),weights = reformat.dat_temp$tau,data=reformat.dat_temp,family=binomial())) %>%
    filter(term=="clusterInflammatory") %>% select(estimate) %>% as.double()
  crs_OR_est[i]  <- tidy(glm(as.formula(paste0("crs24~", lbcl_covar)),weights = tau,data=reformat.dat_temp,family=binomial())) %>%
    filter(term=="clusterInflammatory") %>% select(estimate) %>% as.double()
  icans_OR_est[i]  <- tidy(glm(as.formula(paste0("icans24~", lbcl_covar)),weights = tau,data=reformat.dat_temp,family=binomial())) %>%
    filter(term=="clusterInflammatory") %>% select(estimate) %>% as.double()
  pfs_HR_est[i]  <- tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~", lbcl_covar)), data = reformat.dat_temp, weights=reformat.dat_temp$tau))%>%
    filter(term=="clusterInflammatory") %>% select(estimate) %>% as.double()
  os_HR_est[i]  <- tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~", lbcl_covar)), data = reformat.dat_temp, weights=reformat.dat_temp$tau))%>%
    filter(term=="clusterInflammatory") %>% select(estimate) %>% as.double()

    i <- i+1
  k <- k+1
}

# Convert 95% confidence intervals of bootstrapped estimates to HRs and ORs
cr_OR_CI <- exp(quantile(cr_OR_est, probs = c(.025, .975)))
crs_OR_CI <- exp(quantile(crs_OR_est, probs = c(.025, .975)))
icans_OR_CI <- exp(quantile(icans_OR_est, probs = c(.025, .975)))
pfs_HR_CI <- exp(quantile(pfs_HR_est, probs = c(.025, .975)))
os_HR_CI <- exp(quantile(os_HR_est, probs = c(.025, .975)))

# Calculate concordance for PFS/OS in our derivation cohort using probability of cluster assignment tau as predicted risk estimate (since clustering was unsupervised)
dev_df_inflam <- dev_df %>% filter(cluster=="Inflammatory")
mm_cindex_os <- round(survcomp::concordance.index(dev_df_inflam$tau, surv.time = dev_df_inflam$tt_os_m, surv.event = dev_df_inflam$ev_os)$c.index, 2)
mm_cindex_pfs <- round(survcomp::concordance.index(dev_df_inflam$tau, surv.time = dev_df_inflam$tt_pfs_m, surv.event = dev_df_inflam$ev_pfs)$c.index, 2)

# Combine all bootstrapped inference metrics into one master table. Placeholder NA values are there to facilitate merging with validation metrics table
bootstrap_all <-
as.data.frame(t(bind_cols(
as.data.frame(c(exp(cr_OR_mean[1]), cr_OR_mean[2], "No CR", cr_OR_CI, lbcl_covar, "devcohort_inf")),
as.data.frame(c(exp(crs_OR_mean[1]), crs_OR_mean[2], "CRS", crs_OR_CI, lbcl_covar, "devcohort_inf")),
as.data.frame(c(exp(icans_OR_mean[1]), icans_OR_mean[2], "ICANS", icans_OR_CI, lbcl_covar, "devcohort_inf")),
as.data.frame(c(exp(pfs_HR_mean[1]),pfs_HR_mean[2], "PFS", pfs_HR_CI, lbcl_covar,"devcohort_inf")),
as.data.frame(c(exp(os_HR_mean[1]), os_HR_mean[2], "OS",os_HR_CI, lbcl_covar, "devcohort_inf"))
))) %>% rownames_to_column("test") %>% select(!test)

colnames(bootstrap_all) <- c(
  "expEstimate", "pvalue", "outcome", "low_ci", "high_ci", "covariates",  "analysis"
)

bootstrap_inf <- bootstrap_all %>% mutate(expEstimate=round(as.numeric(expEstimate), 2), low_ci=round(as.numeric(low_ci), 2), high_ci=round(as.numeric(high_ci), 2), pvalue=signif(as.numeric(pvalue), 2))

# Binary outcome bar plots (CR, CRS, ICANS)
fig_df <- dev_df |>
  dplyr::rename(Cluster=cluster)

metric1 <- bootstrap_inf %>% filter(analysis=="devcohort_inf" & outcome=="CRS")
fig1h <- clust_bar_plot(df=fig_df, tp_pre="", list_res=c("CRS"), sz=5, metric=metric1, xmi=1, ymi=0.9)+
  ylim(0, 1.1)

metric2 <- bootstrap_inf %>% filter(analysis=="devcohort_inf" & outcome=="ICANS")
fig1i <- clust_bar_plot(df=fig_df, tp_pre="", list_res=c("ICANS"), sz=5, metric=metric2, xmi=1, ymi=0.9)+
  ylim(0, 1.1)

metric3 <- bootstrap_inf %>% filter(analysis=="devcohort_inf" & outcome=="No CR")
fig1j <- clust_bar_plot(df=fig_df, tp_pre="", list_res=c("CR"), sz=5, metric=metric3, xmi=1, ymi=0.9)+
  ylim(0, 1.1)


# KM survival estimate plots, can use this for figures 1k and 1l

fig1k <- clus_ggsurvplot(df=dev_df %>% mutate(Cluster=cluster),
             event="pfs",
             metric=bootstrap_inf %>% filter(analysis=="devcohort_inf" & outcome=="PFS"),
             wght=FALSE,
             sz=18,
             xmi=16,
             ymi=0.8,
             labadj=0.05,
             labl=""
             )

fig1l <- clus_ggsurvplot(df=dev_df %>% mutate(Cluster=cluster),
                      event="os",
                      metric=bootstrap_inf %>% filter(analysis=="devcohort_inf" & outcome=="OS"),
                      wght=FALSE,
                      sz=18,
                      xmi=16,
                      ymi=0.9,
                      labadj=0.01,
                      labl=""
                      )

