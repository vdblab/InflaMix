rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
load("output/dev_inf_and_coxpredmodels.RData")
source("scripts/functions_constants.R")

library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)
library(survcomp)
library(pROC)

# Select all patients not in the derivation cohort
df_all_chrt <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  # To evaluate InflaMix validation with just six labs (albumin, hgb, crp, ldh, alp, ast) -
  # uncomment the line below, which ignores all other labs (This will give figures 3d-l)
  #mutate_at(vars(starts_with(focus_labs)), ~.*NA) %>%
  filter(cohort!="MSK Development") %>%
  mutate(record_id=record_id)

# Master data frame of all validation cohort data
val_df <- reformat_table_tp(df_all_chrt, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  mutate(everCR_100=factor(everCR_100, levels=c("CR", "Not CR"))) %>%
  mutate(crs24=factor(crs24, levels=c("CRS 0-1", "CRS >1"))) %>%
  mutate(icans24=factor(icans24, levels=c("ICANS 0-1", "ICANS >1"))) %>%
  rename(Cluster=cluster) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref.factor)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  mutate(preld_crpuln_ratio=noLCS_preld_crp/uln_preld_crp) %>%
  mutate(preld_ferruln_ratio=noLCS_preld_ferritin/uln_preld_ferritin) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# For each of the following validation analyses, define the cohort, calculate inference metrics
# Note that the functions will evaluate various inferences that are not relevant or of interest.
# Table 2 of the paper clarifies the relevant analyses.

#### II. MSK LBCL Validation
coh2_df <- val_df %>% filter(cohort=="MSK Validation")
coh2_metrics <- clust_metrics(df=coh2_df, dxres=lymphoma_dxres,
              covar=lbcl_covar, lcovar=lbcl_covar, coxmodels=coxmodels_crpferr,
              res_covar=lbcl_covar, wght=TRUE, aname="msk_validation")

# Calculate AUCs for no CR- need to consider probability of inflammatory cluster assignment as the covariate for these analyses
dev_df_hc_abc <- coh2_df %>% mutate(cluster=Cluster) %>%
  filter(cluster=="Inflammatory") %>%
  group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
  mutate(record_id=record_id)

dev_df_pred <- dev_df_hc_abc
# AUC of univariate LR model predicting no CR considering only clustering
concordance.index(x=dev_df_pred$tau, cl=dev_df_pred$everCR_100)

set.seed(18)
# AUC of multivariate LR model with CRP and ferritin as additional covariates without considering cluster
predicted0 <- predict.glm(coxmodels_crpferr$crpferr_cox0nocr, newdata = dev_df_pred, type="response")
roc02 <- roc(dev_df_pred$everCR_100 ~ predicted0)
# AUC of multivariate LR model with CRP and ferritin as additional covariates considering cluster
predicted <- predict.glm(coxmodels_crpferr$crpferr_coxnocr, newdata = dev_df_pred, type="response")
rpc2 <- roc(dev_df_pred$everCR_100 ~ predicted)
# 95% CI for multivariate AUC considering cluster
ci(rpc2)
# Does addition of clustering significantly improve AUC of the model over other covariates?
roc.test(roc02, rpc2, method="bootstrap")$p.value #gives p-value


##### III. SMC+HMH LBCL Validation
coh3_df <- val_df %>% filter(cohort=="Center Validation")
coh3_metrics <- clust_metrics(df=coh3_df, dxres=lymphoma_dxres,
                                 covar=lbcl_covar, lcovar=lbcl_covar,coxmodels=coxmodels_crpferr,
                                 res_covar=lbcl_covar, wght=TRUE, aname="center_validation")

# Calculate AUCs for no CR- need to consider probability of inflammatory cluster assignment as the covariate for these analyses
dev_df_hc_def <- coh3_df %>% mutate(cluster=Cluster )%>%
  filter(cluster=="Inflammatory") %>%
  group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
  mutate(record_id=record_id)

dev_df_pred <- dev_df_hc_def
# AUC of univariate LR model predicting no CR considering only clustering
concordance.index(x=dev_df_pred$tau, cl=dev_df_pred$everCR_100)

dev_df_pred <- dev_df_hc_abc
concordance.index(x=dev_df_pred$tau, cl=dev_df_pred$everCR_100)

set.seed(18)
# AUC of multivariate LR model without considering cluster
predicted0 <- predict.glm(coxmodels_crpferr$crp_cox0nocr, newdata = dev_df_pred, type="response")
roc0 <- roc(dev_df_pred$everCR_100 ~ predicted0)
# AUC of multivariate LR model considering cluster
predicted <- predict.glm(coxmodels_crpferr$crp_coxnocr, newdata = dev_df_pred, type="response")
rpc <- roc(dev_df_pred$everCR_100 ~ predicted)
# 95% CI for multivariate AUC considering cluster
ci(rpc)
# Does addition of clustering significantly improve AUC of the model over other covariates?
roc.test(roc0, rpc, method='bootstrap')$p.value #gives p-value


#### Cohort IV. MCL + FL Validation
coh4_df <- val_df %>% filter(cohort=="Non-Hodgkin Lymphoma Validation")

coh4_metrics <- clust_metrics(df=coh4_df, dxres=lymphoma_dxres,
                                 covar=lymphoma_covar, lcovar=lbcl_covar,coxmodels=coxmodels_crpferr,
                                 res_covar=lymphoma_covar, wght=TRUE, aname="nhl_validation")

##############################
# Exploratory Analyses by costim domain
df_all_chrt_plusdev <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1)

all_chrt_df <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  mutate(everCR_100=factor(everCR_100, levels=c("CR", "Not CR"))) %>%
  mutate(crs24=factor(crs24, levels=c("CRS 0-1", "CRS >1"))) %>%
  mutate(icans24=factor(icans24, levels=c("ICANS 0-1", "ICANS >1"))) %>%
  rename(Cluster=cluster) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref.factor)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# All patients treated with CD28 CAR-T
cd28_df <- all_chrt_df %>% filter(costim=="CD28")
cd28_metrics <- clust_metrics(df=cd28_df, dxres=lymphoma_dxres,
                                 covar=costim_covar, lcovar=lbcl_covar,coxmodels=coxmodels_crpferr,
                                 res_covar=costim_covar, wght=TRUE, aname="cd28_validation")

# All patients treated with 41BB CAR-T
bb41_df <- all_chrt_df %>% filter(costim=="41BB")
bb41_metrics <- clust_metrics(df=bb41_df, dxres=lymphoma_dxres,
                                 covar=costim_covar, lcovar=lbcl_covar,coxmodels=coxmodels_crpferr,
                                 res_covar=costim_covar, wght=TRUE, aname="41BB_validation")

validation_inf <- bind_rows(
  coh2_metrics,
  coh3_metrics,
  coh4_metrics,
  cd28_metrics,
  bb41_metrics
) %>% mutate(mm_cindex=round(mm_cindex,2))


# Use the below code to plot PFS and OS KM survival estimate curves
# Fig 2c, d, f, g, i, j
# Fig 3d, e, g, h, j, k
df_analysis = coh2_df
anval = "msk_validation"
otcm = "PFS"
tmax=25
qm = 5
sz1=23
pxm = 15
pym = 0.9
oxm = 7
oym = 0.13

ssr_survplot(df=df_analysis,timemax=tmax,qmonth=qm,wght=FALSE,sz=sz1,
                        xmi=pxm,
                        ymi=pym,
                        event=tolower(otcm),
                        metric=validation_inf %>% filter(analysis==anval & outcome==otcm),
                        labl="")


# Use the below code to plot CR bar plots
# Fig 2e, h, k
# Fig 3f, i, l
m1 = validation_inf %>% filter(analysis==anval & outcome=="No CR")
sz1=12
xm = 1
ym = 0.98
labl = ""

plot_grid((clust_bar_plot(df=df_analysis, tp_pre="", list_res=c("CR"), sz=sz1/1.8, xmi=xm, ymi=ym, metric=m1) + ylim(0, 1.1) ), labels=c(labl), scale=1, label_size = 30)
