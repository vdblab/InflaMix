rm(list=ls())
load(bstfun::here_data("output/model.RData"))
load(bstfun::here_data("output/scaledata.RData"))
source("scripts/analysis/functions_constants.R")

library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)

lbcl_covar="cluster+age+primary_ref.factor+costim+bl_mtv"

df_all_chrt <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>% filter(dx_simple.factor=="Large B-cell Lymphoma")

dev_df <- reformat_table_tp(df_all_chrt, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  mutate(everCR_100=factor(everCR_100, levels=c("CR", "Not CR"))) %>%
  mutate(crs24=factor(crs24, levels=c("CRS 0-1", "CRS >1"))) %>%
  mutate(icans24=factor(icans24, levels=c("ICANS 0-1", "ICANS >1"))) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref.factor)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(bl_mtv)) %>%
  mutate(record_id=record_id)

cr_OR_mean <- tidy(glm(as.formula(paste0("everCR_100~", lbcl_covar)), weight=tau, data=dev_df, family=binomial())) %>%
  filter(term=="clusterInflammatory")
pfs_HR_mean  <- tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~", lbcl_covar)), data = dev_df, weights=dev_df$tau))%>%
  filter(term=="clusterInflammatory")
os_HR_mean  <- tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~", lbcl_covar)), data = dev_df, weights=dev_df$tau))%>%
  filter(term=="clusterInflammatory")

