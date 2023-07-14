rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/functions_constants.R")

library(ggalluvial)
library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)
library(ggsci)
library(ggrepel)
library(gridExtra)
library(ggpattern)
library(cowplot)
library(gtable)
library(gridGraphics)
library(ggplot2)

# Filter patients, analysis_type==1 signifies complete clinical metadata
df_all_chrt_plusdev <- df_all %>% left_join(df_labs_all %>% select(!cohort), by="record_id") %>%filter(analysis_type==1)

# Assign clusters using InflaMix at d0, pre-lymphodepletion, and pre-apheresis timepoints.
df_all_chrt_d0 <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  rename(d0_tau=tau)
df_all_chrt_preld <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preld_") %>%
  reframe(record_id, cluster=cluster, preld_tau=tau)
df_all_chrt_preaph <- reformat_table_tp(df_all_chrt_plusdev, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preaph_") %>%
  reframe(record_id, cluster=cluster, preaph_tau=tau)

# Merge cluster assignment probabilities together in wide format
df <- df_all_chrt_d0 %>%
  left_join(df_all_chrt_preaph, by=c("record_id", "cluster")) %>%
  left_join(df_all_chrt_preld, by=c("record_id", "cluster")) %>%
  mutate(ev_os=as.integer(ev_os), ev_pfs=as.integer(ev_pfs)) %>%
  mutate(everCR_100=factor(everCR_100, levels=c("CR", "Not CR"))) %>%
  filter(!is.na(age)) %>%
  filter(!is.na(primary_ref.factor)) %>%
  filter(!is.na(costim)) %>%
  filter(!is.na(noLCS_preld_ldh)) %>%
  filter(!is.na(kps90)) %>%
  mutate(bin_preld_ldh=factor(bin_preld_ldh)) %>%
  mutate(record_id=record_id)

# Restructuring data frame and selecting only the most probable cluster assignment with its associated probability
prealluv <- df %>% group_by(record_id) %>% dplyr::slice(which.max(preaph_tau)) %>% ungroup() %>%
  select(!contains("tau")) %>% mutate(preaph_cluster=cluster) %>% select(!cluster) %>%
  left_join(df %>% select(record_id, cluster, preld_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(preld_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(preld_cluster=cluster) %>% select(!cluster), by="record_id") %>%
  left_join(df %>% select(record_id, cluster, d0_tau) %>% group_by(record_id) %>%
              dplyr::slice(which.max(d0_tau)) %>% ungroup() %>% select(!contains("tau")) %>%
              mutate(d0_cluster=cluster) %>% select(!cluster), by="record_id") %>%
  column_to_rownames("record_id")

# Restruture the data into a format required for alluvial plotting
df_alluvial <- prealluv %>%
  group_by(
    preaph_cluster,
    preld_cluster,
    d0_cluster,
    bridge.factor
  ) %>% mutate(preaph_cluster=factor(preaph_cluster),
               preld_cluster=factor(preld_cluster),
               d0_cluster=factor(d0_cluster)) %>%
  summarize(Freq=n()) %>% ungroup() %>% drop_na()  %>% mutate(fills=factor(ifelse(preaph_cluster=="Inflammatory", 1, 0)))#%>% filter(Freq >=10)

df_alluvial1 <- df_alluvial %>%
  mutate(preaph_cluster=ifelse(preaph_cluster=="Inflammatory", "Infl.", "Non-Infl."),
         preld_cluster=ifelse(preld_cluster=="Inflammatory", "Infl.", "Non-Infl."),
         d0_cluster=ifelse(d0_cluster=="Inflammatory", "Infl.", "Non-Infl."))

# Plot Figure 6a - alluvial plot
pl_alluvial <- ggplot(as.data.frame(df_alluvial1),
                      aes(y = Freq,
                          axis1 = preaph_cluster,
                          axis2 = preld_cluster,
                          axis3 = d0_cluster
                      )) +
  guides(fill=FALSE)+
  geom_alluvium(aes(fill = bridge.factor), width = 1/12) +
  geom_stratum(width = 1/4, color = "black", aes(fill=fills)) +
  geom_label(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_x_discrete(limits = c("Pre-Apheresis",
                              "Pre-Lymphodepletion",
                              "Peri-Infusion"
  ), expand = c(.05, .05)) +
  theme_minimal()+
  scale_color_manual(values=c("#F39B7FFF", "#7E6148FF"))+
  scale_fill_manual(values=c("#1F77B4FF", "#FF7F0EFF", "#00A087FF", "#7E6148FF"))+
  ggtitle("") +
  labs(x = "",
       y = "Number Patients",
       title = ""
  )+
  theme(
    axis.text.x = element_text(angle = 0, size=13, face="bold"),
    axis.text.y = element_text(angle = 0, size=15),
    axis.title.x = element_text(angle = 0, size=18, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
  )

fig6a <-pl_alluvial

# Cluster heatmap in Figure 4b
library(ComplexHeatmap)
dfb1 <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  mutate(record_id=record_id)
dfb <- reformat_table_tp(dfb1, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preaph_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup()

dfc1 <- dfb %>% rename(preaph_cluster=cluster)
dfc2 <- dfb%>% rename(preaph_cluster=cluster)

names(dfc1)[names(dfc1) == 'preaph_albumin'] <- "preaph_Albumin"
names(dfc1)[names(dfc1) == 'preaph_wbc'] <- "preaph_WBC"
names(dfc1)[names(dfc1) == 'preaph_alk'] <- "preaph_ALP"
names(dfc1)[names(dfc1) == 'preaph_ast'] <- "preaph_AST"
names(dfc1)[names(dfc1) == 'preaph_crp'] <- "preaph_CRP"
names(dfc1)[names(dfc1) == 'preaph_ddimer'] <- "preaph_D-dimer"
names(dfc1)[names(dfc1) == 'preaph_ferritin'] <- "preaph_Ferritin"
names(dfc1)[names(dfc1) == 'preaph_hb'] <- "preaph_Hgb"
names(dfc1)[names(dfc1) == 'preaph_il10'] <- "preaph_IL-10"
names(dfc1)[names(dfc1) == 'preaph_il6'] <- "preaph_IL-6"
names(dfc1)[names(dfc1) == 'preaph_ldh'] <- "preaph_LDH"
names(dfc1)[names(dfc1) == 'preaph_plt'] <- "preaph_Plt"
names(dfc1)[names(dfc1) == 'preaph_tbr'] <- "preaph_Tbili"
names(dfc1)[names(dfc1) == 'preaph_tnfa'] <- "preaph_TNFa"

names(dfc2)[names(dfc2) == 'noLCS_preaph_albumin'] <- "noLCS_preaph_Albumin"
names(dfc2)[names(dfc2) == 'noLCS_preaph_wbc'] <- "noLCS_preaph_WBC"
names(dfc2)[names(dfc2) == 'noLCS_preaph_alk'] <- "noLCS_preaph_ALP"
names(dfc2)[names(dfc2) == 'noLCS_preaph_ast'] <- "noLCS_preaph_AST"
names(dfc2)[names(dfc2) == 'noLCS_preaph_crp'] <- "noLCS_preaph_CRP"
names(dfc2)[names(dfc2) == 'noLCS_preaph_ddimer'] <- "noLCS_preaph_D-dimer"
names(dfc2)[names(dfc2) == 'noLCS_preaph_ferritin'] <- "noLCS_preaph_Ferritin"
names(dfc2)[names(dfc2) == 'noLCS_preaph_hb'] <- "noLCS_preaph_Hgb"
names(dfc2)[names(dfc2) == 'noLCS_preaph_il10'] <- "noLCS_preaph_IL-10"
names(dfc2)[names(dfc2) == 'noLCS_preaph_il6'] <- "noLCS_preaph_IL-6"
names(dfc2)[names(dfc2) == 'noLCS_preaph_ldh'] <- "noLCS_preaph_LDH"
names(dfc2)[names(dfc2) == 'noLCS_preaph_plt'] <- "noLCS_preaph_Plts"
names(dfc2)[names(dfc2) == 'noLCS_preaph_tbr'] <- "noLCS_preaph_Tbili"
names(dfc2)[names(dfc2) == 'noLCS_preaph_tnfa'] <- "noLCS_preaph_TNFa"

# Median colorscale values
df_hm <- dfc1  %>% select(c(record_id, preaph_cluster, starts_with("preaph_"))) %>%
  select(!record_id) %>% pivot_longer(!preaph_cluster, names_to = "labs", values_to = "result_value") %>%
  filter(!is.na(preaph_cluster)) %>%
  group_by(preaph_cluster, labs) %>%
  summarize(mean_result=median(result_value, na.rm=TRUE)) %>% ungroup() %>%
  mutate(mean_result=ifelse(is.nan(mean_result), NA, mean_result)) %>%
  mutate(labs=gsub("preaph_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "preaph_cluster", values_from = "mean_result") %>%
  column_to_rownames("labs") %>% as.matrix()
df_hm <- df_hm[order(row.names(df_hm)), ]

# Median values with IQR
df_words <- dfc2  %>% select(c(record_id, preaph_cluster, starts_with("noLCS_preaph_"))) %>%
  select(!record_id) %>% pivot_longer(!preaph_cluster, names_to = "labs", values_to = "result_value") %>%
  filter(!is.na(preaph_cluster)) %>%
  group_by(preaph_cluster, labs) %>%
  summarize(mean_result=median(result_value, na.rm=TRUE),
            low_result=quantile(result_value, probs=0.25, na.rm=TRUE),
            high_result=quantile(result_value, probs=0.75, na.rm=TRUE)) %>% ungroup() %>%
  mutate(mean_result=ifelse(is.nan(mean_result), NA, mean_result)) %>%
  mutate(value=paste0(round(mean_result,2), " (", round(low_result,2), " - ", round(high_result,2), ")")) %>%
  select(!c(mean_result, low_result, high_result)) %>%
  mutate(labs=gsub("noLCS_preaph_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "preaph_cluster", values_from = "value") %>%
  column_to_rownames("labs") %>% as.matrix()
df_words <- df_words[order(row.names(df_words)), ]

df_full_hm <- dfc1 %>% select(c(record_id, preaph_cluster, starts_with("preaph_"))) %>%
  select(!preaph_cluster) %>% column_to_rownames("record_id") %>% t() %>% as.data.frame() %>%
  rownames_to_column("labs") %>%  mutate(labs=gsub("preaph_(.+)", "\\1",labs)) %>% column_to_rownames("labs")
df_full_hm <- df_full_hm[order(row.names(df_full_hm)), ]

dfb_infprob <- reformat_table_tp(dfb1, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="preaph_") %>%
  filter(cluster=="Inflammatory")

dfc3 <- dfb_infprob %>% rename(preaph_cluster=cluster)
names(dfc3)[names(dfc3) == 'preaph_albumin'] <- "preaph_Albumin"
names(dfc3)[names(dfc3) == 'preaph_wbc'] <- "d0_WBC"
names(dfc3)[names(dfc3) == 'preaph_alk'] <- "preaph_ALP"
names(dfc3)[names(dfc3) == 'preaph_ast'] <- "preaph_AST"
names(dfc3)[names(dfc3) == 'preaph_crp'] <- "preaph_CRP"
names(dfc3)[names(dfc3) == 'preaph_ddimer'] <- "preaph_D-dimer"
names(dfc3)[names(dfc3) == 'preaph_ferritin'] <- "preaph_Ferritin"
names(dfc3)[names(dfc3) == 'preaph_hb'] <- "preaph_Hgb"
names(dfc3)[names(dfc3) == 'preaph_il10'] <- "preaph_IL-10"
names(dfc3)[names(dfc3) == 'preaph_il6'] <- "preaph_IL-6"
names(dfc3)[names(dfc3) == 'preaph_ldh'] <- "preaph_LDH"
names(dfc3)[names(dfc3) == 'preaph_plt'] <- "preaph_Plt"
names(dfc3)[names(dfc3) == 'preaph_tbr'] <- "preaph_Tbili"
names(dfc3)[names(dfc3) == 'preaph_tnfa'] <- "preaph_TNFa"

df_full_hm_meta1 <- dfc1 %>% select(c(record_id, preaph_cluster, tau)) %>% column_to_rownames("record_id")%>%
  mutate(preaph_cluster=ifelse(preaph_cluster=="Inflammatory", "#FF7F0EFF", "#1F77B4FF"))

# Generate colorscale for the heatmap as for all heatmaps in the paper
library(circlize)
cp1 <- colorRamp2(c(-10, -1.5, 0, 1.5, 10), c("blue4", "blue3", "#ffffff", "red3", "darkred"))
ro <- as.integer(c(4, 9, 10, 5, 3, 6, 13, 8, 2, 12, 14, 11, 7, 1))

df_rnames <- data.frame(rn=rownames(df_hm))
rownames(df_rnames) <- rownames(df_hm)
hainf = rowAnnotation(foo = anno_text(df_rnames$rn, location=0.6, just="center", gp = gpar(fontsize = 12)))
haninf = rowAnnotation(foo = anno_text(df_rnames$rn, location=0.4, just="center", gp = gpar(fontsize = 12)))
colann_df <- data.frame(a="Inflammatory", cl="#FF7F0EFF")
medhm_inf <- Heatmap(as.matrix(df_hm[,2]), width=unit(4.4, "cm"), col=cp1,
                     show_column_dend = FALSE, show_row_dend=FALSE,
                     show_row_names = FALSE,
                     row_order=ro,
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.text(sprintf("%s", df_words[i, 2]), x, y, gp = gpar(fontsize = 10.5, fontface="bold"))
                     }, show_heatmap_legend = FALSE,
                     right_annotation = hainf,
                     bottom_annotation = HeatmapAnnotation(foo=anno_text("Raw Value\nMedian (IQR)", just = "center", rot=0, location=0.5,
                                                                         height = max_text_height("Non-Inflammatory")*2.5,
                                                                         width = max_text_width("Non-Inflammatory")*2,
                                                                         gp = gpar( fontsize=12, fontface="bold"))),
                     top_annotation=HeatmapAnnotation(foo=anno_text("Inflammatory", just = "center", rot=0, location=0.5,
                                                                    height = max_text_height("Inflammatory")*2.5,
                                                                    gp = gpar(fill ="#FF7F0EFF", col = "white", border = "white", fontsize=12, fontface="bold")))
)

medhm_ninf <- Heatmap(as.matrix(df_hm[,1]), width=unit(4.4, "cm"), col=cp1,
                      column_labels="",
                      column_names_rot = 0,
                      show_column_dend = FALSE, show_row_dend=FALSE,
                      row_order=ro,
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.text(sprintf("%s", df_words[i, 1]), x, y, gp = gpar(fontsize = 10.5, fontface="bold"))
                      }, show_heatmap_legend = FALSE,
                      show_row_names = FALSE,
                      left_annotation = haninf,
                      bottom_annotation = HeatmapAnnotation(foo=anno_text("Raw Value\nMedian (IQR)", just = "center", rot=0, location=0.5,
                                                                          height = max_text_height("Non-Inflammatory")*2.5,
                                                                          width = max_text_width("Non-Inflammatory")*2,
                                                                          gp = gpar( fontsize=12, fontface="bold"))),
                      top_annotation=HeatmapAnnotation(foo=anno_text("Non-\nInflammatory", just = "center", rot=0, location=0.5,
                                                                     height = max_text_height("Non-Inflammatory")*2.5,
                                                                     width = max_text_width("Non-Inflammatory")*2,
                                                                     gp = gpar(fill ="#1F77B4FF", col = "white", border = "white", fontsize=12, fontface="bold"))))

column_ha = HeatmapAnnotation(`Cluster\nAssig.\nProb.` = anno_barplot(df_full_hm_meta1$tau,
                                                                      axis_param=list(gp=gpar(fontsize = 10)),
                                                                      annotation_name_gp= gpar(fontsize = 14),
                                                                      bar_width = 1, height = unit(1.5, "cm"), gp = gpar(fill = df_full_hm_meta1$preaph_cluster, col=df_full_hm_meta1$preaph_cluster)))
dfc3_ord <- data.frame(record_id=colnames(df_full_hm)) %>% left_join(dfc3 %>% select(record_id, inf_tau=tau), by="record_id")
a <- dfc3_ord[rev(order(dfc3_ord$inf_tau)),] %>% rownames_to_column("ordr") %>% select(record_id, ordr, inf_tau) %>% select(!ordr) %>%
  rownames_to_column("ordr")
b <- data.frame(record_id=colnames(df_full_hm)) %>% left_join(a, by="record_id") %>% rownames_to_column("trordr") %>%
  select(ordr, trordr, inf_tau) %>% arrange(-inf_tau)

corehm <- Heatmap(as.matrix(df_full_hm), col=cp1, column_gap=unit(1, "cm"),
                  show_column_names = FALSE,
                  row_order=ro,
                  column_order=as.numeric(b$trordr),
                  cluster_columns = FALSE,
                  row_names_side = "left",
                  show_column_dend = FALSE, show_row_dend=FALSE,
                  top_annotation=column_ha, cluster_rows = FALSE,
                  #show_heatmap_legend = FALSE,
                  heatmap_legend_param = list(direction = "horizontal",
                                              legend_width = unit(15, "cm"),
                                              title = "Normalized and Scaled Value",
                                              title_gp=gpar(fontsize=13, fontface="bold"),
                                              title_position = "topcenter"),
                  column_title="All Patients",
                  column_title_side="top",
                  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

# Figure 4b
ht_list <- medhm_inf + corehm +medhm_ninf
ht_opt$TITLE_PADDING = unit(c(0.2, 0), "cm")
ht_opt$COLUMN_ANNO_PADDING = unit(0.3, "cm")
fig4b <- grid.grabExpr(draw(ht_list, ht_gap = unit(0.3, "cm"),heatmap_legend_side = "bottom"))
fig4b_fin <- plot_grid(fig4b, labels=c("B"), scale=0.93, label_size=30)


####### Inferences for cluster transition outcomes ####

# Dummy variables to allow for rowbinding later on.
crp_cox_cindex_os <- as.data.frame(t(c(crp_coxph_cindex=NA, crp_coxph_cindex_contrast=NA, crp_pvalue_cindex_contrast=NA)))
crp_cox_cindex_pfs<- as.data.frame(t(c(crp_coxph_cindex=NA, crp_coxph_cindex_contrast=NA, crp_pvalue_cindex_contrast=NA)))
crpferr_cox_cindex_os<- as.data.frame(t(c(crpferr_coxph_cindex=NA, crpferr_coxph_cindex_contrast=NA, crpferr_pvalue_cindex_contrast=NA)))
crpferr_cox_cindex_pfs<- as.data.frame(t(c(crpferr_coxph_cindex=NA, crpferr_coxph_cindex_contrast=NA, crpferr_pvalue_cindex_contrast=NA)))

# Identify how patients are transitioning between clusters from preapheresis to preinfusion.
fig6c_df_inf <- prealluv %>%
  filter(!is.na(preaph_cluster)) %>%
  filter(!is.na(d0_cluster)) %>%
  mutate(
    h1_cluster=ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Inflammatory", "Inf -> Inf",
                      ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Non-Inflammatory","Inf -> Non-Inf",
                             ifelse(preaph_cluster=="Non-Inflammatory" & d0_cluster=="Inflammatory","Non-Inf -> Inf",
                                    "Non-Inf -> Non-Inf")))) %>%
  filter(preaph_cluster=="Inflammatory") %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Inf -> Inf","Inf -> Non-Inf"))) %>%
  mutate(h1_cluster=h1_cluster) %>% mutate(Cluster=h1_cluster) %>% rownames_to_column("record_id")

# Logistic regression and Cox PH regression for disease response and survival when transitioning from Inf -> Non. Inf
fig6c_inf_cr <- tidy(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)), data=fig6c_df_inf, family=binomial())) %>%
  filter(term=="ClusterInf -> Non-Inf") %>% select(estimate, p.value) %>% mutate(outcome="No CR") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)),data=fig6c_df_inf,family=binomial()))["ClusterInf -> Non-Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)),data=fig6c_df_inf,family=binomial()))["ClusterInf -> Non-Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))

fig6c_inf_pfs <-  tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data=fig6c_df_inf)) %>%
  filter(term=="ClusterInf -> Non-Inf") %>% select(estimate, p.value) %>% mutate(outcome="PFS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data = fig6c_df_inf))["ClusterInf -> Non-Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data = fig6c_df_inf))["ClusterInf -> Non-Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>% mutate(mm_cindex=NA) %>% bind_cols(crp_cox_cindex_pfs, crpferr_cox_cindex_pfs) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))

fig6c_inf_os <-  tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data=fig6c_df_inf)) %>%
  filter(term=="ClusterInf -> Non-Inf") %>% select(estimate, p.value) %>% mutate(outcome="OS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data = fig6c_df_inf))["ClusterInf -> Non-Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data = fig6c_df_inf))["ClusterInf -> Non-Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>% mutate(mm_cindex=NA) %>% bind_cols(crp_cox_cindex_os, crpferr_cox_cindex_os) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))

# Logistic regression and Cox PH regression for disease response and survival when transitioning from Non.Inf -> Inf
fig6c_df_ninf <- prealluv %>%
  filter(!is.na(preaph_cluster)) %>%
  filter(!is.na(d0_cluster)) %>%
  mutate(
    h1_cluster=ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Inflammatory", "Inf -> Inf",
                      ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Non-Inflammatory","Inf -> Non-Inf",
                             ifelse(preaph_cluster=="Non-Inflammatory" & d0_cluster=="Inflammatory","Non-Inf -> Inf",
                                    "Non-Inf -> Non-Inf")))) %>%
  filter(preaph_cluster=="Non-Inflammatory") %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Non-Inf -> Non-Inf", "Non-Inf -> Inf"))) %>%
  mutate(h1_cluster=h1_cluster) %>% mutate(Cluster=h1_cluster) %>% rownames_to_column("record_id")

fig6c_ninf_cr <- tidy(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)), data=fig6c_df_ninf, family=binomial())) %>%
  filter(term=="ClusterNon-Inf -> Inf") %>% select(estimate, p.value) %>% mutate(outcome="No CR") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)),data=fig6c_df_ninf,family=binomial()))["ClusterNon-Inf -> Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(glm(as.formula(paste0("everCR_100", "~Cluster", trans_covar)),data=fig6c_df_ninf,family=binomial()))["ClusterNon-Inf -> Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))

fig6c_ninf_pfs <-  tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data=fig6c_df_ninf)) %>%
  filter(term=="ClusterNon-Inf -> Inf") %>% select(estimate, p.value) %>% mutate(outcome="PFS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data = fig6c_df_ninf))["ClusterNon-Inf -> Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", trans_covar)), data = fig6c_df_ninf))["ClusterNon-Inf -> Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>% mutate(mm_cindex=NA) %>% bind_cols(crp_cox_cindex_pfs, crpferr_cox_cindex_pfs) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))

fig6c_ninf_os <-  tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data=fig6c_df_ninf)) %>%
  filter(term=="ClusterNon-Inf -> Inf") %>% select(estimate, p.value) %>% mutate(outcome="OS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
  mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data = fig6c_df_ninf))["ClusterNon-Inf -> Inf", "2.5 %"])) %>%
  mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", trans_covar)), data = fig6c_df_ninf))["ClusterNon-Inf -> Inf", "97.5 %"])) %>%
  mutate(covariates=trans_covar) %>% mutate(mm_cindex=NA) %>% bind_cols(crp_cox_cindex_os, crpferr_cox_cindex_os) %>%
  mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2), ) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2))


# Condense all the metrics into one data frame.
fig6cde_metrics <- bind_rows(
  fig6c_inf_cr %>% mutate(analysis="inf_ninf_trans"),
  fig6c_inf_pfs %>% mutate(analysis="inf_ninf_trans"),
  fig6c_inf_os %>% mutate(analysis="inf_ninf_trans"),
  fig6c_ninf_cr %>% mutate(analysis="ninf_inf_trans"),
  fig6c_ninf_pfs %>% mutate(analysis="ninf_inf_trans"),
  fig6c_ninf_os %>% mutate(analysis="ninf_inf_trans")
)

fig6cde_df <- prealluv %>%
  filter(!is.na(preaph_cluster)) %>%
  filter(!is.na(d0_cluster)) %>%
  mutate(
    h1_cluster=ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Inflammatory", "Inf -> Inf",
                      ifelse(preaph_cluster=="Inflammatory" & d0_cluster=="Non-Inflammatory","Inf -> Non-Inf",
                             ifelse(preaph_cluster=="Non-Inflammatory" & d0_cluster=="Inflammatory","Non-Inf -> Inf",
                                    "Non-Inf -> Non-Inf")))) %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Inf -> Inf","Inf -> Non-Inf", "Non-Inf -> Inf", "Non-Inf -> Non-Inf"))) %>%
  mutate(h1_cluster=h1_cluster)


metric_inf <-  fig6cde_metrics %>% filter(outcome=="No CR" & analysis=="inf_ninf_trans")
metric_ninf <-  fig6cde_metrics %>% filter(outcome=="No CR" & analysis=="ninf_inf_trans")

# Plot Figure 6e
fig6e_pre <- fig6cde_df %>% group_by(h1_cluster, everCR_100) %>% reframe(n=n()) %>%
  group_by(h1_cluster) %>% mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>%
  filter(everCR_100=="CR") %>%
  mutate(h1_cluster=factor(h1_cluster, levels=c("Inf -> Inf", "Inf -> Non-Inf", "Non-Inf -> Non-Inf", "Non-Inf -> Inf"))) %>%
  mutate(patfil=factor(h1_cluster, levels=c("Non-Inf -> Inf", "Non-Inf -> Non-Inf", "Inf -> Non-Inf", "Inf -> Inf"))) %>%
  mutate(praph_clus=ifelse(h1_cluster=="Inf -> Inf" | h1_cluster=="Inf -> Non-Inf", "Infl. -> Infl.             Infl. -> Non-Infl.", "Non-Infl. -> Non-Infl.      Non-Infl. -> Infl.")) %>%
  ggplot(aes(x=praph_clus, y=prop,
             fill=h1_cluster, label=paste0(round(prop,2),"\n",n,"/",tot)
  )) +
  geom_bar_pattern(
    aes(fill=h1_cluster, pattern=h1_cluster, pattern_fill=patfil),
    stat="identity", color="black",
    position=position_dodge(),
  ) +
  scale_pattern_fill_manual(values = c("darkorange", "darkorange", "#1261af","#1261af")) +
  scale_pattern_manual(values=c("none", 'stripe', "none", 'stripe')) +
  scale_fill_manual(values = c("darkorange", "darkorange", "#1261af","#1261af")) +
  geom_label(aes(y=prop/2, group=h1_cluster), color="black", fill="white",
             fontface="bold", position = position_dodge(width = 1), size=8) +
  guides(fill = FALSE, pattern=FALSE, fill=FALSE, pattern_fill=FALSE) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(angle = 0, size=13),
    axis.title.x = element_text(angle = 0, size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
  ) +
  labs(x = "",
       y = "Proportion CR",
       title=paste0(""),
  ) + ylim(0, 1) +
  annotate(geom = "text", x=1, y=0.8, label = paste0(
    "Adj. HR (95% CI): ", metric_inf$expEstimate, " (", metric_inf$low_ci, " - ", metric_inf$high_ci, ")\n", "p ",  cpval(metric_inf$pvalue)
  ), size=4.5) +
  annotate(geom = "text", x=2, y=0.8, label = paste0(
    "Adj. HR (95% CI): ", metric_ninf$expEstimate, " (", metric_ninf$low_ci, " - ", metric_ninf$high_ci, ")\n", "p ",  cpval(metric_ninf$pvalue)
  ), size=4.5)
fig6e_fin <- plot_grid(fig6e_pre,labels=c("E"), scale=0.95, label_size = 30)


# Plot survival curves for both sets of cluster transitions for PFS and OS (Figures 6c-d)
fig6c_temp <- survfit(Surv(tt_pfs_m, ev_pfs) ~ h1_cluster, data = fig6cde_df)
names(fig6c_temp$strata) <- c("Infl. -> Infl.","Infl. -> Non-Infl.", "Non-Infl. -> Infl.", "Non-Infl. -> Non-Infl.")

metric_infninf <- fig6cde_metrics %>% filter(outcome=="PFS" & analysis=="inf_ninf_trans")
metric_ninfinf <- fig6cde_metrics %>% filter(outcome=="PFS" & analysis=="ninf_inf_trans")
fig6c_pre <- fig6c_temp %>% ggsurvfit(linetype_aes = TRUE, linewidth=0.75)
fig6c_pre <- fig6c_pre +
  coord_cartesian(xlim = c(0, 25)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent,
    expand = c(0.01, 0)
  ) +
  scale_linetype_manual(values=c("solid", "longdash", "longdash", "solid"))+
  scale_color_manual(values = c("#FF7F0EFF", "darkorange3", "royalblue4","#1F77B4FF")) +
  theme_minimal() +
  add_risktable(risktable_height=0.2, size=5,
                risktable_stats = "n.risk",
                theme = list( theme_risktable_default(axis.text.y.size = 13, plot.title.size = 13), theme(plot.title = element_text(face = "bold")) )
  )+
  add_risktable_strata_symbol(symbol = "\U25CF", size = 25)+
  add_censor_mark() +
  labs(
    title = "",
    y = "PFS Probability",
    x="Months from CAR-T Infusion"
  ) +
  guides(color="legend")+
  theme(
    strip.background  = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    axis.text.x = element_text(angle = 0, size=15),
    axis.text.y = element_text(angle = 0, size=15),
    axis.title.x = element_text(angle = 0, size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"),
    title = element_text(angle = 0, size=18, face="bold"),
    legend.text = element_text(size=16),
    legend.position = "bottom"
  ) +
  annotate(geom = "text", x=15, y=0.9, label = paste0(
    "Infl. -> Non-Infl. Adj. OR (95% CI): ", metric_infninf$expEstimate, " (", metric_infninf$low_ci, " - ", metric_infninf$high_ci, ")\n", "p ",  cpval(metric_infninf$pvalue)
  ), size=6) +
  annotate(geom = "text", x=15, y=0.7, label = paste0(
    "Infl. -> Infl. Adj. OR (95% CI): ", metric_ninfinf$expEstimate, " (", metric_ninfinf$low_ci, " - ", metric_ninfinf$high_ci, ")\n", "p ",  cpval(metric_ninfinf$pvalue)
  ), size=6)

legend <- get_legend(
  fig6c_pre + theme(legend.box.margin = margin(0, 0, 0, 0))
)

fig6d_temp <- survfit(Surv(tt_os_m, ev_os) ~ h1_cluster, data = fig6cde_df)
names(fig6d_temp$strata) <- c("Infl. -> Infl.","Infl. -> Non-Infl.", "Non-Infl. -> Infl.", "Non-Infl. -> Non-Infl.")

metric_infninf <- fig6cde_metrics %>% filter(outcome=="OS" & analysis=="inf_ninf_trans")
metric_ninfinf <- fig6cde_metrics %>% filter(outcome=="OS" & analysis=="ninf_inf_trans")
fig6d_pre <- fig6d_temp %>% ggsurvfit(linetype_aes = TRUE, linewidth=0.75) +
  coord_cartesian(xlim = c(0, 25)) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent,
    expand = c(0.01, 0)
  ) +
  scale_linetype_manual(values=c("solid", "longdash", "longdash", "solid"))+
  scale_color_manual(values = c("#FF7F0EFF", "darkorange3", "royalblue4","#1F77B4FF")) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  add_risktable(risktable_height=0.2, size=5,
                risktable_stats = "n.risk",
                theme = list( theme_risktable_default(axis.text.y.size = 13, plot.title.size = 13), theme(plot.title = element_text(face = "bold")) )
  )+
  add_risktable_strata_symbol(symbol = "\U25CF", size = 25)+
  add_censor_mark() +
  labs(
    title = "",
    y = "OS Probability",
    x="Months from CAR-T Infusion"
  ) +
  guides(color="legend", linetype="legend")+
  theme(
    axis.text.x = element_text(angle = 0, size=15),
    axis.text.y = element_text(angle = 0, size=15),
    axis.title.x = element_text(angle = 0, size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"),
    legend.text = element_text(size=16),
    legend.position="bottom"
  )+
  annotate(geom = "text", x=15, y=0.9, label = paste0(
    "Infl. -> Non-Infl. Adj. OR (95% CI): ", metric_infninf$expEstimate, " (", metric_infninf$low_ci, " - ", metric_infninf$high_ci, ")\n", "p ",  cpval(metric_infninf$pvalue)
  ), size=6) +
  annotate(geom = "text", x=5, y=0.2, label = paste0(
    "Non-Infl. -> Infl. Adj. OR (95% CI): ", metric_ninfinf$expEstimate, " (", metric_ninfinf$low_ci, " - ", metric_ninfinf$high_ci, ")\n", "p ",  cpval(metric_ninfinf$pvalue)
  ), size=6)
