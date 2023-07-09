rm(list=ls())
load("output/model.RData")
load("output/scaledata.RData")
source("scripts/functions_constants.R")

library(colorspace)
library(RColorBrewer)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(cowplot)
library(ComplexHeatmap)

dfb1 <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  filter(cohort=="MSK Development") %>% # Change this filter to generate a heatmap for each cohort
  mutate(record_id=record_id)
dfb <- reformat_table_tp(dfb1, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup()

dfc1 <- dfb %>% rename(d0_cluster=cluster)
dfc2 <- dfb %>% rename(d0_cluster=cluster)

names(dfc1)[names(dfc1) == 'd0_albumin'] <- "d0_Albumin"
names(dfc1)[names(dfc1) == 'd0_wbc'] <- "d0_WBC"
names(dfc1)[names(dfc1) == 'd0_alk'] <- "d0_ALP"
names(dfc1)[names(dfc1) == 'd0_ast'] <- "d0_AST"
names(dfc1)[names(dfc1) == 'd0_crp'] <- "d0_CRP"
names(dfc1)[names(dfc1) == 'd0_ddimer'] <- "d0_Ddimer"
names(dfc1)[names(dfc1) == 'd0_ferritin'] <- "d0_Ferritin"
names(dfc1)[names(dfc1) == 'd0_hb'] <- "d0_Hgb"
names(dfc1)[names(dfc1) == 'd0_il10'] <- "d0_IL-10"
names(dfc1)[names(dfc1) == 'd0_il6'] <- "d0_IL-6"
names(dfc1)[names(dfc1) == 'd0_ldh'] <- "d0_LDH"
names(dfc1)[names(dfc1) == 'd0_plt'] <- "d0_Plt"
names(dfc1)[names(dfc1) == 'd0_tbr'] <- "d0_Tbili"
names(dfc1)[names(dfc1) == 'd0_tnfa'] <- "d0_TNFa"

names(dfc2)[names(dfc2) == 'noLCS_d0_albumin'] <- "noLCS_d0_Albumin"
names(dfc2)[names(dfc2) == 'noLCS_d0_wbc'] <- "noLCS_d0_WBC"
names(dfc2)[names(dfc2) == 'noLCS_d0_alk'] <- "noLCS_d0_ALP"
names(dfc2)[names(dfc2) == 'noLCS_d0_ast'] <- "noLCS_d0_AST"
names(dfc2)[names(dfc2) == 'noLCS_d0_crp'] <- "noLCS_d0_CRP"
names(dfc2)[names(dfc2) == 'noLCS_d0_ddimer'] <- "noLCS_d0_Ddimer"
names(dfc2)[names(dfc2) == 'noLCS_d0_ferritin'] <- "noLCS_d0_Ferritin"
names(dfc2)[names(dfc2) == 'noLCS_d0_hb'] <- "noLCS_d0_Hgb"
names(dfc2)[names(dfc2) == 'noLCS_d0_il10'] <- "noLCS_d0_IL-10"
names(dfc2)[names(dfc2) == 'noLCS_d0_il6'] <- "noLCS_d0_IL-6"
names(dfc2)[names(dfc2) == 'noLCS_d0_ldh'] <- "noLCS_d0_LDH"
names(dfc2)[names(dfc2) == 'noLCS_d0_plt'] <- "noLCS_d0_Plts"
names(dfc2)[names(dfc2) == 'noLCS_d0_tbr'] <- "noLCS_d0_Tbili"
names(dfc2)[names(dfc2) == 'noLCS_d0_tnfa'] <- "noLCS_d0_TNFa"

df_hm <- dfc1  %>% select(c(record_id, d0_cluster, starts_with("d0_"))) %>%
  select(!record_id) %>% pivot_longer(!d0_cluster, names_to = "labs", values_to = "result_value") %>%
  filter(!is.na(d0_cluster)) %>%
  group_by(d0_cluster, labs) %>%
  summarize(mean_result=median(result_value, na.rm=TRUE)) %>% ungroup() %>%
  mutate(mean_result=ifelse(is.nan(mean_result), NA, mean_result)) %>%
  mutate(labs=gsub("d0_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "d0_cluster", values_from = "mean_result") %>%
  column_to_rownames("labs") %>% as.matrix()
df_hm <- df_hm[order(row.names(df_hm)), ]

df_words <- dfc2  %>% select(c(record_id, d0_cluster, starts_with("noLCS_d0_"))) %>%
  select(!record_id) %>% pivot_longer(!d0_cluster, names_to = "labs", values_to = "result_value") %>%
  filter(!is.na(d0_cluster)) %>%
  group_by(d0_cluster, labs) %>%
  summarize(mean_result=median(result_value, na.rm=TRUE),
            low_result=quantile(result_value, probs=0.25, na.rm=TRUE),
            high_result=quantile(result_value, probs=0.75, na.rm=TRUE)) %>% ungroup() %>%
  mutate(mean_result=ifelse(is.nan(mean_result), NA, mean_result)) %>%
  mutate(value=paste0(round(mean_result,2), " (", round(low_result,2), " - ", round(high_result,2), ")")) %>%
  select(!c(mean_result, low_result, high_result)) %>%
  mutate(labs=gsub("noLCS_d0_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "d0_cluster", values_from = "value") %>%
  column_to_rownames("labs") %>% as.matrix()
df_words <- df_words[order(row.names(df_words)), ]

df_full_hm <- dfc1 %>% select(c(record_id, d0_cluster, starts_with("d0_"))) %>%
  select(!d0_cluster) %>% column_to_rownames("record_id") %>% t() %>% as.data.frame() %>%
  rownames_to_column("labs") %>%  mutate(labs=gsub("d0_(.+)", "\\1",labs)) %>% column_to_rownames("labs")
df_full_hm <- df_full_hm[order(row.names(df_full_hm)), ]

dfb_infprob <- reformat_table_tp(dfb1, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  filter(cluster=="Inflammatory")

dfc3 <- dfb_infprob %>% rename(d0_cluster=cluster)
names(dfc3)[names(dfc3) == 'd0_albumin'] <- "d0_Albumin"
names(dfc3)[names(dfc3) == 'd0_wbc'] <- "d0_WBC"
names(dfc3)[names(dfc3) == 'd0_alk'] <- "d0_ALP"
names(dfc3)[names(dfc3) == 'd0_ast'] <- "d0_AST"
names(dfc3)[names(dfc3) == 'd0_crp'] <- "d0_CRP"
names(dfc3)[names(dfc3) == 'd0_ddimer'] <- "d0_Ddimer"
names(dfc3)[names(dfc3) == 'd0_ferritin'] <- "d0_Ferritin"
names(dfc3)[names(dfc3) == 'd0_hb'] <- "d0_Hgb"
names(dfc3)[names(dfc3) == 'd0_il10'] <- "d0_IL-10"
names(dfc3)[names(dfc3) == 'd0_il6'] <- "d0_IL-6"
names(dfc3)[names(dfc3) == 'd0_ldh'] <- "d0_LDH"
names(dfc3)[names(dfc3) == 'd0_plt'] <- "d0_Plt"
names(dfc3)[names(dfc3) == 'd0_tbr'] <- "d0_Tbili"
names(dfc3)[names(dfc3) == 'd0_tnfa'] <- "d0_TNFa"

df_full_hm_meta1 <- dfc1 %>% select(c(record_id, d0_cluster, tau)) %>% column_to_rownames("record_id")%>%
  mutate(d0_cluster=ifelse(d0_cluster=="Inflammatory", "#FF7F0EFF", "#1F77B4FF"))

library(circlize)
cp1 <- colorRamp2(c(-10, -1.5, 0, 1.5, 10), c("blue4", "blue3", "#ffffff", "red3", "darkred"))
ro <- as.integer(c(4, 9, 10, 5, 3, 6, 13, 8, 2, 12, 14, 11, 7, 1))

df_rnames <- data.frame(rn=rownames(df_hm))
rownames(df_rnames) <- rownames(df_hm)
hainf = rowAnnotation(foo = anno_text(df_rnames$rn, location=0.6, just="center", gp = gpar(fontsize = 12)))
haninf = rowAnnotation(foo = anno_text(df_rnames$rn, location=0.4, just="center", gp = gpar(fontsize = 12)))
colann_df <- data.frame(a="Inflammatory", cl="#FF7F0EFF")
medhm_inf <- Heatmap(as.matrix(df_hm[,2]), width=unit(4.1, "cm"), col=cp1,
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

medhm_ninf <- Heatmap(as.matrix(df_hm[,1]), width=unit(4.1, "cm"), col=cp1,
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
                                                                      bar_width = 1, height = unit(1.5, "cm"), gp = gpar(fill = df_full_hm_meta1$d0_cluster)))
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
                  column_title="I. MSK LBCL Derivation Cohort Patients", # Change this title to reflect the cohort
                  column_title_side="top",
                  column_title_gp = gpar(fontsize = 14, fontface = "bold")
)

ht_list <- medhm_inf + corehm +medhm_ninf
ht_opt$TITLE_PADDING = unit(c(0.2, 0), "cm")
ht_opt$COLUMN_ANNO_PADDING = unit(0.3, "cm")
efig3a_d <- grid.grabExpr(draw(ht_list, ht_gap = unit(0.3, "cm"),heatmap_legend_side = "bottom"))
efig3a_d_fin <- plot_grid(efig3a_d, labels=c("C"), scale=0.93, label_size=30)
