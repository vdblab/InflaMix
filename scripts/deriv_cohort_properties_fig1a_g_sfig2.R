rm(list=ls())
load(bstfun::here_data("output/model.RData"))
load(bstfun::here_data("output/scaledata.RData"))
source("scripts/analysis/functions_constants.R")

library(corrplot)
library(umap)
library(ggsci)
library(pheatmap)
library(colorspace)
library(ggpubr)
library(mclust)
library(mvtnorm)
library(tidyverse)
library(tidymodels)
library(ranger)
library(vip)
library(RColorBrewer)
library(gtable)
library(gridExtra)
library(gridGraphics)
library(cowplot)

############ Figure 1a Correlation Plot ##########

# Select dataset
dfa <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  filter(cohort=="MSK Development") %>%
  select(record_id, starts_with("d0_")) %>%
  select(record_id, contains(cluster_labs)) %>%
  drop_na() %>% column_to_rownames("record_id") %>% select(starts_with("d0_"))

# Uncomment these lines to instead evaluate all 16 labs including creatinine and alt to generate supplementary figure 2
# dfa <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
#   filter(analysis_type==1) %>%
#   filter(cohort=="MSK Development") %>%
#   select(record_id, starts_with("d0_")) %>%
#   select(record_id, contains(c(cluster_labs, "creatinine", "alt"))) %>%
#   drop_na() %>% column_to_rownames("record_id") %>% select(starts_with("d0_"))

colnames(dfa) <- c("Albumin",
                      "ALP",
                      "AST",
                      "Ferritin",
                      "Hgb",
                      "LDH",
                      "Plt",
                      "Tbili",
                      "IL-10",
                      "IL-6",
                      "TNFa",
                      "CRP",
                      "Ddimer",
                   "WBC"
)

# Generate correlation matrix
cormat <- as.data.frame(cor(dfa))
rownames(cormat) <- colnames(cormat)

# Evaluate p values
testRes = cor.mtest(dfa, conf.level = 0.95)
testRes1 <- testRes

# Adjust p values via FDR - don't double count values and ignore diagonals.
testResadjust = matrix(p.adjust(testRes[[1]][lower.tri(testRes[[1]], diag=FALSE)], method = "fdr"))
testRes1[[1]][lower.tri(testRes1[[1]], diag=FALSE)] <- testResadjust
diag(testRes1[[1]]) <- rep(1, dim(testRes1[[1]])[1])
testResadjust <- as.matrix(testRes1[[1]])

colnames(testResadjust) <- colnames(cormat)
rownames(testResadjust) <- colnames(testResadjust)

# Function to plot correlation plot
corfun <- function(x){
  corrplot(as.matrix(x),
           method="square",
           order="hclust",
           type="lower",
           tl.col="black",
           tl.srt = 0.5,
           tl.cex=1.3,
           col = rev(COL2('RdBu', 21)),
           title = "\n",
           p.mat = testResadjust,
           sig.level = c(0.99),
           pch.col="black",
           number.cex=1.2,
           pch.cex = 0.000001,
           #insig = "",
           cl.pos="n",
           diag=FALSE
  )$corrPos -> p1
  # Threshold p values to different symbols by function cpval
  p2 <- as.data.frame(p1) %>% mutate(siglab = ifelse(p.value > 0.05, "",
                                                            ifelse(p.value < 0.05 & p.value > 0.01, "*",
                                                                   ifelse(p.value < 0.01 & p.value > 0.001, "**", "***"))))
  text(p1$x, p1$y, paste0(round(p1$corr, 2)),cex = 1.2, pos=3, offset=0.3)
  text(p1$x, p1$y, p2$siglab, pos=1,cex = 1.4, offset=0.5)
  colorlegend(xlim=c(12.5,13.5), ylim=c(5,12), rev(colorRampPalette(brewer.pal(11, "RdBu"))(11)), c(seq(-1,1,.1)), align="l", vertical=TRUE, addlabels=TRUE)
  recordPlot()
}

# Call plotting function
fig1a_pre <- corfun(cormat)
grid.echo()
fig1a_grob <- grid.grab()

# Plot the colorscale legend
matrix.colors <- getGrob(fig1a_grob, gPath("square"), grep = TRUE)[["gp"]][["fill"]]
fig1a_grob <- fig2a_grob %>%
  editGrob(gPath("square"),
           grep = TRUE,
           gp = gpar(col = NA,fill = NA)) %>%
  editGrob(gPath("symbols-rect-1"),
           grep = TRUE,
           gp = gpar(fill = matrix.colors)) %>%
  editGrob(gPath("background"), grep = TRUE,
           gp = gpar(fill = NA))

dev.off()

# Plot the label legend
labelsig_tbl <- data.frame(Label=c("*", "**", "***"), p = c("< 0.05", "< 0.01", "< 0.001")) %>%
  tableGrob(rows=NULL, cols=c("", "FDR-adjusted\np-value"), theme=ttheme_minimal(base_size=16))
g <- gtable_add_grob(labelsig_tbl,
                     grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                     t = 2, b = nrow(labelsig_tbl), l = 1, r = ncol(labelsig_tbl)) %>%
  gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                 t = 1, l = 1, r = ncol(labelsig_tbl))

# Final figure 1a
fig1a <- ggdraw(fig1a_grob) + draw_grob(g, x=0.4, y=0.6, width=0.3, height=0.4)

############ Figure 1b Heatmap #######
dfb1 <- left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  filter(cohort=="MSK Development") %>%
  mutate(record_id=record_id)

dfb <- reformat_table_tp(dfb1, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  group_by(record_id) %>% dplyr::slice(which.max(tau)) %>% ungroup()

# Figure 1b - median heatmap for the derivation cohort
dfc1 <- dfb %>% rename(d0_cluster=cluster)
dfc2 <- dfb %>% rename(d0_cluster=cluster)

names(dfc1)[names(dfc1) == 'd0_albumin'] <- "d0_Albumin"
names(dfc1)[names(dfc1) == 'd0_wbc'] <- "d0_WBC"
names(dfc1)[names(dfc1) == 'd0_alk'] <- "d0_ALP"
names(dfc1)[names(dfc1) == 'd0_ast'] <- "d0_AST"
names(dfc1)[names(dfc1) == 'd0_crp'] <- "d0_CRP"
names(dfc1)[names(dfc1) == 'd0_ddimer'] <- "d0_D-dimer"
names(dfc1)[names(dfc1) == 'd0_ferritin'] <- "d0_Ferritin"
names(dfc1)[names(dfc1) == 'd0_hb'] <- "d0_Hgb"
names(dfc1)[names(dfc1) == 'd0_il10'] <- "d0_IL-10"
names(dfc1)[names(dfc1) == 'd0_il6'] <- "d0_IL-6"
names(dfc1)[names(dfc1) == 'd0_ldh'] <- "d0_LDH"
names(dfc1)[names(dfc1) == 'd0_plt'] <- "d0_Plts"
names(dfc1)[names(dfc1) == 'd0_tbr'] <- "d0_Tbili"
names(dfc1)[names(dfc1) == 'd0_tnfa'] <- "d0_TNFa"

names(dfc2)[names(dfc2) == 'noLCS_d0_albumin'] <- "noLCS_d0_Albumin"
names(dfc2)[names(dfc2) == 'noLCS_d0_wbc'] <- "noLCS_d0_WBC"
names(dfc2)[names(dfc2) == 'noLCS_d0_alk'] <- "noLCS_d0_ALP"
names(dfc2)[names(dfc2) == 'noLCS_d0_ast'] <- "noLCS_d0_AST"
names(dfc2)[names(dfc2) == 'noLCS_d0_crp'] <- "noLCS_d0_CRP"
names(dfc2)[names(dfc2) == 'noLCS_d0_ddimer'] <- "noLCS_d0_D-dimer"
names(dfc2)[names(dfc2) == 'noLCS_d0_ferritin'] <- "noLCS_d0_Ferritin"
names(dfc2)[names(dfc2) == 'noLCS_d0_hb'] <- "noLCS_d0_Hgb"
names(dfc2)[names(dfc2) == 'noLCS_d0_il10'] <- "noLCS_d0_IL-10"
names(dfc2)[names(dfc2) == 'noLCS_d0_il6'] <- "noLCS_d0_IL-6"
names(dfc2)[names(dfc2) == 'noLCS_d0_ldh'] <- "noLCS_d0_LDH"
names(dfc2)[names(dfc2) == 'noLCS_d0_plt'] <- "noLCS_d0_Plts"
names(dfc2)[names(dfc2) == 'noLCS_d0_tbr'] <- "noLCS_d0_Tbili"
names(dfc2)[names(dfc2) == 'noLCS_d0_tnfa'] <- "noLCS_d0_TNFa"

# Data frame for median values by cluster and IQR
df_words <- dfc2  %>% select(c(record_id, d0_cluster, starts_with("noLCS_d0_"))) %>%
  select(!record_id) %>% pivot_longer(!d0_cluster, names_to = "labs", values_to = "result_value") %>%
  group_by(d0_cluster, labs) %>%
  summarize(mean_result=median(result_value),
            low_result=quantile(result_value, probs=0.25),
            high_result=quantile(result_value, probs=0.75)) %>% ungroup() %>%
  mutate(value=paste0(round(mean_result,2), " (", round(low_result,2), " - ", round(high_result,2), ")")) %>%
  select(!c(mean_result, low_result, high_result)) %>%
  mutate(labs=gsub("noLCS_d0_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "d0_cluster", values_from = "value") %>%
  column_to_rownames("labs") %>% as.matrix()

# Data frame for the heatmap of median values by cluster
df_hm <- dfc1  %>% select(c(record_id, d0_cluster, starts_with("d0_"))) %>%
  select(!record_id) %>% pivot_longer(!d0_cluster, names_to = "labs", values_to = "result_value") %>%
  group_by(d0_cluster, labs) %>%
  summarize(mean_result=median(result_value)) %>% ungroup() %>%
  mutate(labs=gsub("d0_(.+)", "\\1",labs))  %>% pivot_wider(id_cols = labs, names_from = "d0_cluster", values_from = "mean_result") %>%
  column_to_rownames("labs") %>% as.matrix()

# Color scale for heatmap - values closer to 0 are weighted higher to better capture the dynamic range
cp1 <- colorRamp2(c(-10, -1.5, 0, 1.5, 10), c("blue4", "blue3", "#ffffff", "red3", "darkred"))
cp11 <- cp1(seq(-1.2, 1.2, 0.1))

# Function to plot the heatmap
clusheat <- function(x, y) {
  fig1b <- pheatmap::pheatmap(x, fontsize_row=12, fontsize_col = 14, cluster_cols = FALSE,
                              legend=FALSE,
                              labels_col = c("Non-Inflammatory", "Inflammatory"), treeheight_row = 0, fontsize = 14,cluster_rows = TRUE,
                              clustering_distance_rows = "correlation",
                              color=cp11,
                              display_numbers = y,
                              angle_col = 0,
                              number_color="black"
  )
}

# Color scale legend
fig1b_leg <- data.frame(r=factor(round(seq(-1.2, 1.2, 0.1),2), levels=c(round(seq(-1.2, 1.2, 0.1), 2)))) %>% mutate(x="") %>% ggplot(aes(y=r, x=x, fill=r)) + geom_tile() +
  scale_fill_manual(values=cp11) + theme_minimal() + guides(fill=FALSE) +
  labs(x = "",
       y = "",
       title = "\n\n") +
  theme(axis.text.y = element_text(angle = 0, size=12))

# Draw heatmap
fig1b_main <- clusheat(df_hm, df_words)[[4]]

fig1b_layout <- rbind(c(NA, 1,1,1,1,1,1,1, 1,1,1,1,2, 2),
                      c(NA, 1,1,1,1,1,1,1, 1,1,1,1,2, 2),
                      c(NA, 1,1,1,1,1,1,1, 1,1,1,1,2, 2),
                      c(NA, 1,1,1,1,1,1,1, 1,1,1,1,NA, NA))

# Set out the colorscale legend and heatmap onto the plot
fig1b <- ggdraw(grid.arrange(fig1b_main, fig1b_leg, layout_matrix=fig1b_layout)) +
  annotate("text", label ="Scaled\nValue", x =0.95, y =0.95, size=5)+
  theme(plot.margin=margin(0, 1.8, 0, 0))


############ Figure 1c UMAP ##########
umap_data <- dfb %>% select(record_id, starts_with("d0_")) %>% select(record_id, contains(cluster_labs))
umap_meta <- dfb %>% select(record_id, cluster, tau)
set.seed(10)
umap_fit <- umap_data %>% column_to_rownames("record_id") %>% umap(n_components=2)
umap_df <- umap_fit$layout %>% as.data.frame() %>%
  dplyr::rename(UMAP1 = "V1",
                UMAP2 = "V2") %>%
  rownames_to_column("record_id") %>% inner_join(umap_meta, by="record_id") %>%
  dplyr::rename(Cluster=cluster) %>% dplyr::mutate(Tau=tau*100)

#
fig1c <- umap_df %>% ggplot(aes(x = UMAP1, y = UMAP2)) +
  ylim(-5, 5)+
  geom_point(aes(color=Cluster, size=Tau), alpha=0.7)+
  labs(x = "UMAP1",
       y = "UMAP2",
       title = "") +
  theme_bw() +
  scale_fill_manual(values=colorclusters(length(dev_pro)))+
  scale_color_manual(values=colorclusters(length(dev_pro)))+
  theme(
    strip.background  = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.grid.minor.x=element_blank(),
    panel.grid.major.x=element_blank(),
    axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1, size=15),
    axis.text.y = element_text(angle = 0, size=15),
    axis.title.x = element_text(angle = 0, size=18, face="bold"),
    axis.title.y = element_text(size=18, face="bold"),
    title = element_text(angle = 0, size=18, face="bold"),
    legend.text = element_text(angle = 0, size=12),
    legend.title= element_text(size=12),
    legend.position = c(0.15, 0.25)) +
  scale_size_continuous(range  = c(0.1, 6),
                        limits = c(0, 100),
                        breaks = c(0, 50, 100))+
  guides(color=guide_legend(title="", title.position = "top", override.aes = list(size = 6)),
         size=guide_legend(title="Cluster Probability", title.position = "top")
        )

############ Figures 1d-f Tumor burden metrics #####

# Scale pre-lymphodepletion LDH to an assay with an ULN = 250
dfe <- dfb %>% rename(d0_cluster=cluster) %>% select(d0_cluster, bl_suvmax, bl_mtv, noLCS_preld_ldh,uln_preld_ldh) %>%
  mutate(preld_ldh=250*(noLCS_preld_ldh/uln_preld_ldh)) %>% select(!c(noLCS_preld_ldh, uln_preld_ldh))

# Module for calculating p-values by Wilcoxon rank sum for MTV, LDH, and SUVmax (just select preld_ldh or bl_mtv instead of bl_suvmax )
dfe3 <- dfb %>% rename(d0_cluster=cluster) %>% select(record_id, d0_cluster, bl_suvmax, bl_mtv, noLCS_preld_ldh,uln_preld_ldh) %>%
  mutate(preld_ldh=250*(noLCS_preld_ldh/uln_preld_ldh)) %>% select(!c(noLCS_preld_ldh, uln_preld_ldh)) %>% select(record_id, d0_cluster, tmb=bl_suvmax)
dfe3_inf <- dfe3 %>% filter(d0_cluster=="Inflammatory")
dfe3_ninf <- dfe3 %>% filter(d0_cluster=="Non-Inflammatory")
wilcox.test(dfe3_inf$tmb, dfe3_ninf$tmb)
# The three p-values for MTV, SUVmax, and LDH can be adjusted with the following line:
p.adjust(c(0.000173, 7.973e-8, 0.0001388), method="fdr")

dfe1 <- dfe %>% pivot_longer(!d0_cluster, names_to = "feature", values_to = "value")

# Plotting function
dxburd_plot <- function(x, feat, tl, lg){
  pl <- x %>% filter(feature==feat) %>% ggplot(aes(x=d0_cluster, y=log(value+1, 10), fill=d0_cluster), x=1) +
    scale_fill_manual(values=colorclusters(length(dev_pro)))+
    scale_color_manual(values=colorclusters(length(dev_pro)))+
    geom_boxplot(aes(color=d0_cluster, color = after_scale(darken(color, .1, space = "HLS")),
                     fill = after_scale(desaturate(lighten(color, .8), .4))), width = .8,
                 outlier.shape = NA) +
    geom_point(
      aes(
        color = d0_cluster,
        color = after_scale(darken(color, .1, space = "HLS",)),
        fill = d0_cluster,
        fill = after_scale(lighten(color, .1, space = "HLS",))
      ),
      stroke = .4,
      size = 2,
      position=position_jitterdodge(jitter.width = 1,dodge.width = 0.3)
    ) +  theme_bw() +
    labs(x = "",
         y = paste0("Log10-Normalized ", tl),
         title = "",
         fill="Cluster",
         color="Cluster"
    )+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, size=15, face="bold"),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=15, face="bold", vjust=1.2),
      title = element_text(angle = 0, size=18, face="bold"),
      legend.text = element_text(angle = 0, size=14, face="bold"),
      legend.position = lg
    )
  return(pl)
}

# Decoy plot to grab a common colorscale legend
fig2e_suv <- dxburd_plot(dfe1, "bl_suvmax", "SUV Max", "bottom")
legend <- get_legend(fig2e_suv)

# Actual plots
fig1d_suv <-  dxburd_plot(dfe1, "bl_suvmax", "SUVmax", "FALSE") + annotate("text", label ="FDR Adj. p < 0.001", x =1.10, y = 1.7, size=4)
fig1e_mtv <-  dxburd_plot(dfe1, "bl_mtv", "Metabolic Tumor Volume", "FALSE")+ annotate("text", label ="FDR Adj. p < 0.001", x =1.10, y = 3.8, size=4)
fig1f_ldh <-  dxburd_plot(dfe1, "preld_ldh", "Lactate Dehydrogenase", "FALSE")+ annotate("text", label ="FDR Adj. p < 0.001", x =1.10, y = 3.9, size=4)

top <- plot_grid(fig1d_ldh, NULL,
                 fig1e_mtv, NULL,
                 fig1f_suv,
                 labels=c("D", "", "E", "", "F"),
                 rel_widths = c(1, 0.1,
                                1, 0.1,
                                1),
                 nrow=1,
                 scale=1,
                 label_size = 22)

fig1d_f <- plot_grid(top, NULL, legend,
                     labels=c("", "", ""),
                     rel_heights = c(1, 0, 0.05),
                     ncol=1, scale=1
)

############ Figure 1g VIMP plot ######

rf_clus <- dfb %>% column_to_rownames("record_id") %>% select(cluster, starts_with("d0_"))

# Save models, aucs, accuracies, and importance values across different test-train iterations
rf100_models <- tibble()
rf100_aucs <- tibble()
rf100_accuracies <- tibble()
rf100_vimps <- data.frame(
  Variable = c("test"),
  Importance= c(0),
  model=c(0)
)
model_auc <- 0
bestmodelauc <- 0

set.seed(NULL)

# Tidy models framework for tuning random forest models adapted from blogpost by Julia Silge https://juliasilge.com/blog/sf-trees-random-tuning/
# Across 100 iterations of test-train splits of the derivation cohort
for(i in 1:100) {
  old_time = Sys.time()
  set.seed(i)

  # 80-20 split, stratify split by cluster assignment
  clus_split <- initial_split(rf_clus, strata = cluster, prop=0.8)
  clus_train <- training(clus_split)
  clus_test <- testing(clus_split)

  # Predictor variables are labs - predicting cluster assignment
  clus_rec <- recipe(cluster ~ ., data = clus_train)
  clus_juice <- juice(prep(clus_rec))

  tune_spec <- rand_forest(
    mtry = tune(),
    trees = 1000,
    min_n = tune()
  ) %>%
    set_mode("classification") %>%
    set_engine("ranger")

  tune_wf <- workflow() %>%
    add_recipe(clus_rec) %>%
    add_model(tune_spec)

  # 10 repeats of 5 fold cross validation
  set.seed(123+i)
  clus_folds <- vfold_cv(clus_train, strata=cluster, v=5, repeats = 10)
  doParallel::registerDoParallel()

  set.seed(234+i)
  tune_res <- tune_grid(
    tune_wf,
    resamples = clus_folds,
    grid = 20
  )

  # Initial tuning over resamples for cross-folds
  b <- tune_res %>%
    collect_metrics() %>%
    filter(.metric == "roc_auc") %>%
    select(mean, min_n, mtry) %>%
    pivot_longer(min_n:mtry,
                 values_to = "value",
                 names_to = "parameter"
    ) %>%
    ggplot(aes(value, mean, color = parameter)) +
    geom_point(show.legend = FALSE) +
    facet_wrap(~parameter, scales = "free_x") +
    labs(x = NULL, y = "AUC")

  # Dynamic range for hyperparameter space for tuning
  minn_low <- b$data$value[which.max(b$data$mean)]-5
  minn_high <- b$data$value[which.max(b$data$mean)]+5
  mtry_low <-  b$data$value[which.max(b$data$mean)+1]-3
  mtry_high <-  b$data$value[which.max(b$data$mean)+1]+3

  if(minn_low <= 0){minn_low <- 1; minn_high <- minn_high+abs(minn_low)+1}
  if(mtry_low <= 0){mtry_low <- 1; mtry_high <- mtry_high+abs(mtry_low)+1}

  rf_grid <- grid_regular(
    mtry(range = c(mtry_low, mtry_high)),
    min_n(range = c(minn_low, minn_high)),
    levels = 5
  )

  set.seed(345+i)
  regular_res <- tune_grid(
    tune_wf,
    resamples = clus_folds,
    grid = rf_grid
  )

  # Get the AUC from the best model
  best_auc <- select_best(regular_res, "roc_auc")

  final_rf <- finalize_model(
    tune_spec,
    best_auc
  )

  # Extract variable importances from best model
  vip_rf <-
    final_rf %>%
    set_engine("ranger", importance = "permutation") %>%
    fit(cluster ~ .,
        data = clus_juice
    ) %>%
    vip(geom = "point",num_features = 18)

  rf100_vimps <- bind_rows(
    rf100_vimps,
    data.frame(
      Variable = vip_rf$data$Variable,
      Importance= vip_rf$data$Importance,
      model=(i)
    )
  )

  rf100_models <- c(rf100_models, final_rf)

  final_wf <- workflow() %>%
    add_recipe(clus_rec) %>%
    add_model(final_rf)

  final_res_rf <- final_wf %>%
    last_fit(clus_split)

  model_auc <- as.double(final_res_rf$.metrics[[1]][2, 3])
  model_accuracy <- as.double(final_res_rf$.metrics[[1]][1, 3])
  rf100_aucs <- c(rf100_aucs, model_auc)
  rf100_accuracies <- c(rf100_accuracies, model_accuracy)

  # Save parameters and features of the best model from this 80-20 split iteration
  if(model_auc > bestmodelauc){
    bestmodel <- final_rf
    bestmodelauc <- model_auc
    besttrainset <- as.data.frame(clus_train)
    bestaccuracy <- model_accuracy
    bestfinal_res_rf <- final_res_rf
  }
  print(i)
  print(Sys.time() - old_time)
}

rf100_vimps1 <- rf100_vimps[-1,] %>% left_join(
  rf100_vimps[-1,] %>% group_by(Variable) %>% summarize(mean=mean(Importance)) %>%
    ungroup() %>% mutate(rank=rank(-mean)) %>% select(!mean), by="Variable"
) %>% mutate(nonst_labs = ifelse(Variable %in% paste0("d0_", c("il6", "il10", "tnfa", "fibrinogen", "ddimer", "ferritin", "crp", "ldh")),
                                 "Non-Standard Lab", "CBC and CMP"))

rf100_vimps1[rf100_vimps1$Variable=="d0_albumin",]$Variable <- "Albumin"
rf100_vimps1[rf100_vimps1$Variable=="d0_wbc",]$Variable  <- "WBC"
rf100_vimps1[rf100_vimps1$Variable=="d0_alk",]$Variable  <- "ALP"
rf100_vimps1[rf100_vimps1$Variable=="d0_ast",]$Variable  <- "AST"
rf100_vimps1[rf100_vimps1$Variable=="d0_crp",]$Variable  <- "CRP"
rf100_vimps1[rf100_vimps1$Variable=="d0_ddimer",]$Variable  <- "Ddimer"
rf100_vimps1[rf100_vimps1$Variable=="d0_ferritin",]$Variable  <- "Ferritin"
rf100_vimps1[rf100_vimps1$Variable=="d0_hb",]$Variable  <- "Hgb"
rf100_vimps1[rf100_vimps1$Variable=="d0_il10",]$Variable  <- "IL10"
rf100_vimps1[rf100_vimps1$Variable=="d0_il6",]$Variable  <- "IL6"
rf100_vimps1[rf100_vimps1$Variable=="d0_ldh",]$Variable  <- "LDH"
rf100_vimps1[rf100_vimps1$Variable=="d0_plt",]$Variable  <- "Plts"
rf100_vimps1[rf100_vimps1$Variable=="d0_tbr",]$Variable  <- "Tbili"
rf100_vimps1[rf100_vimps1$Variable=="d0_tnfa",]$Variable  <- "TNFa"

# Plot variable importance distributions
rf100_vimps2 <- rf100_vimps1 %>% mutate(Variable=fct_reorder(factor(Variable), -rank))
fig1g <- rf100_vimps2 %>% ggplot(aes(x=Variable, y=Importance, color=nonst_labs )) +
  geom_boxplot(aes(color=nonst_labs, color = after_scale(darken(color, .1, space = "HLS")),
                   fill = after_scale(desaturate(lighten(color, .8), .4))), outlier.shape = NA) +
  geom_jitter(alpha=0.4, width = 0.3) +
  coord_flip() +
  labs(y = "Variable Importance",
       x = "Laboratory Variable",
  )+
  guides(color = guide_legend(title = "")) +
  theme_minimal()+
  theme(
    axis.text.x = element_text(angle = 0, size=14),
    axis.text.y = element_text(angle = 0, size=13),
    axis.title.x = element_text(angle = 0, size=15, face="bold"),
    axis.title.y = element_text(size=15, face="bold"),
    title = element_text(angle = 0, size=15, face="bold"),
    legend.text = element_text(angle = 0, size=12, face="bold"),
    legend.key.size = unit(c(2,2), "cm"),
    legend.position = c(.8, 0.4)
  )



