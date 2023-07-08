library(tidyverse)
library(mclust)
library(ggsci)
library(ggrepel)

rm(list=ls())
load(bstfun::here_data("analysis/scaledata.RData"))
load(bstfun::here_data("analysis/model.RData"))
source("scripts/analysis/functions_constants.R")

featdist_df <- left_join(df_all, df_labs_all %>% select(!cohort) %>% select(c(record_id, starts_with("d0_"))), by="record_id") %>%
  filter(analysis_type==1)  %>%
  filter(cohort=="MSK Development")

wcox <- reformat_table_tp(df=featdist_df, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>% select(c(cluster, contains("d0_"))) %>%
  pivot_longer(!cluster, names_to = "labs", values_to = "value") %>%
  filter(labs %in% c("d0_il6", "d0_crp", "d0_ldh", "d0_tbr", "d0_wbc", "d0_hb")) %>%
  group_by(labs) %>% do(w = wilcox.test(value~cluster, data=., paired=FALSE)) %>%
  summarize(labs, wilcox = w$p.value) %>% ungroup() %>% mutate(wilcox=p.adjust(wilcox, method="fdr"))

efig2df <- reformat_table_tp(df=featdist_df, gmm_mu=dev_mu, gmm_pro=dev_pro, gmm_sig=dev_sig, inflammClus = dev_inflammClus, least_inflammClus = dev_LEASTinflammClus, tp="d0_") %>%
  group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>% select(c(cluster, contains("d0_"))) %>%
  pivot_longer(!cluster, names_to = "labs", values_to = "value") %>%
  filter(labs %in% c("d0_il6", "d0_crp", "d0_ldh", "d0_tbr", "d0_wbc", "d0_hb")) %>%
  group_by(labs, cluster) %>%
  mutate(value.show = as.numeric(  # so ggplot doesn't complain about alpha being discrete
    between(value,
            quantile(value)[2] - 1.5*IQR(value),
            quantile(value)[4] + 1.5*IQR(value)))) %>% ungroup() %>%
  filter(value.show != 0) %>%
  left_join(wcox, by="labs") %>%
  mutate(imp_labs=ifelse(labs %in% c("d0_il6", "d0_crp", "d0_ldh"), "Important Lab", "Non-Important Lab"))

efig2df[efig2df$labs=="d0_il6",]$labs <- "Interleukin-6"
efig2df[efig2df$labs=="d0_crp",]$labs <- "C-Reactive Protein"
efig2df[efig2df$labs=="d0_ldh",]$labs <- "Lactate Dehydrogenase"
efig2df[efig2df$labs=="d0_hb",]$labs <- "Hemoblogin"
efig2df[efig2df$labs=="d0_wbc",]$labs <- "White Blood Cells"
efig2df[efig2df$labs=="d0_tbr",]$labs <- "Total Bilirubin"

lablist <- unique(efig2df$labs)
lablist <- factor(lablist, levels=c("Interleukin-6",
                                    "C-Reactive Protein",
                                    "Lactate Dehydrogenase",
                                    "Hemoblogin",
                                    "White Blood Cells",
                                    "Total Bilirubin"

))

labcomp <- function(x, ll, i, corner) {
  y1 <- ""
  if(corner==TRUE){
    y1 <- "Center-Scaled Value"
  }
  y <-
    x  %>% filter(labs==lablist[i]) %>%
    mutate(labs=paste0(labs, "\np = ", signif(wilcox,digits = 2))) %>%
    ggplot(aes(x=cluster, y=value, color=cluster))+ geom_boxplot(outlier.shape=NULL, outlier.color = NULL) +
    geom_point(aes(x=cluster), position=position_jitterdodge(jitter.width = 1))+
    theme_minimal() +
    facet_wrap(.~labs, scales = "free_y") +
    scale_fill_manual(values=colorclusters(length(dev_pro)))+
    scale_color_manual(values=colorclusters(length(dev_pro)))+
    labs(
         y=y1,
         color="",
         x=""
         ) +
    #guides(color=guide_legend(title="", title.position = "", override.aes = list(size = 6))) +
    theme(
      strip.text = element_text(size=14, face="bold"),
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, size=13),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size=20, face="bold"),
      legend.text  = element_text(size=20, face="bold"),
    ) +
    guides(color=guide_legend(title="",nrow = 1 , override.aes = list(size = 6)))


  return(y)
}

efig2b <- labcomp(x=efig2df, ll=lablist, i=1, corner=TRUE)

legend <- get_legend(
  # create some space to the left of the legend
  efig2b
)

efig2b <- labcomp(x=efig2df, ll=lablist, i=6, corner=TRUE) + theme(legend.position="none")
efig2c <- labcomp(efig2df,lablist, i=1, FALSE)  + theme(legend.position="none")
efig2d <- labcomp(efig2df,lablist, i=3, FALSE) + theme(legend.position="none")
efig2e <- labcomp(efig2df,lablist, i=2, TRUE) + theme(legend.position="none")
efig2f <- labcomp(efig2df,lablist, i=5, FALSE) + theme(legend.position="none")
efig2g <- labcomp(efig2df,lablist, i=4, FALSE) + theme(legend.position="none")


top <- plot_grid(efig2b, NULL,
                 efig2c, NULL,
                 efig2d,
          labels=c("b", "", "c", "", "d"),
          rel_widths = c(1, 0.05,
                         1, 0.05,
                         1),
          nrow=1,
          scale=0.95,
          label_size = 22)

bottom <- plot_grid(efig2e, NULL,
                    efig2f, NULL,
                    efig2g,
                    labels=c("e", "", "f", "", "g"),
                    rel_widths = c(1, 0.05,
                                   1, 0.05,
                                   1),
                    nrow=1,
                    scale=0.95,
                    label_size = 22)
title1 <- ggdraw() +
  draw_label(
    "Important by Random Forest",
    fontface = 'bold',
    x = 0.35,
    hjust = 0,
    size=20
  )

title2 <- ggdraw() +
  draw_label(
    "Less Important by Random Forest",
    fontface = 'bold',
    x = 0.35,
    hjust = 0,
    size=20
  )

efig2_bg <- plot_grid(top, NULL, bottom, labels=c("", "", ""), rel_heights = c(1, 0.025, 1), ncol=1)
