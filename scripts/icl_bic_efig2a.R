library(tidyverse)
library(mclust)
library(ggsci)
library(ggrepel)
rm(list=ls())
load("output/scaledata.RData")
source("scripts/functions_constants.R")

# Preprocess and Generate Model from Development Cohort
df_dev <-  left_join(df_all, df_labs_all %>% select(!cohort), by="record_id") %>%
  filter(analysis_type==1) %>%
  filter(cohort=="MSK Development") %>%
  select(record_id, starts_with("d0_")) %>%
  select(record_id, contains(cluster_labs)) %>%
  drop_na() %>% column_to_rownames("record_id")

# Finite Mixture Modeling
BIC <- mclustBIC(df_dev %>% select(starts_with("d0_")))
ICL <- mclustICL(df_dev %>% select(starts_with("d0_")))
bicl <- bind_rows(
  as.data.frame(BIC[1:dim(BIC)[1], 1:dim(BIC)[2]]) %>%
  mutate(nclus=1:dim(BIC)[1]) %>% pivot_longer(!nclus, names_to = "shape", values_to = "BICL") %>%
    mutate(ic="BIC"),
  as.data.frame(ICL[1:dim(ICL)[1], 1:dim(ICL)[2]]) %>%
    mutate(nclus=1:dim(ICL)[1]) %>% pivot_longer(!nclus, names_to = "shape", values_to = "BICL") %>%
    mutate(ic="ICL"),
) %>% filter(ic=="ICL") %>% filter(!is.na(BICL)) %>%
  group_by(shape) %>% mutate(mnclus=max(nclus, na.rm=TRUE)) %>% ungroup() %>%
  mutate(nclus_label=ifelse(nclus==mnclus, shape, NA)) %>%
  mutate(shshape=factor(ifelse(shape=="VVV", "VVV", "Other")))

efig2a <- bicl %>%
  ggplot(aes(x=nclus, y=BICL, color=shape, label=nclus_label, shape=shshape)) + geom_point(aes( size=shshape)) + geom_line() +
  theme_minimal() + geom_label_repel() +
  labs(x = "Number of Clusters",
       y = "Integrated Complete Likelihood\nBayesian Information Criterion (ICL-BIC)",
       title = ""
  )+
  scale_x_discrete(limits=1:9, labels = 1:9) +
  theme(
    axis.text.x = element_text(angle = 0, size=13, face="bold"),
    axis.text.y = element_text(size=13, face="bold"),
    axis.title.x = element_text(angle = 0, size=18, face="bold", vjust=-1),
    axis.title.y = element_text(size=18, face="bold"),
    title = element_text(angle = 0, size=18, face="bold"),
    legend.text = element_text(angle = 0, size=14),
    legend.position = "none"
  )+
  guides(color=FALSE, label=FALSE, shape=FALSE)

