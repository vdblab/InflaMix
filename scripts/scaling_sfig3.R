rm(list=ls())
library(tidyverse)
library(ggsci)
library(moments)
library(ggrepel)

# Metadata and laboratory data for all patients across all treatment centers
df_all <- read.csv("data/df_all.csv")
labs_long_raw <- read.csv("data/labs_long_raw.csv")
source("scripts/functions_constants.R")

lbcl_dx <- c("DLBCL NOS", "High-grade B-cell lymphoma with MYC and BCL2 and/or BCL6 rearrangement",
             "High-grade B-cell lymphoma, NOS", "Primary mediastinal B-cell lymphoma",
             "T-cell rich DLBCL", "Primary cutaneous DLBCL, leg type", "Intravascular large B-cell lymphoma",
             "EBV-positive DLBCL", "DLBCL associated with chronic inflammation",
             "ALK positive large B-cell lymphoma", "Plasmablastic lymphoma",
             "Primary effusion lymphoma", "Large B-cell lymphoma arising in HHV8-associated multicentric Castleman Disease")

dev_products <- c("Axicabtagene ciloleucel", "Tisagenlecleucel", "Lisocabtagene maraleucel")


cluster_labs = c(
  "albumin",
  "alk",
  "ast",
  "ferritin",
  "hb",
  "ldh",
  "plt",
  "tbr",
  "il10",
  "il6",
  "tnfa",
  "crp",
  "ddimer",
  "wbc"
)

# raw = untransformed lab values, uln= upper limit of normal
# Normalize labs by upper limit of normal
df_labs_all_temp <- labs_long_raw %>% mutate(tr_uln=raw/uln) %>%
  select(!c(uln, raw)) %>%
  pivot_wider(
    names_from=c("day", "lab"),
    values_from = "tr_uln")


# Define cohorts.

# Derivation Cohort
cohort_dev <- (df_all %>%
                 filter(center=="MSK") %>% # treated at MSK
                 filter(dx.factor %in%  lbcl_dx) %>% # Large B-cell lymphoma diagnosis
                 filter(car_t_product_simple.factor %in% dev_products) %>% # Commercial products only
                 filter(dt_car_t == 1) %>% # Only patients treated before January 1, 2022
                 left_join(
                   df_labs_all_temp %>%
                     select(record_id, starts_with("d0_")) %>%
                     select(record_id, contains(cluster_labs)), by="record_id") %>%
                 select(record_id, starts_with("d0")) %>%
                 drop_na() # Only patients with completely available pre-infusion labs used for clustering
)$record_id

cohort_templbcl <- (df_all %>%
                      filter(center=="MSK") %>%
                      filter(dx.factor %in%  lbcl_dx) %>%
                      filter(!c(record_id %in% cohort_dev))
)$record_id # MSK patients with LBCL not in the derivation cohort

cohort_extlbcl <- (df_all %>%
                     filter(center!="MSK") %>%
                     filter(dx.factor %in%  lbcl_dx)
)$record_id # Patients at SMC or HMH with a LBCL diagnosis

cohort_mclfl <- (df_all %>%
                   filter(dx_simple.factor %in%  c("Follicular Lymphoma", "Mantle Cell Lymphoma"))
)$record_id # Patients with MCL or FL

# Annotate the cohort
df_labs_all_long1 <- labs_long_raw %>%
  mutate(cohort=ifelse(record_id %in% cohort_dev, "dev",
                       ifelse(record_id %in% cohort_templbcl, "templbcl",
                              ifelse(record_id %in% cohort_extlbcl, "extlbcl", "mclfl"))))

# Evaluate the skew for labs in the derivation cohort
# Compare the skew after log transformation
skewtest1 <- df_labs_all_long1 %>%
  mutate(tr_uln = raw/uln) %>%
  mutate(tr_uln_log=log(raw+1e-3, base=10)-log(uln, base=10)) %>%
  filter(lab!="fibrinogen") %>%
  filter(day=="d0") %>%
  filter(cohort %in% c("dev")) %>% group_by(cohort, lab) %>%
  summarize(skew=skewness(tr_uln), logskew=skewness(tr_uln_log)) %>% ungroup() %>%
  mutate(skewdiff=logskew/skew) %>% mutate(skewthresh=factor(ifelse(skew>1 & (skewdiff>0 & skewdiff < 0.10), "Yes", "No")))

# Graphical comparison of skew before and after log transformation for each lab
sfig3 <- ggplot(data=skewtest1, aes(x=skew, y=logskew, label=lab, color=skewthresh)) + geom_point(size=5) + geom_label_repel(size=6)+
  labs(x="Skew",
       y="Log Skew",
       color="Skew > 1 and Reduced > 90%\nBy Log Transformation")+
  scale_color_npg()+
  theme_minimal()

# We will log transform ferritin and crp based on our analysis of labs in the derivation cohort
log_labs <- c("ferritin","crp")

# Distribution parameters by lab in the derivation cohort
sc_devcoh_d0 <- df_labs_all_long1 %>% mutate(tr_uln=raw/uln) %>%
  filter(lab %in% cluster_labs) %>%
  filter(record_id %in% cohort_dev) %>% filter(day=="d0") %>%
  mutate(tr_uln=ifelse(lab %in% log_labs, log(tr_uln, 10), tr_uln)) %>% # log labs
  group_by(lab) %>% mutate(
                           tr_uln_mean=mean(tr_uln, na.rm=TRUE),
                           tr_uln_sd=sd(tr_uln, na.rm=TRUE)
  ) %>% ungroup() %>%
  select(lab, tr_uln_mean, tr_uln_sd) %>%
  distinct()

# For all patients across all cohorts - normalize labs by ULN and then scale them to laboratory distributions from the derivation cohort
df_labs_all_long2 <- df_labs_all_long1 %>%
  filter(lab %in% cluster_labs) %>%
  mutate(tr_uln=raw/uln) %>%
  mutate(tr_uln=ifelse(lab %in% log_labs, log(tr_uln, 10), tr_uln)) %>% # log labs
  left_join(sc_devcoh_d0, by="lab") %>%
  mutate(scld = (tr_uln-tr_uln_mean)/tr_uln_sd ) %>%
  select(!c(tr_uln,tr_uln_mean, tr_uln_sd))

# Re-organize the laboratory tables into a wide-format
df_labs_all_long <- df_labs_all_long2 %>%
  mutate(abuln=ifelse(raw > uln, 1, 0))

# Wide format for normalized and scaled data
df_labs_all1 <- df_labs_all_long %>%
  select(!c(uln, raw, abuln)) %>%
  pivot_wider(
    names_from=c("day", "lab"),
    values_from = "scld")

# Wide format for raw lab values - merge with upper limits of normal.
df_labs_all2 <- df_labs_all_long %>%
  select(!c(uln, scld, abuln)) %>%
  pivot_wider(
    names_from=c("day", "lab"),
    values_from = "raw"
  ) %>% column_to_rownames("record_id") %>% dplyr::rename_all(~paste0("noLCS_", .)) %>%
  rownames_to_column("record_id") %>%
  left_join(df_labs_all_long%>%filter(lab=="ldh" & day=="preld") %>% select(record_id, abuln), by="record_id") %>%
  rename(bin_preld_ldh=abuln) %>%
  left_join(df_labs_all_long%>%filter(lab %in% cluster_labs & day=="preld") %>% select(record_id, lab, uln) %>%
              pivot_wider(names_prefix = "uln_preld_", names_from = "lab", values_from = "uln"), by="record_id") %>%
  left_join(df_labs_all_long%>%filter(lab %in% cluster_labs & day=="d0") %>% select(record_id, lab, uln) %>%
              pivot_wider(names_prefix = "uln_d0_", names_from = "lab", values_from = "uln"), by="record_id")

# Combine scaled, raw, and uln labs all together in a wide format.
# Center validation = SMC+HMH LBCL Validation
# Non-Hodgkin Lymphoma Validation = MCL+FL Validation
df_labs_all <- inner_join(df_labs_all1, df_labs_all2, by="record_id") %>%
  mutate(cohort=ifelse(cohort=="dev", "MSK Development", ifelse(
    cohort=="templbcl", "MSK Validation", ifelse(
      cohort=="extlbcl", "Center Validation", "Non-Hodgkin Lymphoma Validation")))) %>%
  mutate(cohort=factor(cohort, levels=c("MSK Development", "MSK Validation", "Center Validation", "Non-Hodgkin Lymphoma Validation"))) %>%
  mutate(record_id=record_id)

# Rename centers and cohorts for the metadata table
df_all <- df_all %>%
  mutate(center=ifelse(center=="MSK", "Memorial Sloan Kettering", ifelse(
    center=="SHE", "Sheba Medical Center", ifelse(
      center=="HMH", "Hackensack Meridian Health", center
    )
  ))) %>%
  mutate(cohort=ifelse(cohort=="dev", "MSK Development", ifelse(
    cohort=="templbcl", "MSK Validation", ifelse(
      cohort=="extlbcl", "Center Validation", "Non-Hodgkin Lymphoma Validation")))) %>%
  mutate(cohort=factor(cohort, levels=c("MSK Development", "MSK Validation", "Center Validation", "Non-Hodgkin Lymphoma Validation"))) %>%
  mutate(record_id=record_id)

# Save processed/scaled data objects for downstream analysis.
save(sfig3, file="figures/sfig3.RData")
save(df_labs_all, df_all, file="output/scaledata.RData")
