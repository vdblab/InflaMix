library(tidyverse)
library(mclust)
library(mvtnorm)
library(survival)
library(broom)
library(ggsurvfit)

# Obtain color scales
colorclusters <- function(nclus){
  fullspectrum=c("#9467BDFF", "#2CA02CFF", "#8C564BFF", "#E377C2FF", "#7F7F7FFF", "#BCBD22FF")
  colorscale=c("#1F77B4FF",  "#FF7F0EFF")
  if(nclus>2){
    colorscale=c("#1F77B4FF", fullspectrum[1:(nclus-2)], "#FF7F0EFF")
  }
  return(colorscale)
}

# Function to assign a cluster from partial labs
assign_clust <- function(x, mu, sig, pro){
  misfeats <- colnames(x)[!c(colnames(x) %in% rownames(as.data.frame(t(x)) %>% filter(!is.na(.))))]
  x1 <- x %>% select(!contains(misfeats))
  mu1 <- mu[!(rownames(mu) %in% misfeats),]
  sig1 <- sig[!(rownames(sig) %in% misfeats),!(rownames(sig) %in% misfeats),]
  cnum <- length(pro)
  densfun <- 1:cnum
  if(sum(dim(x1))==2){
    mu1 <- as.matrix(data.frame(t(mu1)))
    colnames(mu1) <- NULL
    rownames(mu1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    sig1 <- array(sig1, dim=c(1,1,2))
    colnames(sig1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    rownames(sig1) <- colnames(x)[!c(colnames(x) %in% misfeats)]
    return(NaN)
  }
  for(i in 1:cnum){
    densfun[i] <- dmvnorm(x1, mean = mu1[,i], sigma = sig1[,,i], log = FALSE)
  }
  return( (pro * densfun) / sum(pro * densfun) )
}

# Function to reformat data table to include cluster assignments with tau probabilities for any timepoint
reformat_table_tp <- function(df, gmm_mu, gmm_sig, gmm_pro, inflammClus, least_inflammClus, tp){
  df_2clust_d0labs <- df %>% select(record_id, starts_with(tp)) %>%
    select(record_id, contains(cluster_labs)) %>% column_to_rownames("record_id") %>%
    rename_all(~stringr::str_replace(.,paste0("^", tp),"")) %>%
    rename_with(~ paste0("d0_", .))
  reformat.dat <- data.frame(cluster=1, tau=0.5, record_id="id_test")
  for(k in 1:dim(df)[1]){
    if(num_na <- sum(is.na(df_2clust_d0labs[k,]))==dim(df_2clust_d0labs)[2]){
      reformat.dat <- reformat.dat %>% bind_rows(data.frame(
        cluster = 1:length(gmm_pro),
        tau = NA,
        record_id=rownames(df_2clust_d0labs)[k]
      ))
    } else(
      reformat.dat <- reformat.dat %>% bind_rows(data.frame(
        cluster = 1:length(gmm_pro),
        tau = assign_clust(x=df_2clust_d0labs[k,], gmm_mu, gmm_sig, gmm_pro),
        record_id=rownames(df_2clust_d0labs)[k]
      ))
    )
  }
  if(length(gmm_pro)>2){
    clus_recode <- reformat.dat %>% select(cluster) %>% filter(!(cluster %in% c(inflammClus, least_inflammClus))) %>%
      distinct() %>% mutate(newclus=1:n())
    reformat.dat <- reformat.dat[-1,] %>% left_join(clus_recode, by="cluster") %>%
      inner_join(df, by="record_id") %>%
      mutate(cluster=factor(ifelse(cluster==inflammClus, "Inflammatory",
                                   ifelse(cluster==least_inflammClus,
                                          "Non-Inflammatory",
                                          ifelse(as.integer(newclus)==1, "Neutral Cluster",
                                                 paste0("Neutral Cluster-", as.integer(newclus) ))))))%>%
      mutate(cluster=factor(cluster, levels=rev(levels(cluster)))) %>%
      mutate(cluster=relevel(cluster, "Non-Inflammatory")) %>%
      filter(tau!=0) %>% filter(!is.nan(tau))
  } else{
    reformat.dat <- reformat.dat[-1,] %>%
      inner_join(df, by="record_id") %>%
      mutate(cluster=factor(ifelse(cluster==inflammClus, "Inflammatory","Non-Inflammatory"))) %>%
      mutate(cluster=factor(cluster, levels=rev(levels(cluster)))) %>%
      mutate(cluster=relevel(cluster, "Non-Inflammatory")) %>%
      filter(tau!=0) %>% filter(!is.nan(tau))
  }
  return(reformat.dat)
}

# Plot barplots for binary outcomes
clust_bar_plot <- function(df, tp_pre, list_res, metric, xmi, ymi, sz) {
  colorscale <- colorclusters(length(unique(df$Cluster)))
  df_restox <- bind_rows(
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, everCR_100) %>% reframe(n=n()) %>%
      filter(!is.na(everCR_100)) %>% group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(everCR_100=="CR") %>%
      rename(feat=everCR_100) %>% mutate(feat_name="CR"),
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, crs24) %>% reframe(n=n()) %>%
      filter(!is.na(crs24)) %>%group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(crs24=="CRS >1") %>%
      rename(feat=crs24) %>% mutate(feat_name="CRS"),
    df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      group_by(Cluster, icans24) %>% reframe(n=n()) %>%
      filter(!is.na(icans24)) %>%group_by(Cluster) %>%
      mutate(prop=n/sum(n), tot=sum(n)) %>% ungroup() %>% filter(icans24=="ICANS >1") %>%
      rename(feat=icans24) %>% mutate(feat_name="ICANS")
  ) %>% filter(feat_name %in% list_res)

  trp <- df_restox %>%
    ggplot(aes(x=feat_name, y=prop, fill=Cluster)) +
    geom_bar(stat="identity", color="black",
             position=position_dodge())+
    geom_label(aes(y=prop/2, group=Cluster, label=paste0(round(prop,2),"\n",n,"/",tot)), color="black", fill="white",
               fontface="bold", position = position_dodge(width = 0.9), size=sz) +
    guides(fill = FALSE) +
    theme_minimal()+
    theme(
      axis.text.x = element_blank(),
      axis.text.y = element_text(angle = 0, size=15),
      axis.title.x = element_text(angle = 0, size=15, face="bold"),
      axis.title.y = element_text(size=16, face="bold"),
    ) +
    labs(x = "Non-Infl.               Infl.",
         y = paste0("Proportion ", list_res),
         title="",
    )+
    scale_fill_manual(values=colorscale) +
    annotate(geom = "text", x=xmi, y=ymi, label = paste0(
      "Adj. OR (95% CI): ", metric$expEstimate, " (", metric$low_ci, " - ", metric$high_ci, ")\n", "p ",  cpval(metric$pvalue)
    ), size=sz)
  return(trp)
}

# Plot KM curves for survival outcomes
ssr_survplot <- function(df, event, timemax, qmonth, metric, wght, sz, xmi, ymi, labl){

  if(wght==FALSE){
    df <- df %>% group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(tau=1)
  }

  if(event=="pfs"){
    event <- "PFS"
    sp <- survfit(Surv(tt_pfs_m, ev_pfs) ~ Cluster, data = df, weights=tau)
  } else {
    if(event=="os") {
      event <- "OS"
      sp <- survfit(Surv(tt_os_m, ev_os) ~ Cluster, data = df, weights=tau)
    } else {event <- ""}
  }

  shapes <- data.frame(cens=1:4, cshape=c(3, 3, 3, 3))
  sdf <- data.frame(time = sp$time,
                    prob = sp$surv,
                    arisk= round(sp$n.risk),
                    cens= round(sp$n.censor),
                    strata = c(rep(names(sp$strata)[1], sp$strata[1]),
                               rep(names(sp$strata)[2], sp$strata[2])),
                    events=sp$n.event
  ) %>% mutate(cens=ifelse(cens==0, NA, cens))%>% mutate(cens=ifelse(time>timemax, NA, cens)) %>%
    left_join(shapes, by="cens")


  sdf_inf <- sdf %>% filter(!grepl("Non-", strata))
  sdf_inf <- bind_rows(
    data.frame(time=0, prob=1, arisk=sdf_inf$arisk[1], cens=NA, strata=sdf_inf$strata[1], events=0, cshape=NA),
    sdf_inf
  )

  if(max(sdf_inf$time)<timemax){
    lastcens=1
    if(is.na(sdf_inf[which.max(sdf_inf$time),]$cens)){lastcens=0}
    sdf_inf <- bind_rows(sdf_inf,
                     data_frame(time=timemax+0.1,
                                prob=NA,
                                arisk=min(sdf_inf$arisk)-(lastcens+sdf_inf[which.max(sdf_inf$time),]$events),
                                cens=NA,
                                strata=sdf_inf[which.max(sdf_inf$time),]$strata,
                                events=sdf_inf[which.max(sdf_inf$time),]$events,
                                cshape=NA
                                ))
  }

  sdf_ninf <- sdf %>% filter(grepl("Non-", strata))
  sdf_ninf <- bind_rows(
    data.frame(time=0, prob=1, arisk=sdf_ninf$arisk[1], cens=NA, strata=sdf_ninf$strata[1], events=0, cshape=NA),
    sdf_ninf
  )
  if(max(sdf_ninf$time)<timemax){
    lastcens=1
    if(is.na(sdf_ninf[which.max(sdf_ninf$time),]$cens)){lastcens=0}
    sdf_ninf <- bind_rows(sdf_ninf,
                     data_frame(time=timemax+0.1,
                                prob=NA,
                                arisk=min(sdf_ninf$arisk)-(lastcens+sdf_ninf[which.max(sdf_ninf$time),]$events),
                                cens=NA,
                                strata=sdf_ninf[which.max(sdf_ninf$time),]$strata,
                                events=sdf_ninf[which.max(sdf_ninf$time),]$events,
                                cshape=NA
                     ))
  }

  sdf <- bind_rows(sdf_ninf, sdf_inf)

  timeseries <- seq(0, timemax, by=timemax/qmonth)
  tmIND1 <- 1:length(timeseries)
  tmIND2 <- 1:length(timeseries)
  for(i in 2:length(timeseries)){
    x1 <- timeseries[i]-sdf[sdf$strata==names(sp$strata)[1],]$time
    x1 <- abs((x1 < 0) * x1)
    x1[x1==0] <- 1e8
    x2 <- timeseries[i]-sdf[sdf$strata==names(sp$strata)[2],]$time
    x2 <- abs((x2 < 0) * x2)
    x2[x2==0] <- 1e8

    tmIND1[i] <- which.min(x1)
    tmIND2[i] <- which.min(x2)
  }

  timedat <- sp$time
  strat1 <- sp$strata[1]
  strat2 <- sp$strata[2]
  timind1 <- which.min(abs(timedat[1:strat1]-20))
  timind2 <- which.min(abs(timedat[-(1:strat1)]-20))+strat1
  survind1 <- sp$surv[timind1]+0.1
  survind2 <- sp$surv[timind2]+0.1

  survpl <- ggplot(sdf, aes(time, prob, color=strata)) + geom_point(aes(shape=cshape), size=4) + geom_step(linewidth=1.5) +
    scale_shape_identity() +
    coord_cartesian(ylim = c(0, 1.05), expand = FALSE, clip = "off") +
    scale_color_manual(values=rev(colorclusters(2)))+
    theme_bw()+
    theme(
      axis.text.x = element_text(angle = 0, size=18),
      axis.text.y = element_text(angle = 0, size=18),
      axis.title.x = element_text(angle = 0, size=18, face="bold", vjust=-2),
      axis.title.y = element_text(size=18, face="bold", vjust=2),
      title = element_text(angle = 0, size=18, face="bold"),
      legend.text = element_text(angle = 0, size=14),
      legend.position = "none",
      plot.margin = unit(c(1, 1, 5, 1), "lines")
    ) +
    annotate(geom = "text", x = 0,
             y = -0.17, label = "At Risk", size =6, fontface="bold") +
    annotate(geom = "text", x = timeseries,
             y = -0.24, label = sdf$arisk[sdf$strata==names(sp$strata)[1]][tmIND1], size =7, color=colorclusters(2)[1])+

    annotate(geom = "text", x = timeseries,
             y = -0.31, label = sdf$arisk[sdf$strata==names(sp$strata)[2]][tmIND2], size =7, color=colorclusters(2)[2])+
    annotate("label", x=20, y=survind1, label="Non-Inflammatory", size=sz/3, color=colorclusters(2)[1])+
    annotate("label", x=20, y=survind2, label="Inflammatory", size=(sz/3)+0.5, color=colorclusters(2)[2])+
    scale_x_continuous(limits=c(0, 24), expand=c(-10,0), oob=scales::squish) +
    labs(
      title ="",
      y = paste0(event, " Probability"),
      x="Months from Infusion"
    ) +
    annotate(geom = "text", x=xmi, y=ymi, label = paste0(
      "Adj. HR (95% CI): ", metric$expEstimate, " (", metric$low_ci, " - ", metric$high_ci, ")\n", "p ",  cpval(metric$pvalue)
      ), size=8.5)

  pl_fin <- plot_grid(survpl, labels=c(labl), scale=1, label_size = 30)
  return(pl_fin)
}

# Calculate Inferences for Validation
clust_metrics <- function(df, dxres, covar, lcovar, res_covar, coxmodels, aname, wght){

  if(!grepl("trans", aname)){
  df_inflam <- df %>% filter(Cluster=="Inflammatory")
  mm_cindex_os_temp <- survcomp::concordance.index(df_inflam$tau, surv.time = df_inflam$tt_os_m, surv.event = df_inflam$ev_os, alternative = "g")
  mm_cindex_os <- data.frame(cindex=mm_cindex_os_temp$c.index,
                    lowci=mm_cindex_os_temp$lower,
                    hici=mm_cindex_os_temp$upper,
                    pval=mm_cindex_os_temp$p.value)

  mm_cindex_pfs_temp <- survcomp::concordance.index(df_inflam$tau, surv.time = df_inflam$tt_pfs_m, surv.event = df_inflam$ev_pfs, alternative = "g")
  mm_cindex_pfs <- data.frame(cindex=mm_cindex_pfs_temp$c.index,
                    lowci=mm_cindex_pfs_temp$lower,
                    hici=mm_cindex_pfs_temp$upper,
                    pval=mm_cindex_pfs_temp$p.value)
  }

  if(covar==lcovar){
    dev_df_hc <- coxmodels$oldf
    contr <- c(-1, 1)

    df_hc_crp <- df %>%
      mutate(cluster=Cluster) %>%
      filter(cluster=="Inflammatory") %>%
      group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(record_id=record_id)

    diff1 <- concordance(coxmodels$crp_cox0os, coxmodels$crp_coxos, newdata=df_hc_crp)
    dtest <- contr %*% coef(diff1)
    dvar <- contr %*% vcov(diff1) %*% contr
    #sd=sqrt(dvar)
    #z=dtest/sqrt(dvar)
    crp_cox_cindex_os <- as.data.frame(t(c(cindex=round(coef(diff1)[2],2), crp_coxph_cindex_contrast=signif(dtest,2),
                                           crp_cindex_lowci = round(coef(diff1)[2]-1.96*sqrt(diff1$var[2,2]),2),
                                           crp_cindex_hici = round(coef(diff1)[2]+1.96*sqrt(diff1$var[2,2]),2),
                                           crp_pvalue_cindex_contrast=signif(2*(1-pnorm(dtest/sqrt(dvar))),2))))
    colnames(crp_cox_cindex_os)[1] <- "crp_coxph_cindex"
    colnames(crp_cox_cindex_os)[3] <- "crp_cindex_lowci"
    colnames(crp_cox_cindex_os)[4] <- "crp_cindex_hici"

    diff1 <- concordance(coxmodels$crp_cox0pfs, coxmodels$crp_coxpfs, newdata=df_hc_crp)
    dtest <- contr %*% coef(diff1)
    dvar <- contr %*% vcov(diff1) %*% contr
    crp_cox_cindex_pfs <- as.data.frame(t(c(cindex=round(coef(diff1)[2],2), crp_coxph_cindex_contrast=signif(dtest,2),
                                            crp_cindex_lowci = round(coef(diff1)[2]-1.96*sqrt(diff1$var[2,2]),2),
                                            crp_cindex_hici = round(coef(diff1)[2]+1.96*sqrt(diff1$var[2,2]),2),
                                            crp_pvalue_cindex_contrast=signif(2*(1-pnorm(dtest/sqrt(dvar))),2))))
    colnames(crp_cox_cindex_pfs)[1] <- "crp_coxph_cindex"
    colnames(crp_cox_cindex_pfs)[3] <- "crp_cindex_lowci"
    colnames(crp_cox_cindex_pfs)[4] <- "crp_cindex_hici"

    df_hc_crpferr <- df    %>%
      mutate(cluster=Cluster) %>%
      filter(cluster=="Inflammatory") %>%
      group_by(record_id) %>% slice(which.max(tau)) %>% ungroup() %>%
      mutate(record_id=record_id)

    diff1 <- concordance(coxmodels$crpferr_cox0os, coxmodels$crpferr_coxos, newdata=df_hc_crpferr)
    dtest <- contr %*% coef(diff1)
    dvar <- contr %*% vcov(diff1) %*% contr
    crpferr_cox_cindex_os <- as.data.frame(t(c(cindex=round(coef(diff1)[2],2), crpferr_coxph_cindex_contrast=signif(dtest,2),
                                               crpferr_cindex_lowci = round(coef(diff1)[2]-1.96*sqrt(diff1$var[2,2]),2),
                                               crpferr_cindex_hici = round(coef(diff1)[2]+1.96*sqrt(diff1$var[2,2]),2),
                                               crpferr_pvalue_cindex_contrast=signif(2*(1-pnorm(dtest/sqrt(dvar))),2))))
    colnames(crpferr_cox_cindex_os)[1] <- "crpferr_coxph_cindex"
    colnames(crpferr_cox_cindex_os)[3] <- "crpferr_cindex_lowci"
    colnames(crpferr_cox_cindex_os)[4] <- "crpferr_cindex_hici"

    diff1 <- concordance(coxmodels$crpferr_cox0pfs, coxmodels$crpferr_coxpfs, newdata=df_hc_crpferr)
    dtest <- contr %*% coef(diff1)
    dvar <- contr %*% vcov(diff1) %*% contr
    crpferr_cox_cindex_pfs <- as.data.frame(t(c(cindex=round(coef(diff1)[2],2), crpferr_coxph_cindex_contrast=signif(dtest,2),
                                                crpferr_cindex_lowci = round(coef(diff1)[2]-1.96*sqrt(diff1$var[2,2]),2),
                                                crpferr_cindex_hici = round(coef(diff1)[2]+1.96*sqrt(diff1$var[2,2]),2),
                                                crpferr_pvalue_cindex_contrast=signif(2*(1-pnorm(dtest/sqrt(dvar))),2))))
    colnames(crpferr_cox_cindex_pfs)[1] <- "crpferr_coxph_cindex"
    colnames(crpferr_cox_cindex_pfs)[3] <- "crpferr_cindex_lowci"
    colnames(crpferr_cox_cindex_pfs)[4] <- "crpferr_cindex_hici"
  } else{
    crp_cox_cindex_os <- as.data.frame(t(c(crp_coxph_cindex=NA, crp_coxph_cindex_contrast=NA, crp_cindex_lowci=NA, crp_cindex_hici=NA,crp_pvalue_cindex_contrast=NA)))
    crp_cox_cindex_pfs<- as.data.frame(t(c(crp_coxph_cindex=NA, crp_coxph_cindex_contrast=NA, crp_cindex_lowci=NA, crp_cindex_hici=NA,crp_pvalue_cindex_contrast=NA)))
    crpferr_cox_cindex_os<- as.data.frame(t(c(crpferr_coxph_cindex=NA, crpferr_coxph_cindex_contrast=NA, crpferr_cindex_lowci=NA, crpferr_cindex_hici=NA,crpferr_pvalue_cindex_contrast=NA)))
    crpferr_cox_cindex_pfs<- as.data.frame(t(c(crpferr_coxph_cindex=NA, crpferr_coxph_cindex_contrast=NA, crpferr_cindex_lowci=NA, crpferr_cindex_hici=NA,crpferr_pvalue_cindex_contrast=NA)))
  }

  if(wght==FALSE){
    if(!grepl("trans", aname)){
      df <- df %>% group_by(record_id) %>%
        dplyr::slice(which.max(tau)) %>%ungroup() %>% mutate(tau=1)
    }
    df <- df %>% mutate(tau=1)
  } else {
    df <- df
  }

  df_metrics <- bind_rows(
    tidy(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial())) %>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="No CR") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial()))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(glm(as.formula(paste0(dxres, "~Cluster", res_covar)),weights = tau,data=df,family=binomial()))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=res_covar),

    tidy(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))%>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="PFS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_pfs_m, ev_pfs)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=covar) %>% mutate(mm_cindex=mm_cindex_pfs$cindex, mm_cindex_low=mm_cindex_pfs$lowci, mm_cindex_hi=mm_cindex_pfs$hici, mm_cindex_p=mm_cindex_pfs$pval) %>% bind_cols(crp_cox_cindex_pfs, crpferr_cox_cindex_pfs),

    tidy(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))%>%
      filter(term=="ClusterInflammatory") %>% select(estimate, p.value) %>% mutate(outcome="OS") %>% rename(expEstimate=estimate, pvalue=p.value) %>%
      mutate(low_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "2.5 %"])) %>%
      mutate(high_ci=exp(confint(coxph(as.formula(paste0("Surv(tt_os_m, ev_os)~Cluster", covar)), data = df, weights=df$tau))["ClusterInflammatory", "97.5 %"])) %>%
      mutate(covariates=covar) %>% mutate(mm_cindex=mm_cindex_os$cindex, mm_cindex_low=mm_cindex_os$lowci, mm_cindex_hi=mm_cindex_os$hici, mm_cindex_p=mm_cindex_os$pval) %>% bind_cols(crp_cox_cindex_os, crpferr_cox_cindex_os)

  )  %>% mutate(low_ci=round(low_ci, 2), high_ci=round(high_ci, 2)) %>% mutate(expEstimate=round(exp(expEstimate), 2), pvalue=signif(pvalue, 2)) %>% mutate(analysis=aname)
  return(df_metrics)
}

# Generate Table Grob
metric_grob_gen <- function(df, metricsi, oh, sz) {
  if(oh=="OR"){
    ohtext = "Adj. OR (95% CI)"
  } else {ohtext ="Adj. HR (95% CI)"}

  dfi <- data.frame(cluster=c("Non-Infl.", "Infl."),
                    Adj_OR = c("Reference", paste0(metricsi$expEstimate, " (", metricsi$low_ci, " - ", metricsi$high_ci, ")")),
                    pvalue= c("", cpval(metricsi$pvalue))) %>%
    tableGrob(rows=NULL, cols=c("", ohtext, "p value"), theme=ttheme_minimal(base_size=sz,
                                                                             core = list(padding=unit(c(2, 2), "mm"))))
  gi <- gtable_add_grob(dfi,
                        grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                        t = 2, b = nrow(dfi), l = 1, r = ncol(dfi)) %>%
    gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 2)),
                    t = 1, l = 1, r = ncol(dfi)) %>%
    gtable_add_grob(grobs = rectGrob(gp = gpar(fill = NA, lwd = 3)),
                    t = 1, b=3, l = 2, r = ncol(dfi))
}

cpval <- function(x) {
  y <- as.character(x)
  if(x >  0.1){y <- ">0.1"}
  if(x <  0.01){y <- "<0.01"}
  if(x <  0.001){y <- "<0.001"}
  if(x <  0.0001){y <- "<0.0001"}
  return(y)
}

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

focus_labs <- paste0("d0_", cluster_labs)[!c(paste0("d0_", cluster_labs) %in% c(
  "d0_albumin",
  "d0_alk",
  "d0_ast",
  "d0_hb",
  "d0_ldh",
  "d0_crp",
  ""
))]

trans_covar <- "+age+bridge.factor+bin_preld_ldh+dx_simple.factor+costim+primary_ref.factor"
lbcl_covar="+age+costim+ bin_preld_ldh+primary_ref.factor"
lymphoma_covar="+ age + bin_preld_ldh + dx_simple.factor+primary_ref.factor"
costim_covar="+ age + primary_ref.factor + bin_preld_ldh  + car_t_product_simple.factor"
product_covar="+ age + primary_ref.factor + bin_preld_ldh  + dx_simple.factor"
lymphoma_dxres = "everCR_100"
