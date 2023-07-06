######################################################
# simulation figures
######################################################
source("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photometry_sim_fLME_fn_multi.R") # simulations
wd <- "/Users/loewingergc/Desktop/Research/photometry_sims" # simulation_delays" # sim_delay2   #"/Users/loewingergc/Desktop/Research/science_simulates"
setwd(wd)
library(dplyr)
library(latex2exp)
library(kableExtra)
library(tidyverse)

nVec <- c(4:8) # just compare subset
fixed_smooth <- TRUE
knots_div <- 2
indp <- c(FALSE) 
resid_var_subj <- TRUE
fix_scl <- 1
fix_shift <- 0
resid_scl <- 5 
rwdLen_vec <- 2
ls_rmse <- ls_rmse_avg <- ls_rmse_max <- ls_ci <- ls_ci_avg <- ls_ci_max <- vector(length = length(nVec), "list")
ls_power <- ls_time <- ls_rmse_avg
itrs <- 1000 # number of simulation iterations
cnt <- 0
rmseMat <- matrix(nc = 4, nr = itrs)
pre_process <- TRUE # only needed first time after downloading from cluster
ri <- FALSE
total_boots <- 1000
boots_div <- 10
boots <- total_boots / boots_div # graphs show time in terms of 10 boots

if(pre_process){
  
  # colnames
  colnm <- c("beta_CI_joint", "beta_CI_naive", "fLME_beta_rmse", "fLME_time",
                         "fLME_avgSignif", "fLME_Joint_avgSignif", "fLME_avgbeta_rmse", "fLME_avgbeta_CI_incl",
                         "perm_avgSignif", "fLME_bias", "LME_avgbeta_rmse", "LME_avgbeta_CI_incl",
                         "perm_avgSignif_0", "perm_avgSignif_2", "ttest_avgbeta_rmse", "ttest_avgbeta_CI_incl",
                         "perm_rmse", "perm_CI_incl",
                         "beta_CI_joint_noIntrcpt", "beta_CI_naive_noIntrcpt", "fLME_beta_rmse_noIntrcpt",
                         "beta_CI_boot", "beta_CI_boot_noIntrcpt", "fLME_Maxbeta_CI_incl_boot", "fLME_reward_CI_incl_boot", "perm_CI_incl_reward",
                         "beta_CI_boot_naive", "beta_CI_boot_noIntrcpt_naive", "pffr_reward_CI_incl", "fLME_reward_CI_incl_boot_naive",
                         "fLME_bias", "fLME_simple_bias",
                         "fLME_Boot_avgSignif", "fLME_Boot_Joint_avgSignif", "fLME_auc_signif", "LME_max_signif", "LME_auc_signif",
                         "ttst_max_signif", "ttst_auc_signif", "perm_max_signif", "perm_auc_signif",
                         "fLME_boot_auc_signif_CF", "fLME_boot_auc_CI_incl", "fLME_boot_auc_signif",
                         "fLME_orig_time", "fLME_sim_time", "fLME_boot_time", "n_trials",
                         "perm_avgbeta_rmse", "perm_avgbeta_CI_incl", "perm_avgSignif_1", "pffr_avgSignif",
                         "fLME_boot_auc_unCorrect_signif", "fLME_avgbeta_CI_incl_boot_unCorrect",
                         "pffr_CI_joint",  "pffr_beta_rmse", "pffr_time", "pffr_CI_joint_noIntrcpt","pffr_beta_rmse_noIntrcpt", "pffr_bias_noIntrcpt",
                         "pffr_avgbeta_CI_incl",  "pffr_auc_signif", "pffr_avgbeta_rmse",
                         "beta_CI_joint0", "beta_CI_boot_joint0", "pffr_CI_joint0", 
                         "fLME_reward_avgCIs",  "fLME_reward_joint_avgCIs", rep("NA", 14))
  
  # extra name for sample size corrected bootstrap CIs
  bootstrap_col_nms <- c(25, 30, 22, 27, 23, 65, 28, 33, 34, 44, 53, 42, 43, 54)
  colnm[69:82] <- paste0(colnm[bootstrap_col_nms], "_ss")

  # pre-process longer reward-period 
  for(n in nVec){
    for(rwdLen in rwdLen_vec){
      
      flNm <-  paste0("photometry_sims_science__", 
                      "randInt_", ri,
                      "_trl_trm_1fixed_smooth_", fixed_smooth,
                      "_knots_div_", knots_div,
                      "_resid_indep_", indp,
                      "_resid_subj_", resid_var_subj,
                      "_resid_scl_", resid_scl, 
                      "_fix_scl_", fix_scl,
                      "_fix_sft_", fix_shift,
                      "_rwdLen_", rwdLen,
                      "_n_trgt_", n)  
      
      filePrefix <- paste0(flNm, "_")
      
      fls <- list.files(wd) 
      fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
      
      if(length(fls_nms) > 0 )     process_files(fileNm=flNm, save.folder = wd, itrs = itrs, colnm = colnm)
      
    }
  }
  
}


for(n in nVec){
  for(rwdLen in rwdLen_vec){

    flNm <-  paste0("photometry_sims_science__", 
                    "randInt_", ri,
                    "_trl_trm_1fixed_smooth_", fixed_smooth,
                    "_knots_div_", knots_div,
                    "_resid_indep_", indp,
                    "_resid_subj_", resid_var_subj,
                    "_resid_scl_", resid_scl, 
                    "_fix_scl_", fix_scl,
                    "_fix_sft_", fix_shift,
                    "_rwdLen_", rwdLen,
                    "_n_trgt_", n)  

    if(file.exists(flNm)){
      # check to see if file exists
      
      cnt <- cnt + 1
      d <- read.csv(flNm) 
      
      # remove NAs
      ind <- apply(d, 1, function(x) all(is.na(x)))
      d <- d[!ind,]
      
      # rmse acorss entire functional domain
      rmseMat <- data.frame("FLMM" = d$fLME_beta_rmse_noIntrcpt,
                            "Perm" = d$perm_rmse,
                            "pffr" = d$pffr_beta_rmse_noIntrcpt)
      
      names(rmseMat) <- gsub(x = names(rmseMat), pattern = "\\.", replacement = "-")
      
      # joint CIs
      max_rmseMat <- data.frame("FLMM" = 1*(d$beta_CI_joint_noIntrcpt),
                                "Perm" = 1*(d$perm_CI_incl),
                                "Boot" = 1*(d$beta_CI_boot_noIntrcpt_ss),
                                "pffr" = 1*(d$pffr_CI_joint_noIntrcpt ))
      
      names(max_rmseMat) <- gsub(x = names(max_rmseMat), pattern = "\\.", replacement = "-")
      
      # RMSE of the AUC
      avg_rmseMat <- data.frame("FLMM" = d$fLME_avgbeta_rmse,
                            "LME" = d$LME_avgbeta_rmse,
                            "t-test" = d$ttest_avgbeta_rmse,
                            "perm" = d$perm_avgbeta_rmse,
                            "pffr" = d$pffr_avgbeta_rmse)
      
      names(avg_rmseMat) <- gsub(x = names(avg_rmseMat), pattern = "\\.", replacement = "-")
      
      # joint CI no intercept
      ciMat <- data.frame("FLMM" = 1*I(d$beta_CI_joint_noIntrcpt == 1),
                          "Perm" = 1*I(d$perm_CI_incl == 1),
                          "Boot" = 1*I(d$beta_CI_boot_noIntrcpt_ss == 1),
                          "pffr" = 1*I(d$pffr_CI_joint_noIntrcpt == 1))
      
      names(ciMat) <- gsub(x = names(ciMat), pattern = "\\.", replacement = "-")
      
      # average over reward period CI inclusion
      avg_ciMat <- data.frame("FLMM" = d$fLME_reward_avgCIs, # average CI inclusion of reward period
                              "fLME_auc" = d$fLME_avgbeta_CI_incl, # AUC CI inclusion (not average over reward period)
                              "LME" = d$LME_avgbeta_CI_incl,
                              "t-test" = d$ttest_avgbeta_CI_incl,
                              "Boot" = d$fLME_reward_CI_incl_boot_naive_ss, # average CI inclusion over reward period
                              "Perm" = d$perm_CI_incl_reward, # average CI inclusion over reward period
                              "Perm_auc" = d$perm_avgbeta_CI_incl, # AUC CI inclusion (not average over reward period)
                              "Boot_auc" = d$fLME_boot_auc_CI_incl_ss, # AUC CI inclusion (not average over reward period)
                              "pffr" = d$pffr_reward_CI_incl)
      
      names(avg_ciMat) <- gsub(x = names(avg_ciMat), pattern = "\\.", replacement = "-")
      
      # pointwise CIs (all no intercept so can compare with perm)
      max_ciMat<- data.frame("FLMM" = d$beta_CI_naive_noIntrcpt,
                                  "Perm" = d$perm_CI_incl,
                                  "Boot" = d$beta_CI_boot_noIntrcpt_naive_ss,
                                  "pffr" = d$pffr_CI_joint_noIntrcpt)
      
      names(max_ciMat) <- gsub(x = names(max_ciMat), pattern = "\\.", replacement = "-")
      
      power_Mat <- data.frame("fLME_avgofSig" = d$fLME_avgSignif,
                              "fLME_cf" = d$fLME_auc_signif,
                              "LME" = d$LME_auc_signif,
                              "t-test" = d$ttst_auc_signif,
                              "Boot_aos" = d$fLME_Boot_avgSignif_ss,
                              "Boot_joint_aos" = d$fLME_Boot_Joint_avgSignif_ss,
                              "boot_cf" = d$fLME_boot_auc_signif_CF_ss,
                              "Perm" = d$perm_auc_signif,
                              "Perm_aoS_half" = d$perm_avgSignif,
                              "Perm_aoS1" = d$perm_avgSignif_1,
                              "pffr" = d$pffr_auc_signif,
                              "pffr_aoS" = d$pffr_avgSignif)
      
      names(power_Mat) <- gsub(x = names(power_Mat), pattern = "\\.", replacement = "-")
      
      time_Mat <- data.frame("FLMM" = d$fLME_time,
                              "Boot" = d$fLME_boot_time / boots,
                              "pffr" = d$pffr_time)
     # "trials" = d$n_trials,
      
      names(time_Mat) <- gsub(x = names(time_Mat), pattern = "\\.", replacement = "-")
      
      ls_time[[cnt]] <- cbind( 
        gather(time_Mat), 
        n, rwdLen, d$n_trials
      )
      
      ls_rmse[[cnt]] <- cbind( 
        gather(rmseMat), 
        n, rwdLen
      )
      
      ls_rmse_max[[cnt]] <- cbind( 
        gather(max_rmseMat), 
        n, rwdLen
      )
      
      ls_rmse_avg[[cnt]] <- cbind( 
        gather(avg_rmseMat), 
        n, rwdLen
      )
      
      ls_ci[[cnt]] <- cbind( 
        gather(ciMat), 
        n, rwdLen
      )
      
      ls_ci_avg[[cnt]] <- cbind( 
        gather(avg_ciMat), 
        n, rwdLen
      )
      
      ls_ci_max[[cnt]] <- cbind( 
        gather(max_ciMat), 
        n, rwdLen
      )
      
      ls_power[[cnt]] <- cbind( 
        gather(power_Mat), 
        n, rwdLen
      )
      
      d1 <- d
      rm(d)
    }
    
    
  }
}


###########################################
# set colors
levs <- c("FLMM", "LME", "t-test", "Perm", "Perm 1/2", "Boot")
myColors <- setNames( c("#ca0020", "#0868ac", "darkgray", "#E69F00", "#525252", "lightgrey"), levs) # , 

#########################
# rmse of average
#########################
###########################################
# factors
dat <- do.call(rbind, ls_rmse_avg)
dat$n <- as.factor(dat$n)
dat$key <- ifelse(dat$key == "perm", "Perm", dat$key)
###########################################
plt_rmse = 
  dat %>% tibble %>%  
  dplyr::filter(key %in% levs,
                # n %in% c(4,6,8,10)) %>%
                n %in% seq(4,8)) %>%
  ggplot(aes( y = value, x = n, fill = key )) +
  geom_boxplot(fatten = 3) + 
  ylab("Cue: \u03B2 RMSE" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "Perm 1/2", "Boot")] ) +
  scale_color_manual(values = myColors[!levs %in% c( "Perm 1/2", "Boot")] ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Method"))  

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_rmse_avg.pdf",
        plot = plt_rmse,
        width = 12,
        height = 6)

##################################################
# Joint CIs across functional domain
##################################################
###########################################
dat <- do.call(rbind, ls_ci)
grp_incl <- c("FLMM", "Perm", "Boot") # include these methods
###########################################
joint_ci = 
  dat %>% tibble %>%  
  dplyr::filter(value <= 1) %>% # some erroneous large values
  dplyr::filter(key %in% grp_incl) %>%
  dplyr::group_by(key, n, rwdLen) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes( y = value, x = as.numeric(n), fill = key )) +
  geom_line( lwd = 1.5,
             aes(colour = key) ) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1), # 0.5
             alpha = 0.7) + #
  ylab("Joint 95% CI Coverage" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[levs %in% grp_incl] ) +
  scale_color_manual(values = myColors[levs %in% grp_incl] ) +
  theme_classic(base_size = 12) +
  scale_x_continuous(labels=seq(4, max(nVec), by = 2), 
                     breaks=seq(4, max(nVec), by = 2),
                     limits = c( min(nVec), max(nVec) )) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(colour= guide_legend(title="Method"))  

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_Joint_CIs.pdf",
        plot = joint_ci,
        width = 12,
        height = 6)
###########################################


######################################################################################
# AUC Pointwise CIs -- Confidence intervals for the AUC measures
# Perm and fLME are average of CI inclusions over reward period (not AUCs)
######################################################################################
dat <- do.call(rbind, ls_ci_avg)
grp_incl <- c("FLMM", "LME", "Perm", "t-test", "Boot") # include these methods
###########################################

# make binary data into line plots
avg_dat <- dat %>% tibble %>%  
  dplyr::filter(key %in% c( "LME", "t-test"),
                n %in% seq(4,8),
                value <= 1) %>%
  dplyr::group_by(key, n, rwdLen) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE))

# continuous data
box_dat <- dat %>% tibble %>%  
  dplyr::group_by(key, n) %>%
  dplyr::filter(key %in% levs[!levs %in% c( "LME", "t-test")],
                n %in% seq(4,8) )

plt_box = 
  ggplot() +
  # geom_boxplot(data = box_dat, 
  #              aes(y = value, x = as.factor(n), fill = key), 
  #              color = "black", fatten = 3) + 
  geom_line(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
            lwd = 1.5) + 
  geom_point(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
             lwd = 1.5) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + #
  ylab("Cue: Pointwise 95% CI Coverage" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test", "Perm 1/2")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Pointwise"),
         color = guide_legend(title = "Proportion")) 

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_CIs_of_avg.pdf",
        plot = plt_box,
        width = 12,
        height = 6)
###########################################

##################################################
# power cue-period
##################################################
###########################################
dat <- do.call(rbind, ls_power)
dat$key <- ifelse(dat$key == "Perm_aoS_half", "Perm 1/2", dat$key)
dat$key <- ifelse(dat$key == "Perm_aoS1", "Perm", dat$key)
dat$key <- ifelse(dat$key == "fLME_avgofSig", "FLMM", dat$key)
dat$key <- ifelse(dat$key == "Boot_aos", "Boot", dat$key)
grp_incl <- c("FLMM", "LME", "Perm", "t-test", "Boot") # include these methods
###########################################

# make binary data into line plots
avg_dat <- dat %>% tibble %>%  
              dplyr::filter(key %in% c( "LME", "t-test"),
                            n %in% seq(4,8),
                            value <= 1) %>%
              dplyr::group_by(key, n, rwdLen) %>%
              dplyr::summarise(value = mean(value, na.rm = TRUE))

# continuous data
box_dat <- dat %>% tibble %>%  
  dplyr::group_by(key, n) %>%
  dplyr::filter(key %in% levs[!levs %in% c( "LME", "t-test")],
                n %in% seq(4,8))
  
plt_box = 
  ggplot() +
  # geom_boxplot(data = box_dat, 
  #              aes(y = value, x = as.factor(n), fill = key), 
  #              color = "black", fatten = 3) + 
  geom_line(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
            lwd = 1.5) + 
  geom_point(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
            lwd = 1.5) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + 
  ylab("Cue: Proportion Significant" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Signif"),
             color = guide_legend(title = "Prop Signif")) 

plt_violin = 
  ggplot() +
  facet_wrap( ~ rwdLen, nrow = 1) +
  geom_violin(data = box_dat, 
              aes(y = value, x = as.factor(n), fill = key), 
              color = "black") + 
  geom_line(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
            lwd = 1.5) + 
  geom_point(data = avg_dat, aes(y = value, x = as.factor(n), colour = key, group = key),
             lwd = 1.5) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + 
  ylab("Proportion Significant" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Signif"),
         color = guide_legend(title = "Prop Signif")) 



library(patchwork)
set.seed(25)
plot_combined <- plt_box + plt_violin +
  plot_layout(ncol = 1,
              nrow = 2,
              byrow = NULL,
              widths = NULL,
              heights = 10,
              guides = NULL)

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_avg_power.pdf",
        plot = plt_box,
        width = 12,
        height = 6)
###########################################


###########################################
# time
###########################################
# factors
dat <- do.call(rbind, ls_time)
dat$n <- as.factor(dat$n)
dat$key <- ifelse(dat$key == "perm", "Perm", dat$key)
grp_incl <- c("FLMM") 

# covnert sample sizes into trials (averaged)
dd <- dat %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(trials = as.factor(signif( mean(`d$n_trials`, na.rm = TRUE), digits = 0) ))
dat <- left_join(dat, dd, by = "n")
dat %>% group_by(n) %>% dplyr::filter(key == "FLMM") %>% summarise(cd=n())
dat %>% group_by(n) %>% dplyr::filter(key == "FLMM") %>% summarise(cd=mean(value))

#################################################################################################
#----------------------------
# comparisons with pffr
#----------------------------

###########################################
# set colors
levs <- c("FLMM", "pffr") #"LME", "t-test", "Perm", "Perm 1/2", "Boot")
myColors <- setNames( c("#ca0020", "#0868ac", "darkgray", "#E69F00", "#525252", "lightgrey"), levs) # , 

#########################
# rmse of average
#########################
###########################################
# factors
dat <- do.call(rbind, ls_rmse_avg)
dat$n <- as.factor(dat$n)
dat$key <- ifelse(dat$key == "perm", "Perm", dat$key)
###########################################
plt_rmse = 
  dat %>% tibble %>%  
  dplyr::filter(key %in% levs,
                # n %in% c(4,6,8,10)) %>%
                n %in% seq(4,8)) %>%
  ggplot(aes( y = value, x = n, fill = key )) +
  geom_boxplot(fatten = 3) + 
  ylab("Cue: Fixed-Effects RMSE") + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "Perm 1/2", "Boot")] ) +
  scale_color_manual(values = myColors[!levs %in% c( "Perm 1/2", "Boot")] ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Method"))  

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_rmse_avg_pffr.pdf",
        plot = plt_rmse,
        width = 12,
        height = 6)

##################################################
# Joint CIs across functional domain
##################################################
###########################################
dat <- do.call(rbind, ls_ci)
grp_incl <- c("FLMM", "pffr") # include these methods
###########################################
joint_ci = 
  dat %>% tibble %>%  
  dplyr::filter(value <= 1) %>% # some erroneous large values
  dplyr::filter(key %in% grp_incl) %>%
  dplyr::group_by(key, n, rwdLen) %>%
  dplyr::summarise(value = mean(value, na.rm = TRUE)) %>%
  ggplot(aes( y = value, x = as.numeric(n), fill = key )) +
  geom_line( lwd = 1.5,
             aes(colour = key) ) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1), # 0.5
             alpha = 0.7) + #
  ylab("Joint 95% CI Coverage" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[levs %in% grp_incl] ) +
  scale_color_manual(values = myColors[levs %in% grp_incl] ) +
  theme_classic(base_size = 12) +
  scale_x_continuous(labels=seq(4, max(nVec), by = 1), 
                     breaks=seq(4, max(nVec), by = 1),
                     limits = c( min(nVec), max(nVec) )) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 3 GCL
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(colour= guide_legend(title="Method"))  

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_Joint_CIs_pffr.pdf",
        plot = joint_ci,
        width = 12,
        height = 6)
###########################################


######################################################################################
# AUC Pointwise CIs -- Confidence intervals for the AUC measures
# Perm and fLME are average of CI inclusions over reward period (not AUCs)
######################################################################################
dat <- do.call(rbind, ls_ci_avg)
grp_incl <- c("FLMM", "pffr") # include these methods
###########################################

# continuous data
box_dat <- dat %>% tibble %>%  
  dplyr::group_by(key, n) %>%
  dplyr::filter(key %in% levs[!levs %in% c( "LME", "t-test")],
                n %in% seq(4,8) )

plt_box = 
  ggplot() +
  geom_boxplot(data = box_dat, 
               aes(y = value, x = as.factor(n), fill = key), 
               color = "black", fatten = 3) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + #
  ylab("Cue: Pointwise 95% CI Coverage" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test", "Perm 1/2")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Pointwise"),
         color = guide_legend(title = "Proportion")) 

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_CIs_of_avg_pffr.pdf",
        plot = plt_box,
        width = 12,
        height = 6)
###########################################

##################################################
# power cue-period
##################################################
###########################################
dat <- do.call(rbind, ls_power)
dat$key <- ifelse(dat$key == "fLME_avgofSig", "FLMM", dat$key)
dat$key <- ifelse(dat$key == "pffr_aoS", "pffr", dat$key)
grp_incl <- c("FLMM", "pffr") # include these methods
###########################################

# continuous data
box_dat <- dat %>% tibble %>%  
  dplyr::group_by(key, n) %>%
  dplyr::filter(key %in% levs[!levs %in% c( "LME", "t-test")],
                n %in% seq(4,8))

plt_box = 
  ggplot() +
  geom_boxplot(data = box_dat, 
               aes(y = value, x = as.factor(n), fill = key), 
               color = "black", fatten = 3) + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + 
  ylab("Cue: Proportion Significant" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Signif"),
         color = guide_legend(title = "Prop Signif")) 

plt_violin = 
  ggplot() +
  facet_wrap( ~ rwdLen, nrow = 1) +
  geom_violin(data = box_dat, 
              aes(y = value, x = as.factor(n), fill = key), 
              color = "black") + 
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(1),
             alpha = 0.7) + 
  ylab("Proportion Significant" ) + 
  xlab( "Sample Size" ) + 
  scale_fill_manual(values = myColors[!levs %in% c( "LME", "t-test")] ) +
  scale_color_manual(values = myColors[levs %in% c( "LME", "t-test")] ) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(0, 1) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), 
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Avg Signif"),
         color = guide_legend(title = "Prop Signif")) 



library(patchwork)
set.seed(25)
plot_combined <- plt_box + plt_violin +
  plot_layout(ncol = 1,
              nrow = 2,
              byrow = NULL,
              widths = NULL,
              heights = 10,
              guides = NULL)

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_avg_power_pffr.pdf",
        plot = plt_box,
        width = 12,
        height = 6)
###########################################


###########################################
# time
###########################################
# factors
dat <- do.call(rbind, ls_time)
dat$n <- as.factor(dat$n)
dat$key <- ifelse(dat$key == "perm", "Perm", dat$key)
grp_incl <- c("FLMM", "pffr") 

# covnert sample sizes into trials (averaged)
dd <- dat %>%
  dplyr::group_by(n) %>%
  dplyr::summarise(trials = as.factor(signif( mean(`d$n_trials`, na.rm = TRUE), digits = 0) ))
dat <- left_join(dat, dd, by = "n")
dat %>% group_by(n) %>% dplyr::filter(key == "FLMM") %>% summarise(cd=n())
dat %>% group_by(n) %>% dplyr::filter(key == "FLMM") %>% summarise(cd=mean(value))

###########################################

plt_time = 
  dat %>% tibble %>%  
  dplyr::filter(key %in% grp_incl,
                n %in% c(4:8)) %>%
  ggplot(aes( y = value * 60, x = trials)) +
  geom_boxplot(aes(x = trials, fill = key)) +
  ylab("Time to fit (sec)" ) + 
  xlab( "Total Trials in Dataset" ) + 
  scale_fill_manual(values = myColors[levs %in% grp_incl] ) +
  scale_color_manual(values = myColors[levs %in% grp_incl] ) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(2, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(fill= guide_legend(title="Method"))  

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_time.pdf",
        plot = plt_time,
        width = 12,
        height = 6)
