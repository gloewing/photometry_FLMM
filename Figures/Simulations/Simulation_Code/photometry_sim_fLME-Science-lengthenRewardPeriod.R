# written by Gabe Loewinger 4/4/23
# simulations using Nature paper's data for data-driven sims
# this R code is saved in: /gpfs/gsfs8/users/loewingergc/photometry_fglmm/code
# reward period lengthened
list.of.packages <- c("Rfast")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  chooseCRANmirror(ind=75)
  install.packages(new.packages, dependencies = TRUE)
} 

library(lme4) ## mixed models
library(refund) ## fpca.face 
library(dplyr) ## organize lapply results
library(progress) ## display progress bar
library(mgcv) ## smoothing in step 2
library(mvtnorm) ## joint CI
library(parallel) ## mcapply
library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(data.table)
library(tidyverse)
library(Rfast)

# params matrix
# simParams <- expand.grid(TRUE, c(FALSE), FALSE, c(4:8), c(1) )
# colnames(simParams) <- c("fix_smooth", "resid_var_indep", "resid_var_subj", "n", "beta_mult")
# write.csv( simParams, "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final/Code/Utility Functions/photo_params.csv", row.names = FALSE)

cluserInd <- TRUE # whether running on computer or on cluster
num_cores <- 10 # parallelize

# save and load
if(cluserInd){
  # if on cluster
  args = commandArgs(TRUE)
  runNum <- as.integer( as.numeric(args[1]) )
  iter <- as.integer( Sys.getenv('SLURM_ARRAY_TASK_ID') ) # seed index from array id
   
  # this R code is saved in: /gpfs/gsfs8/users/loewingergc/photometry_fglmm/code
  # bash in: /gpfs/gsfs8/users/loewingergc/bash
   
  # paths
  wd <- "/gpfs/gsfs8/users/loewingergc/code_utils/"
  data_path <-"/gpfs/gsfs8/users/loewingergc/photometry_fglmm/data" # path to original Matlab files for data
  save.folder <- "/gpfs/gsfs8/users/loewingergc/photometry_fglmm/sims_delays" # photometry sims
  source(paste0(wd, "photometry_sim_fLME_fn_multi.R")) # simulations
  source(paste0(wd, "fui.R")) # function code
  simParams <- read.csv(paste0(data_path,"/photo_params.csv") )
   
}else{
  runNum <- 5
  wd <- "~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/"
  iter <- 9
  data_path <-"/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper" # path to original Matlab files for data
  source("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photometry_sim_fLME_fn_multi.R") # simulations
  source(paste0(wd, "fui.R")) # function code
  source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM.R')
  source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
  
  simParams <- expand.grid(TRUE, c(FALSE), FALSE, c(4:10), c(1) )
  colnames(simParams) <- c("fix_smooth", "resid_var_indep", "resid_var_subj", "n", "beta_mult")
}
  
# simulation parameters
fixed_smooth <- simParams[runNum, 1] #TRUE # whether to use smoothed coefficients as truth or not
target_Hz <- 15
reps <- 1000 # replicates of simulations
boots <- 5000 # for naive bootstrap permutation test
knots_div <- 2 # number of functional samples divisor (so if it is 1, there is a knot for every obsveration on fn domain)
resid_var_indep <- simParams[runNum, 2] # FALSE # assume residual variance is indepedent across fn. domain as in paper
resid_var_subj <- TRUE # simParams[runNum, 3] #TRUE # residual variance subject specific or global
n_sim <- simParams[runNum, 4] 
fixed_effects_scale <- 1 #simParams[runNum, 5]
random_effects_scale <- 1
resid_scale <- 5 # 
comparison_methods <- TRUE # whether to fit other models besides fLME
trial_trim <- 1 # take every other trial to reduce computational burden
rand_int <- FALSE
fixed_shift <- 0
flme_boot <- TRUE  
cue_period_length <- 2 # length of period of cue period for 
 
fileNm <- paste0("photometry_sims_science_",
                 "_randInt_", rand_int,
                 "_trl_trm_", trial_trim,
                   "fixed_smooth_", fixed_smooth,
                   "_knots_div_", knots_div,
                   "_resid_indep_", resid_var_indep,
                   "_resid_subj_", resid_var_subj,
                 "_resid_scl_", resid_scale,
                  "_fix_scl_", fixed_effects_scale,
                 "_fix_sft_",fixed_shift,
                 "_rwdLen_", cue_period_length,
                   "_n_trgt_", n_sim)

# results matrix
res_mat <- matrix(NA, ncol = 82, nrow = reps)
colnames(res_mat) <- c("beta_CI_joint", "beta_CI_naive", "fLME_beta_rmse", "fLME_time",
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
colnames(res_mat)[69:82] <- paste0(colnames(res_mat)[bootstrap_col_nms], "_ss")

#################################
# original data variables
#################################

# day when cue duration changes
mouseList <- c('HJ-FP-M2', 'HJ-FP-M3', 'HJ-FP-M4', 'HJ-FP-F1', 'HJ-FP-F2', 'HJ-FP-M6', 'HJ-FP-M7') # c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-M6','HJ-FP-M7', 'HJ-FP-F1','HJ-FP-F2')
session_max <- c(32, 27, 37, 22, 27, 22, 23) # last day of each condition
cue_change <- c(29, 24, 32, 19, 24, 19, 20) # c(29, 24, 32, 19, 20, 19, 24)
csPlus_vec <- c(15, 15, 16, 16, 15, 15, 16) # c(15, 15, 16, 15, 16, 16, 15)
cue_mat <- as.data.frame(cbind(mouseList, cue_change, csPlus_vec))
colnames(cue_mat) <- c("id", "day", "cs")
# cbind(dat$eventlog$nonsolenoidflag, dat$eventlog$eventtime, dat$eventlog$eventindex)

n <- length(mouseList)
preSessions <- 3 # number of days before change day
session_min <- cue_change - (preSessions + 1) # assumes this folder only has files we need
pre_min_tm <- 1 # pre_lick period
post_min_tm <- 5 # post lick period
aic_period_length <- 5 # used for AIC
trial_num <- 100 # trials per session
reward_period <- seq(pre_min_tm* target_Hz+1, cue_period_length*target_Hz)  # reward period used in paper

# samples before after
Hz <- 1 / 0.008
pre_samps <- round(pre_min_tm * Hz) # samples before lick
post_samps <- round(post_min_tm * Hz) # samples after lick

#################################################
if( !cluserInd ){
  
  path <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_4/data/"
  
  # event indices
  index_lick = 5
  index_bgdrw = 10 # reward
  index_s1 <- 15 # cue1
  index_s2 <- 16 # cue 2
  index_sessionend = 0
  index_trial_start <- 12 # for both cs+ and cs-
  
  dat_lst <- vector(length = 5 * n, "list")
  cnt <- 0
  
  for(i in 1:n){
    
    sess <- as.integer(cue_mat$day[i]) # change of cue
    sess_vec <- seq(session_min[i], session_max[i])
    id <- cue_mat$id[i] # id name
    index_csPlus <- as.integer(cue_mat$cs[i] ) # index of cs+
    index_csMinus <- ifelse(index_csPlus == index_s1, index_s2, index_s1)
    cnt_i <- 0
    for(s in sess_vec){
      suff <- ifelse(s == sess, "-8s.nwb", ".nwb")  
      
      # file name
      fl_nm <- paste0("sub-", id, "_ses-Day", s, suff)
      
      #  check if file exists
      fl_list <- list.files(path = path, pattern = fl_nm)
      
      # if exists
      if(length(fl_list) == 1){ 
        
        # read ttl data
        dat <- rhdf5::h5read(file = paste0(path, fl_nm),
                             name = "/acquisition")
        
        # event time ttls
        ttls <- as.data.frame( cbind(dat$eventlog$eventtime, dat$eventlog$eventindex) )
        colnames(ttls) <- c("time", "event")
        
        # read photometry use df/f
        dat_photo <- rhdf5::h5read(file = paste0(path, fl_nm),
                                   name = "/processing/photometry/dff")
        
        dat_photo <- as.data.frame( do.call(cbind, dat_photo) )
        dat_photo <- dat_photo[rowSums(is.na(dat_photo)) < nrow(dat_photo), ] # remove NAs
        
        #########################
        # data 
        #########################
        num_trials <- length(ttls$time[ttls$event == index_trial_start])
        
        cnt <- cnt + 1 # counter for list
        cnt_i <- cnt_i + 1
        
        # align times of behavior to photometry times
        ttls <- time_align(time_truth = dat_photo$timestamps, data = ttls, name = "time") #align times in ttls to the times in photometry data
        
        # find cue length for each CS+
        ttls_cs <- ttls %>% as_tibble() %>% dplyr::filter(event == index_csPlus ) # CS+ or reward
        ttls_rew <- ttls %>% as_tibble() %>% dplyr::filter(event == index_bgdrw ) # CS+ or reward
        ttls_csMinus <- ttls %>% as_tibble() %>% dplyr::filter(event == index_csMinus ) # CS-
        
        # save "reward" times that are closest to cs+
        rew_idx_ttl <- lapply(ttls_cs$time, function(x) which( (ttls_rew$time - x) > 0  )[1] )
        rew_idx_ttl <- do.call(c, rew_idx_ttl )
        ttls_rew <- ttls_rew[unique(rew_idx_ttl), ]
        
        # ensure number of rewards == number of CS+
        if(length(ttls_rew$time) != length(ttls_cs$time)){
          print(paste("i", i, "cnt", cnt_i))
          break
        }   
        reward_lengths <- round(ttls_rew$time - ttls_cs$time, 2) # find time between CS+ and reward
        
        # indices of CS+/CS-
        idx_ttls_cs <- sapply(ttls_cs$time, function(x) which(dat_photo$timestamps == x))  #which(ttls$event == index_csPlus ) # CS+ or reward
        idx_ttls_csMinus <- sapply(ttls_csMinus$time, function(x) which(dat_photo$timestamps == x)) # CS+ or reward
        
        #
        idx_mat <- data.frame( idx = c(idx_ttls_cs, idx_ttls_csMinus), # indices
                               onset = c(ttls_cs$time, ttls_csMinus$time), # time when cue starts
                               reward_delay = c(reward_lengths, rep(NA, length(idx_ttls_csMinus) )), # delay to reward
                               csPlus = c(rep(1, length(reward_lengths)), rep(0, length(idx_ttls_csMinus) )) ) # cs+/cs-
        rm(ttls_rew)
        
        # center around start of CS+
        idx <- t( sapply(idx_mat$idx, function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for trials based on CS start times
        dat <- as.data.frame( apply(idx, 2, function(x) dat_photo$data[x]) ) # photometry trials
        L <- ncol(dat)
        
        # downsample photometry
        by_num <- round( Hz / target_Hz )
        photo_idx_pre <- sort(-seq(-pre_samps, -1, by = by_num)) # seq(1, pre_samps, by = by_num)
        photo_idx_post <- seq(pre_samps + by_num, L, by = by_num)
        photo_idx <- unique( c(photo_idx_pre, photo_idx_post) )
        dat <- dat[,photo_idx]
        L <- ncol(dat)
        colnames(dat) <- paste0("photometry.", 1:L)
        
        # combine data
        trials <- seq(1, num_trials )
        sess_info <- as.data.frame( cbind(id, cnt_i - (preSessions + 1), c(1:length(reward_lengths), rep(NA, length(ttls_csMinus$time))), 
                                          idx_mat[,3:4] ) )
        colnames(sess_info) <- c("id", "session", "trial", "delay", "cs")
        
        dat_lst[[cnt]] <- as.data.frame( cbind(sess_info, dat) )
        
        rm(dat, trials, sess_info, ttls, dat_photo, idx)
        
      }
    }
  }
  
  # cocnatenate results
  dat <- do.call(rbind, dat_lst)
  
  #####################################################################################
  # consider centering iri, trial, etc, so intercept is interpretable and more likely to be significant
  
  # photometry sample length
  out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
  L <- length(out_index)
  nknots_min <- round(L/4)
  cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates
  preSamps <- pre_min_tm * target_Hz
  reward_period_idx <- seq(preSamps + 1, aic_period_length * target_Hz + preSamps) # for AIC
  
  # remove NAs
  dat <- dat[complete.cases(dat[,out_index]),] # remove trials with missing photometry data
  dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
  dat2 <- dat
  rm(dat_lst)
  #########################################################################################
  # creat AUC variables for model comparison
  # Supplement:
  # Experiments 2-5: To analyze learning-dependent changes of behavior, the cumulative sum of anticipatory lick rate across
  # trials was calculated. Lick rate during the baseline period (1 s window before CS+ onset) was subtracted from lick rate
  # during 2 s from CS+ onset. Similarly, to analyze learning-dependent dynamics of dopamine response, AUC of ΔF/F during
  # 2 s from CS+ onset was normalized by AUC during baseline period.
  # ------ # used the above to calculate "dat_AUC_cue" below
  #-------------------------
  # To test if the omission response appears when cue duration changes from 2 s to 8 s (Experiment 3), in 8 s cue trials, we
  # calculated baseline subtracted AUC of ΔF/F for 3 to 4 s from cue onset, which was the reward period before the cue duration
  # change. This was compared to the omission response in extinction (Experiment 4), in which we used the baseline subtracted
  # AUC of ΔF/F for 1 s after the offset of trace interval (i.e., the original reward time).
  photo_dat <- dat[, -cov_idx] # photometry data
  
  # preSamps is the 2 seconds before cue onset
  cue_period <- seq(preSamps + 1, preSamps + 2 * target_Hz)
  cue_post <- seq(max(cue_period) + 1, preSamps + 5 * target_Hz)
  
  dat_AUC_pre <- rowMeans( photo_dat[, 1:preSamps] ) # AUC before reward period
  dat_AUC_cue <- rowMeans(photo_dat[, cue_period]) # AUC during 2 second post CS+ onset
  dat_AUC_post <- rowMeans(photo_dat[, cue_post]) # AUC  -- ΔF/F for 3 to 4 s from cue onset (from above)
  # 
  # # scale covariates to make contrasts more interpretable
  dat$dat_AUC_pre <- dat_AUC_pre
  dat$dat_AUC_post <- dat_AUC_post - dat_AUC_pre
  dat$dat_AUC_cue <- dat_AUC_cue - dat_AUC_pre
  
  dat$delay <- ifelse(dat$delay > 5, 1, 0) # 1 means it is long
  #########################################################################################
  
  dat3 <- dat %>%
    as_tibble() %>%
    dplyr::filter( #session == 0 | session == -1 & trial %in% 31:50,
      session %in% c(-1,0), # day 0 is when behavior changes
      cs == 1) %>% # CS+ trials only
    dplyr::rename(cue = delay) %>%
    as.data.frame()
  
  # write.csv(dat3,
  #           "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper/science_delay_sim_data.csv",
  #           row.names = FALSE)
  # dat <- read.csv("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper/science_delay_sim_data.csv")

  rm(dat)
  
}else{
  
  dat <- read.csv(paste0(data_path, "/science_delay_sim_data.csv") )# "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper/science_delay_sim_data.csv") #%>%
  
}

# trial trim to reduce computational burden
trls_reduced <- unique(dat$trial)[as.integer(round( seq(1, length(unique(dat$trial)), by = trial_trim)   ))]
dat <- dat[dat$trial %in% trls_reduced, ] # reduce number of trials analzed
dat$id <- as.numeric(as.factor(dat$id)) # make 1:n

# number of knots
L <- ncol(dat[, (grep(paste0("^", "photometry"), colnames(dat)))])
nknots_min <- round(L/knots_div) # L = 100
nknots_min_cov <- round(L/4)
# order data
dat_photo <- dat[order(dat$id, dat$session, dat$trial, decreasing = FALSE), ] #%>%
n <- length(unique(dat_photo$id)) # number of subjects in (potentially) reduced sample

dat_photo$trial <- scale(dat_photo$trial) # important for stability of model

rm(dat)

############################################
# fit original model for ground truth
############################################
# fit model
set.seed(1)
timeStart <- Sys.time()

fit <- fui(formula = photometry ~ cue + 
               (cue | id), 
             data = dat_photo, 
             family = "gaussian", 
             var = TRUE, 
             parallel = FALSE, 
             silent = FALSE,
             analytic = TRUE,
             nknots_min = nknots_min,
             nknots_min_cov = nknots_min_cov,
             smooth_method = "GCV.Cp",
             splines = "tp",
             design_mat = TRUE,
             G_return = TRUE,
             residuals = TRUE,
             subj_ID = "id",
             num_cores = num_cores)

timeEnd <- Sys.time()

res_mat[iter,45] <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))

# number of covariates
p <- nrow(fit$betaHat)

##################
# simulate data
##################
set.seed(iter)
beta_idx <- 2 # this corresponds to cue effect --only evaluate performance of models on these to make comparable with t-tests

dat_sim <- photo_stimulate(X = NULL, # design matrix
                           Z = NULL,
                           include = dat_photo["session"],
                           N = n_sim,
                           model = fit, # fLME model object
                           fixed_smooth = fixed_smooth, # whether to use smoothed betas for true fixed effects
                           IDs = dat_photo$id,  # col of IDs (N x 1 vector)
                           resid_var_indep = resid_var_indep, # assume residual variance is indepedent across fn. domain as in paper
                           resid_var_subj = resid_var_subj, # calculate error variance in subject-specific fashion as opposed to globally
                           resid_var_scale = resid_scale,
                           fixed_effects_draw = FALSE,
                           fixed_effects_scale = fixed_effects_scale,
                           random_effects_scale = random_effects_scale,
                           rand_int = rand_int,
                           rand_int_sigma = FALSE, #rand_int_sigma,
                           ID_level_re = TRUE, # indicates that random effects are subject-level only (do not vary within subject)
                           fixed_shift = fixed_shift, # amount to shift fixed effect functional coefficient
                           fixed_shift_idx = 2 )

#--------------------------------------------------------------------------------------------------------------------
#########################################################################
# analyze data with various approaches and assess 95% CI coverage 
#########################################################################

#----------------------------------------------------------
# A) functional LME
#----------------------------------------------------------
rm(fit)
##############################
# fit model to simulated data
##############################
# fit model and get model CIs
timeStart <- Sys.time()
fit <- fui(photometry ~ cue +
               (cue | ID), 
             data = dat_sim$data, 
             family = "gaussian", 
             var = TRUE, 
             parallel = FALSE, 
             silent = FALSE,
             analytic = TRUE,
             nknots_min = nknots_min,
             nknots_min_cov = nknots_min_cov,
             smooth_method = "GCV.Cp",
             splines = "tp",
             design_mat = TRUE,
             subj_ID = "ID",
             num_cores = num_cores)

timeEnd <- Sys.time()

res_mat[iter,4] <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))

##############################
# generate 95% CIs
##############################
lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = nrow(fit$betaHat), ncol = ncol(fit$betaHat)) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # joint CIs
  lower.joint[r,] = fit$betaHat[r,] - fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  upper.joint[r,] = fit$betaHat[r,] + fit$qn[r]*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
  
  # naive CIs
  lower[r,] = fit$betaHat[r,] - 1.96*sqrt(diag(fit$betaHat.var[,,r]))
  upper[r,] = fit$betaHat[r,] + 1.96*sqrt(diag(fit$betaHat.var[,,r]))
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

# average of CI inclusion over reward period
res_mat[iter, 67] <- mean( naive_incl[beta_idx, reward_period] )  # whether pointwise CI contains estimate
res_mat[iter, 68] <- mean( joint_incl[beta_idx, reward_period] )  # whether joint CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,1] <- mean(joint_incl)
res_mat[iter,2] <- mean(naive_incl)
res_mat[iter,3] <- sqrt( mean( (dat_sim$beta - fit$betaHat)^2 ) ) # rmse of model fit

# evaluate performance w/o intercept (so only cue effect here since one covariate)
res_mat[iter,19] <- mean(joint_incl[-1])
res_mat[iter,64] <- mean(joint_incl[1])
res_mat[iter,20] <- mean(naive_incl[-1])
res_mat[iter,21] <- sqrt( mean( (dat_sim$beta[-1,] - fit$betaHat[-1,])^2 ) ) # rmse of model fit
res_mat[iter,10] <- mean(dat_sim$beta[beta_idx,] - fit$betaHat[beta_idx,]) # bias of cue term

# average of significance
res_mat[iter, 5] <- mean( 1 * I(upper[beta_idx, reward_period] * lower[beta_idx, reward_period] > 0) )
res_mat[iter, 6] <- mean( 1 * I(upper.joint[beta_idx, reward_period] * lower.joint[beta_idx, reward_period] > 0) )

# sample size of sims (i.e., number of trials analyzed)
res_mat[iter,48] <- nrow(dat_sim$data)
#--------------------------------------------------------------------------------------------------------------------

#----------------------------------------------------------
# B) summary statistics for cue analysis with linear models
#----------------------------------------------------------
if(comparison_methods){
  
  photo_indices <- grep(paste0("^", "photometry"), colnames(dat_sim$data))
  
  Y <- dat_sim$data[, photo_indices] # photometry data
  
  # find summary stats in simulated data
  peak_ampl <- apply(Y[, reward_period], 1, max) # peak amplitude starting at cue onset
  peak_time <- apply(Y[, reward_period], 1, which.max) / target_Hz # time point of amplitude starting at cue onset
  avg_signal <- rowMeans(Y[, reward_period]) # average signal during cue
  
  # concatenate data
  dat_ss <- cbind(dat_sim$data[c("ID", "cue")], peak_ampl, peak_time, avg_signal)
  
  #############################
  # LME for summary statistics
  #############################
  # lme for average amplitude (use nlme for p-values): estimates identical with lmer (lme4) and lme(nlme))
  mod_auc <- nlme::lme(avg_signal ~ cue,
                       random = ~cue | ID, 
                       data = dat_ss,
                       control = lmeControl(maxIter = 5000,
                                            msMaxIter = 5000,
                                            niterEM = 5000) )
  
  ###########################################################
  # paired-samples t test of average summary statistics
  ###########################################################
  # average summary stats
  
  # average within animal across trials
  dat_ss_anova <- dat_ss %>% as_tibble() %>%
    group_by(ID, cue) %>%
    summarise(peak_ampl_mean = mean(peak_ampl),
              avg_ampl_mean = mean(avg_signal)) %>%
    as.data.frame()

  t_test_auc <- t.test(formula = avg_ampl_mean ~ cue,
                       data = dat_ss_anova,
                       alternative = "two.sided",
                       mu = 0, 
                       paired = TRUE,   
                       conf.level = 0.95)
  
  ##################################
  # summary stats comparison
  ##################################
  # *** functional LME ***
  # RMSE for fLME mean difference
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  res_mat[iter, 7] <- sqrt( mean(  (beta_avg - mean(fit$betaHat[beta_idx, reward_period]) )^2) )
  
  # variance of average
  len <- rep(1/length(reward_period), length(reward_period))
  V <- fit$betaHat.var[reward_period,reward_period,beta_idx]
  V <- eigenval_trim(V)   # trim eigenvalues to ensure matrix is PSD
  avg_var <- t(len) %*% V %*% len 
  
  # significance of auc
  avg_up <- mean(fit$betaHat[beta_idx, reward_period]) + 1.96 * sqrt(avg_var)
  avg_low <- mean(fit$betaHat[beta_idx, reward_period]) - 1.96 * sqrt(avg_var)
  
  # closed form expression
  res_mat[iter, 35] <- 1 * I(avg_up * avg_low > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
  res_mat[iter, 8] <- mean( I( beta_avg <= avg_up) & beta_avg >= mean(avg_low) )  # whether CI contains estimate
  
  # RMSE for non-functional LME mean difference
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  res_mat[iter, 11] <- sqrt( mean( (beta_avg - mod_auc$coefficients$fixed[beta_idx] )^2) )
  peak_upper <- mod_auc$coefficients$fixed[beta_idx] + 1.96 * summary(mod_auc)$tTable[,2]
  peak_lower <- mod_auc$coefficients$fixed[beta_idx] - 1.96 * summary(mod_auc)$tTable[,2]
  res_mat[iter, 12] <- mean( I( beta_avg <= peak_upper & beta_avg >= peak_lower) ) # whether CI contains estimate
  res_mat[iter, 37] <- I(summary(mod_auc)$tTable[beta_idx, 5] <= 0.05)# p-value
  rm(beta_avg, V)
  
  # mean difference t-test
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  res_mat[iter, 15] <- sqrt( mean( (beta_avg - t_test_auc$estimate )^2) )
  peak_upper <- t_test_auc$conf.int[2]
  peak_lower <- t_test_auc$conf.int[1]
  res_mat[iter, 16] <- mean( I( beta_avg <= peak_upper & beta_avg >= peak_lower) ) # whether CI contains estimate
  res_mat[iter, 39] <- I(t_test_auc$p.value <= 0.05)
  rm(beta_avg)
  
  ##################################
  # permutation test
  ##################################
  # compare over entire trial (including baseline/pre-cue)
  perm_reward_period_idx <- photo_indices #[seq(pre_min_tm * target_Hz + 1, length(photo_indices) )] # indices acounting for non-photometry variables (e.g., ID, cue)
  L <- length(perm_reward_period_idx)
  perm_mat <- matrix(nrow = boots, ncol = L)
  perm_vec <- vector(length = boots)
  ids_samp <- unique(dat_sim$data$ID)
  n <- length(ids_samp)
  cue_vals <- unique(dat_sim$data$cue)
  idx_list <- lapply(ids_samp, function(x) which(dat_sim$data$ID == x))
  
  for(j in 1:boots){
    samp <- sample.int(n, size = n, replace = TRUE)
    idx <- do.call(c, lapply(1:n, function(x) idx_list[[samp[x]]]))
    d <- dat_sim$data[idx, ]
    perm_mat[j,] <- colMeans(d[d$cue == cue_vals[2], perm_reward_period_idx]) - 
      colMeans(d[d$cue == cue_vals[1], perm_reward_period_idx])
    
    perm_vec[j] <- mean(perm_mat[j,reward_period])
  }
  
  rm(d, idx, samp, idx_list, ids_samp)
  
  # CIs
  d <- dat_sim$data
  beta_perm <- colMeans(as.matrix(d[d$cue == cue_vals[2], perm_reward_period_idx])) - 
                colMeans(as.matrix(d[d$cue == cue_vals[1], perm_reward_period_idx])) # average difference
  perm_CIs <- apply(perm_mat, 2, function(x) quantile(x, probs = c(0.025, 0.975)) ) #matrixStats::colSds(perm_mat)
  perm_low <- perm_CIs[1,] * sqrt(n / (n-1))
  perm_up <- perm_CIs[2,] * sqrt(n / (n-1))
  
  res_mat[iter, 17] <- sqrt( mean( (dat_sim$beta[beta_idx,] - beta_perm )^2) )
  res_mat[iter, 18] <- mean( I(dat_sim$beta[beta_idx,] <= perm_up & 
                                 dat_sim$beta[beta_idx,] >= perm_low) ) # whether CI contains estimate
  res_mat[iter, 26] <- mean( I(dat_sim$beta[beta_idx,reward_period] <= perm_up[reward_period] & 
                                 dat_sim$beta[beta_idx,reward_period] >= perm_low[reward_period] ) ) # whether CI contains estimate over reward period
  
  
  # mean difference
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  mean_diff <- mean(as.matrix(d[d$cue == cue_vals[2], perm_reward_period_idx])) - 
                  mean(as.matrix(d[d$cue == cue_vals[1], perm_reward_period_idx]))
  auc_up <- quantile(perm_vec, probs = c(0.025, 0.975))[2] * sqrt(n / (n-1))
  auc_low <- quantile(perm_vec, probs = c(0.025, 0.975))[1] * sqrt(n / (n-1))
  res_mat[iter, 49] <- sqrt( mean(  (beta_avg - mean_diff )^2) )
  res_mat[iter, 50] <- mean( I( beta_avg <= auc_up & beta_avg >= auc_low ) )  # whether CI contains estimate
  
  # significance of auc
  res_mat[iter, 41] <- 1 * I(auc_up * auc_low > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
  
  # average of significance
  perm_signif_vec <- 1 * I(perm_up[reward_period] * perm_low[reward_period] > 0)
  res_mat[iter, 13] <- mean(perm_signif_vec) # naive without correction of consecutive threshold
  res_mat[iter, 9] <-  mean( cluster_finder(r = perm_signif_vec, min.length = round(target_Hz/4)) ) # consecutive threshold
  res_mat[iter, 14] <- mean( cluster_finder(r = perm_signif_vec, min.length = round(target_Hz/2)) ) # consecutive threshold
  res_mat[iter, 51] <- mean( cluster_finder(r = perm_signif_vec, min.length = round(target_Hz)) ) # consecutive threshold
  rm(d)
  #--------------------------------------------------------------------------------------------------------------------
  
  ###################################
  # bootstrap fLME
  ###################################
  rm(fit)
  ##############################
  # fit model to simulated data
  ##############################
  if(flme_boot){  
    set.seed(1)
    num_boots <- 1000 
    # fit model and get model CIs
    timeStart <- Sys.time()
    fit_boot <- fui(photometry ~ cue +
                        (cue | ID), 
                      data = dat_sim$data, 
                      family = "gaussian", 
                      var = TRUE, 
                      parallel = FALSE, 
                      silent = TRUE,
                      analytic = FALSE,
                      nknots_min = nknots_min,
                      nknots_min_cov = nknots_min_cov,
                      smooth_method = "GCV.Cp",
                      splines = "tp",
                      design_mat = FALSE,
                      subj_ID = "ID",
                      boot_type = "reb",
                      num_boots = num_boots,
                      num_cores = num_cores)
    timeEnd <- Sys.time() 
    
    res_mat[iter,47] <- as.numeric(difftime(timeEnd, timeStart, units = "mins"))
    
    ##############################
    # generate 95% CIs
    ##############################
    lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = nrow(fit_boot$betaHat), ncol = ncol(fit_boot$betaHat)) # initialize CIs
    joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
    lower.avg <- upper.avg <- lower # CIs for averaging
    
    # iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
    for(r in 1:p){
      # joint CIs
      lower.joint[r,] = fit_boot$betaHat[r,] - fit_boot$qn[r]*sqrt(diag(fit_boot$betaHat.var[,,r]))
      upper.joint[r,] = fit_boot$betaHat[r,] + fit_boot$qn[r]*sqrt(diag(fit_boot$betaHat.var[,,r]))
      
      # check whether estimate in CI
      joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
      
      # naive CIs
      lower[r,] = fit_boot$betaHat[r,] - 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]))
      upper[r,] = fit_boot$betaHat[r,] + 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]))
      
      # naive CIs for averaging (remove the multiplying the variance by 1.2 correction for correlation)
      lower.avg[r,] = fit_boot$betaHat[r,] - 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]) / 1.2)
      upper.avg[r,] = fit_boot$betaHat[r,] + 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]) / 1.2)
      
      # check whether estimate in CI
      naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
    }
    
    # average of CI inclusion over reward period
    res_mat[iter, 25] <- mean( joint_incl[beta_idx, reward_period] )  # whether joint CI contains estimate
    res_mat[iter, 30] <- mean( naive_incl[beta_idx, reward_period] )  # whether pointwise CI contains estimate
    
    # calculate regression coefficient inclusion in 95% CI (across functional domain)
    joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
    naive_incl <- rowSums(naive_incl) / ncol(naive_incl)
    
    # average CI coverage across regression coefficients
    res_mat[iter, 22] <- mean(joint_incl)
    res_mat[iter, 27] <- mean(naive_incl)
  
    # don't evaluate performance on intercept (so only cue effect here since one covariate)
    res_mat[iter,23] <- mean(joint_incl[-1])
    res_mat[iter,65] <- mean(joint_incl[1])
    res_mat[iter,28] <- mean(naive_incl[-1])
    
    # average of significance
    res_mat[iter, 33] <- mean( 1 * I(upper[beta_idx, reward_period] * lower[beta_idx, reward_period] > 0) )
    res_mat[iter, 34] <- mean( 1 * I(upper.joint[beta_idx, reward_period] * lower.joint[beta_idx, reward_period] > 0) )
    
    # significance of AUC
    res_mat[iter, 44] <- 1 * I(mean(upper[beta_idx, reward_period]) * mean(lower[beta_idx, reward_period]) > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
    res_mat[iter, 53] <- 1 * I(mean(upper.avg[beta_idx, reward_period]) * mean(lower.avg[beta_idx, reward_period]) > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
    
    # closed form expression for variance of average of betas
    len <- rep(1/length(reward_period), length(reward_period))
    V <- fit_boot$betaHat.var[reward_period,reward_period,beta_idx]
    V <- eigenval_trim(V)   # trim eigenvalues to ensure matrix is PSD
    avg_var <- t(len) %*% V %*% len 
    # CIs
    avg_up <- mean(fit_boot$betaHat[beta_idx, reward_period]) + 1.96 * sqrt(avg_var)
    avg_low <- mean(fit_boot$betaHat[beta_idx, reward_period]) - 1.96 * sqrt(avg_var)
    # significance of average
    res_mat[iter, 42] <- 1 * I(avg_up * avg_low > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
    res_mat[iter, 43] <- mean( I( beta_avg <= avg_up) & beta_avg >= mean(avg_low) )  # whether CI contains estimate
  
  # AUC beta
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  res_mat[iter, 54] <- mean( I( beta_avg <= mean(upper.avg[beta_idx, reward_period]) & beta_avg >= mean(lower.avg[beta_idx, reward_period]) ) ) # whether CI contains estimate
  
  rm(beta_avg)
  
  ###############################################
  # generate 95% CIs with sample size correction
  ###############################################
  lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = nrow(fit_boot$betaHat), ncol = ncol(fit_boot$betaHat)) # initialize CIs
  joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix
  lower.avg <- upper.avg <- lower # CIs for averaging
  n <- length(unique(dat_sim$data$ID))
  
  # iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
  for(r in 1:p){
    # joint CIs
    lower.joint[r,] = fit_boot$betaHat[r,] - fit_boot$qn[r]*sqrt(diag(fit_boot$betaHat.var[,,r])) * sqrt(n / (n-1))
    upper.joint[r,] = fit_boot$betaHat[r,] + fit_boot$qn[r]*sqrt(diag(fit_boot$betaHat.var[,,r])) * sqrt(n / (n-1))
    
    # check whether estimate in CI
    joint_incl[r,] <- dat_sim$beta[r,] <= upper.joint[r,] & dat_sim$beta[r,] >= lower.joint[r,]
    
    # naive CIs
    lower[r,] = fit_boot$betaHat[r,] - 1.96*sqrt(diag(fit_boot$betaHat.var[,,r])) * sqrt(n / (n-1))
    upper[r,] = fit_boot$betaHat[r,] + 1.96*sqrt(diag(fit_boot$betaHat.var[,,r])) * sqrt(n / (n-1))
    
    # naive CIs for averaging (remove the multiplying the variance by 1.2 correction for correlation)
    lower.avg[r,] = fit_boot$betaHat[r,] - 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]) / 1.2) * sqrt(n / (n-1))
    upper.avg[r,] = fit_boot$betaHat[r,] + 1.96*sqrt(diag(fit_boot$betaHat.var[,,r]) / 1.2) * sqrt(n / (n-1))
    
    # check whether estimate in CI
    naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
  }
  
  # average of CI inclusion over reward period
  res_mat[iter, 69] <- mean( joint_incl[beta_idx, reward_period] )  # whether joint CI contains estimate
  res_mat[iter, 70] <- mean( naive_incl[beta_idx, reward_period] )  # whether pointwise CI contains estimate
  
  # calculate regression coefficient inclusion in 95% CI (across functional domain)
  joint_incl <- rowSums(joint_incl) / ncol(joint_incl)
  naive_incl <- rowSums(naive_incl) / ncol(naive_incl)
  
  # average CI coverage across regression coefficients
  res_mat[iter, 71] <- mean(joint_incl)
  res_mat[iter, 72] <- mean(naive_incl)
  
  # don't evaluate performance on intercept (so only cue effect here since one covariate)
  res_mat[iter,73] <- mean(joint_incl[-1])
  res_mat[iter,74] <- mean(joint_incl[1])
  res_mat[iter,75] <- mean(naive_incl[-1])
  
  # average of significance
  res_mat[iter, 76] <- mean( 1 * I(upper[beta_idx, reward_period] * lower[beta_idx, reward_period] > 0) )
  res_mat[iter, 77] <- mean( 1 * I(upper.joint[beta_idx, reward_period] * lower.joint[beta_idx, reward_period] > 0) )
  
  # significance of AUC
  res_mat[iter, 78] <- 1 * I(mean(upper[beta_idx, reward_period]) * mean(lower[beta_idx, reward_period]) > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
  res_mat[iter, 79] <- 1 * I(mean(upper.avg[beta_idx, reward_period]) * mean(lower.avg[beta_idx, reward_period]) > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
  
  # closed form expression for variance of average of betas
  len <- rep(1/length(reward_period), length(reward_period))
  V <- fit_boot$betaHat.var[reward_period,reward_period,beta_idx]
  V <- eigenval_trim(V)   # trim eigenvalues to ensure matrix is PSD
  avg_var <- t(len) %*% V %*% len 
  # CIs
  avg_up <- mean(fit_boot$betaHat[beta_idx, reward_period]) + 1.96 * sqrt(avg_var) * sqrt(n / (n-1))
  avg_low <- mean(fit_boot$betaHat[beta_idx, reward_period]) - 1.96 * sqrt(avg_var) * sqrt(n / (n-1))
  # significance of average
  res_mat[iter, 80] <- 1 * I(avg_up * avg_low > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
  beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
  res_mat[iter, 81] <- mean( I( beta_avg <= avg_up) & beta_avg >= mean(avg_low) )  # whether CI contains estimate
  
  # AUC beta
  res_mat[iter, 82] <- mean( I( beta_avg <= mean(upper.avg[beta_idx, reward_period]) & beta_avg >= mean(lower.avg[beta_idx, reward_period]) ) ) # whether CI contains estimate
  
  rm(beta_avg)
  }
}
#-------------------------------------------------------------------
# refund
#-------------------------------------------------------------------
dat <- tibble(photometry = Y,
              cue = as.numeric(dat_sim$data$cue),
              ID = as.factor(dat_sim$data$ID))


timeStart <- Sys.time()
fit_pffr <- refund::pffr(formula = photometry ~ cue + s(ID, bs="re") + s(ID, cue, bs="re"),
                         algorithm = "bam",
                         family = "gaussian",
                         bs.yindex = list(bs = "ps", k = nknots_min, m = c(2, 1)),
                         data = dat)

timeEnd <- Sys.time()

##############################
# generate 95% CIs
##############################
# pffr 
coef_pffr <- coef(fit_pffr, n1 = L)
betaHat_pffr <- betaHat_pffr.var <- dat_sim$beta
betaHat_pffr[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$value)
betaHat_pffr[2,] <- as.vector(coef_pffr$smterms$`cue(yindex)`$value)
betaHat_pffr.var[1,] <- as.vector(coef_pffr$smterms$`Intercept(yindex)`$se^2)
betaHat_pffr.var[2,] <- as.vector(coef_pffr$smterms$`cue(yindex)`$se^2)

lower <- upper <- lower.joint <- upper.joint <- matrix(NA, nrow = nrow(dat_sim$beta), ncol = ncol(dat_sim$beta)) # initialize CIs
joint_incl <- naive_incl <- lower # initialize inclusion indicator matrix

# iterate across regression coefficients, calculate PIs and determine inclusion across fn. domain
for(r in 1:p){
  # naive CIs
  lower[r,] = betaHat_pffr[r,] - 1.96*sqrt(betaHat_pffr.var[r,])
  upper[r,] = betaHat_pffr[r,] + 1.96*sqrt(betaHat_pffr.var[r,])
  
  # check whether estimate in CI
  naive_incl[r,] <- dat_sim$beta[r,] <= upper[r,] & dat_sim$beta[r,] >= lower[r,]
}

# average CI inclusion over reward period
res_mat[iter, 29] <- mean( naive_incl[beta_idx, reward_period] )  # whether pointwise CI contains estimate

# calculate regression coefficient inclusion in 95% CI (across functional domain)
naive_incl <- rowSums(naive_incl) / ncol(naive_incl)

# average CI coverage across regression coefficients
res_mat[iter,55] <- mean(naive_incl)
res_mat[iter,56] <- sqrt( mean( (dat_sim$beta - betaHat_pffr)^2 ) ) # rmse of model fit
res_mat[iter,57] <- as.numeric( difftime(timeEnd, timeStart, units="mins") ) # time to fit model

# don't evaluate performance on intercept (so only cue effect here since one covariate)
res_mat[iter,58] <- mean(naive_incl[-1])
res_mat[iter,66] <- mean(naive_incl[1])
res_mat[iter,59] <- sqrt( mean( (dat_sim$beta[-1,] - betaHat_pffr[-1,])^2 ) ) # rmse of model fit
res_mat[iter,60] <- mean(dat_sim$beta[beta_idx,] - betaHat_pffr[beta_idx,]) # bias of cue term

# average beta
beta_avg <- mean( dat_sim$beta[beta_idx, reward_period] ) # average beta for cue difference
res_mat[iter, 61] <- mean( I( beta_avg <= mean(upper[beta_idx, reward_period]) & beta_avg >= mean(lower[beta_idx, reward_period]) ) ) # whether CI contains estimate
res_mat[iter, 62] <- 1 * I(mean(upper[beta_idx, reward_period]) * mean(lower[beta_idx, reward_period]) > 0) # if upper and lower CIs are either both positive or negative, then doesn't contain 0 so significant
res_mat[iter,63] <- sqrt( mean( (beta_avg - mean(betaHat_pffr[beta_idx,reward_period]))^2 ) ) # rmse of model fit

# average of significance
res_mat[iter, 52] <- mean( 1 * I(upper[beta_idx, reward_period] * lower[beta_idx, reward_period] > 0) )
#--------------------------------------------------------------------------------------------------------------------

####################################
# save to CSV on cluster
####################################
# saveFn_new <- function(file, fileNm, iterNum, iters, save.folder = NA)
  saveFn_new(file = res_mat,
       fileNm = fileNm,
       iterNum = iter,
       iters = reps,
       save.folder = save.folder)

Sys.sleep(15 * 60) # ensure no issues with cluster