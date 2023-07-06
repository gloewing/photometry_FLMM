# simulations using Science paper's data for data-driven sims
# this R code is saved in: /gpfs/gsfs8/users/loewingergc/photometry_fglmm/code

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

# params matrix
# simParams <- expand.grid(TRUE, c(TRUE, FALSE), TRUE, c(4:12) )
# colnames(simParams) <- c("fix_smooth", "resid_var_indep", "resid_var_subj", "n")
# write.csv( simParams, "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photo_params.csv", row.names = FALSE)

cluserInd <- FALSE # whether running on computer or on cluster
num_cores <- 10 # parallelize

simParams <- read.csv( "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photo_params.csv" )

runNum <- 1
iter <- 2
data_path <-"/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper" # path to original Matlab files for data
source("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photometry_sim_fLME_fn_multi.R") # simulations
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_wrap.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM.R')
wd <- "~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/"

# simulation parameters
fixed_smooth <- simParams[runNum, 1] #TRUE # whether to use smoothed coefficients as truth or not
reps <- 100 # replicates of simulations
boots <- 5000 # for naive bootstrap permutation test
knots_div <- 2 # number of functional samples divisor (so if it is 1, there is a knot for every obsveration on fn domain)
resid_var_indep <- simParams[runNum, 2] # FALSE # assume residual variance is indepedent across fn. domain as in paper
resid_var_subj <- simParams[runNum, 3] #TRUE # residual variance subject specific or global
n_sim <- 7 #simParams[runNum, 4] #8
comparison_methods <- FALSE # whether to fit other models besides fLME
trial_trim <- 3 # take every other trial to reduce computational burden

fileNm <- paste0("photometry_sims_science_lick_boot_",
                 "_trl_trm_", trial_trim,
                   "fixed_smooth_", fixed_smooth,
                   "_knots_div_", knots_div,
                   "_resid_indep_", resid_var_indep,
                   "_resid_subj_", resid_var_subj,
                   "_n_trgt_", n_sim
                   )

  
res_mat <- matrix(NA, ncol = 32, nrow = reps)
colnames(res_mat) <- c("beta_CI_joint", "beta_CI_naive", "fLME_beta_rmse", "fLME_time",
                        "fLME_Maxbeta_rmse", "fLME_Maxbeta_CI_incl", "fLME_avgbeta_rmse", "fLME_avgbeta_CI_incl",
                        "LME_Maxbeta_rmse", "LME_Maxbeta_CI_incl", "LME_avgbeta_rmse", "LME_avgbeta_CI_incl",
                        "ttest_Maxbeta_rmse", "ttest_Maxbeta_CI_incl", "ttest_avgbeta_rmse", "ttest_avgbeta_CI_incl",
                        "perm_rmse", "perm_CI_incl",
                        "beta_CI_joint_noIntrcpt", "beta_CI_naive_noIntrcpt", "fLME_beta_rmse_noIntrcpt",
                        "beta_CI_boot", "beta_CI_boot_noIntrcpt", "fLME_Maxbeta_CI_incl_boot", "fLME_avgbeta_CI_incl_boot", "fLME_boot_time",
                        "beta_CI_boot_naive", "beta_CI_boot_noIntrcpt_naive", "fLME_Maxbeta_CI_incl_boot_naive", "fLME_avgbeta_CI_incl_boot_naive",
                       "fLME_bias", "fLME_simple_bias")

#################################
# original data variables
#################################

mouseList = c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-F1','HJ-FP-F2','HJ-FP-M6','HJ-FP-M7','HJ-FP-M8')
n <- length(mouseList)
session_max <- 11 # maximum session number (varies between animals)
session_min <- 1
min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
pre_reward_period_length <- 1.5 # length of pre-reward period
pre_lick_reward_period <- 0.5
reward_period_length <- 1.5
pre_min_tm <- pre_reward_period_length + pre_lick_reward_period # pre_lick period
post_lick_reward_period <- 1 # 1 second after 
post_reward_length <- 1.5 # 2 seconds after reward period ends
post_min_tm <- post_reward_length + post_lick_reward_period # post lick period
total_trial_length <- pre_min_tm + post_min_tm # total trial length
trial_num <- 100 # trials per session
ibi_length <- 1 # 1 second inter-lick-interval for bout calculator
trial_length <- post_min_tm + pre_min_tm # total length of session in seconds

# samples before after
Hz <- 1 / 0.008
target_Hz <- 25
pre_samps <- round(pre_min_tm * Hz) # samples before lick
post_samps <- round(post_min_tm * Hz) # samples after lick

if( !cluserInd ){
  # fGLMM code
  source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
  source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_smooth.R") # function code
  source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_multi.R") # function code
  source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_wrap.R") # function code
  
  path <- "~/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_3/"
  
  # event indices
  index_lick = 5
  index_bgdrw = 7
  index_sessionend = 0
  
  dat_lst <- vector(length = session_max * n, "list")
  cnt <- 0
  
  path <- "~/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_3/"
  mouseList = c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-F1','HJ-FP-F2','HJ-FP-M6','HJ-FP-M7','HJ-FP-M8')
  n <- length(mouseList)
  session_max <- 11 # maximum session number (varies between animals)
  session_min <- 1
  min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
  pre_reward_period_length <- 1.5 # length of pre-reward period
  pre_lick_reward_period <- 0.5
  reward_period_length <- 1.5
  pre_min_tm <- pre_reward_period_length + pre_lick_reward_period # pre_lick period
  post_lick_reward_period <- 1 # 1 second after 
  post_reward_length <- 1.5 # 2 seconds after reward period ends
  post_min_tm <- post_reward_length + post_lick_reward_period # post lick period
  total_trial_length <- pre_min_tm + post_min_tm # total trial length
  trial_num <- 100 # trials per session
  ibi_length <- 1 # 1 second inter-lick-interval for bout calculator
  trial_length <- post_min_tm + pre_min_tm # total length of session in seconds
  
  # samples before after
  Hz <- 1 / 0.008
  target_Hz <- 25
  pre_samps <- round(pre_min_tm * Hz) # samples before lick
  post_samps <- round(post_min_tm * Hz) # samples after lick
  
  # event indices
  index_lick = 5
  index_bgdrw = 7
  index_sessionend = 0
  
  dat_lst <- vector(length = session_max * n, "list")
  cnt <- 0
  
  for(i in 1:n){
    for(s in session_min:session_max){
      
      id <- mouseList[i] # id name
      
      # file name
      fl_nm <- paste0("sub-", id, "_ses-Day", s, ".nwb")
      
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
        dat_photo <- dat_photo[!is.na(dat_photo$data), ] # remove NAs
        
        #########################
        # inter trial interval 
        #########################
        # make sure IRI is long enough 
        iri_length <- mean(diff(ttls$time[ttls$event == index_bgdrw]))
        num_trials <- length(ttls$time[ttls$event == index_bgdrw])
        
        # -- first iri is difference between when photometry starts and first actual IRI just the time in trial when it occurs
        iri <- c(ttls$time[ttls$event == index_bgdrw][1],  diff(ttls$time[ttls$event == index_bgdrw]) ) # IRIs (assuming first event is removed) -- use IRIs before aligning time to photometry 
        #iri <- diff( ttls$time[ttls$event == index_bgdrw] )  # IRIs (assuming first event is removed) -- use IRIs before aligning time to photometry 
        iri <- round(iri, 4)
        
        # variance of mean of 100 trials is (100 / (1/12)^2 ) / 100^2 = 1.44, so 8 seconds is > 3 * sd less than mean
        if(iri_length > 8){
          
          cnt <- cnt + 1 # counter for list
          
          # align times of behavior to photometry times
          ttls <- time_align(time_truth = dat_photo$timestamps, data = ttls, name = "time") #align times in ttls to the times in photometry data
          
          reward_times <- ttls$time[ which(ttls$event == index_bgdrw) ] # reward times from behavior
          
          # number of behavior trials before photometry starts
          trials_removed <- sum( I( ttls$time <= min(dat_photo$timestamps)  &  ttls$event == index_bgdrw ) ) # indicator vector to remove trials before photometry
          
          # reward times
          reward_times <- ttls$time[ which(ttls$event == index_bgdrw) ] # reward times from behavior
          reward_idx <- sapply(reward_times, function(x) which(x == dat_photo$timestamps) ) # reward indices in photometry
          trial_iri <-  which(diff(reward_times) < min_tm) # which trials to remove
          trial_iri <-  unique(c(trial_iri, trial_iri + 1)) # remove trials after those too because their pre-peroid is not long enough
          incl_ind <- rep(TRUE, num_trials)
          incl_ind[trial_iri] <- FALSE
          trial_iri <- incl_ind
          if(trials_removed > 0)       trial_iri[1:trials_removed] <- FALSE # remove initial trials if photometry starts after reward
          
          # determine whether enough recorded photometry for pre-cue peroid for each reward
          photo_start_time <- min(dat_photo$timestamps)
          photo_stat <- I( reward_times -  min(dat_photo$timestamps) > pre_min_tm)
          trial_iri <- I(trial_iri * photo_stat == 1)
          
          # lick times
          lick_times <- ttls$time[ which(ttls$event == index_lick) ] # reward times from behavior
          lick_idx <- sapply(lick_times, function(x) which(x == dat_photo$timestamps) ) # lick indices in photometry
          tot_licks <- sapply(reward_times, function(x)  sum(lick_times > x & lick_times <= (x + min_tm) ) ) # number of licks in post-reward period
          first_lick <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1] - x ) # time-to-first lick latency in post-reward period
          trial_iri <- I(trial_iri * I(tot_licks > 0) == 1) # make sure there is at least 1 lick
          
          # consummatory licks
          first_lick_session <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1]) # session time of first lick in post-reward period in overall session time
          lick_bouts <- bout_finder_trial(x = lick_times, Hz = Hz, ibi = ibi_length, bout_start_time = first_lick_session, include = trial_iri)
          
          # first lick index
          first_lick_idx <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1] ) # number of licks in post-reward period
          first_lick_idx <- sapply(first_lick_idx, function(x) which(x == dat_photo$timestamps))
          first_lick_idx <- sapply(first_lick_idx, function(x) ifelse(length(x) > 0, x, NA)) # make sure to have NAs in empty places
          first_lick <- round(first_lick, 4)  
          
          # center around licks
          idx <- t( sapply(first_lick_idx[trial_iri], function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for trials
          dat <- as.data.frame( apply(idx, 2, function(x) dat_photo$data[x]) ) # photometry trials
          L <- ncol(dat)
          
          # ******************************
          # normalize trials based on pre-baseline
          dat <- dat - rowMeans(dat[,1:(post_reward_length * Hz)])
          # ******************************
          
          # downsample photometry
          by_num <- round( Hz / target_Hz )
          photo_idx_pre <- sort(-seq(-pre_samps, -1, by = by_num)) # seq(1, pre_samps, by = by_num)
          photo_idx_post <- seq(pre_samps + by_num, L, by = by_num)
          photo_idx <- unique( c(photo_idx_pre, photo_idx_post) )
          dat <- dat[,photo_idx]
          L <- ncol(dat)
          colnames(dat) <- paste0("photometry.", 1:L)
          
          # new samples before/after
          pre_samps_target <- length(photo_idx_pre) # samples before lick with new sampling rate
          post_samps_target <- length(photo_idx_post) # samples after lick with new sampling rate
          
          # combine data
          trials <- seq(1, num_trials )[trial_iri]
          sess_info <- as.data.frame( cbind(id, s, trials, tot_licks[trial_iri], first_lick[trial_iri], iri[trial_iri], lick_bouts[trial_iri,3:4]) )
          colnames(sess_info) <- c("id", "session", "trial", "licks", "lick_time", "iri", "bout_events", "bout_rate")
          
          dat_lst[[cnt]] <- as.data.frame( cbind(sess_info, dat) )
          
          rm(dat, trials, sess_info, ttls, dat_photo, idx, tot_licks, first_lick, iri, trials_removed, lick_bouts)

        }
        
      }
      
    }
    
  }
  
  # cocnatenate results
  dat <- do.call(rbind, dat_lst)
  dat <- dat[!is.na(dat$bout_rate),] # remove trials with problematic bouts (1.4% of trials)

  dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
  dat2 <- dat
  rm(dat_lst)
  #####################################################################################
  
  # photometry sample length
  out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
  L <- length(out_index)
  nknots_min <- round(L/2)
  cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates
  
  #########################################################################################
  # creat AUC variables for EDA
  
  # session_max <- 11 # maximum session number (varies between animals)
  # session_min <- 1
  # min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
  # pre_reward_period_length <- 1.5 # length of pre-reward period
  # pre_lick_reward_period <- 0.5
  # reward_period_length <- 1.5
  # pre_min_tm <- pre_reward_period_length + pre_lick_reward_period # pre_lick period
  # post_lick_reward_period <- 1 # 1 second after 
  # post_reward_length <- 1.5 # 2 seconds after reward period ends
  # post_min_tm <- post_reward_length + post_lick_reward_period # post lick period
  # total_trial_length <- pre_min_tm + post_min_tm # total trial length
  # trial_num <- 100 # trials per session
  # ibi_length <- 1 # 1 second inter-lick-interval for bout calculator
  # trial_length <- post_min_tm + pre_min_tm # total length of session in seconds
  
  
  photo_dat <- dat[, -cov_idx] # photometry data
  pre_lick_reward_period_samps <- round(target_Hz * pre_reward_period_length) # samples of pre_cue period
  reward_period_max_idx <- round( (pre_reward_period_length + reward_period_length) * target_Hz ) # max_idx of post_cue period
  dat_AUC_pre <- rowMeans( photo_dat[, 1:pre_lick_reward_period_samps] ) # AUC before reward period
  dat_AUC_post <- rowMeans(photo_dat[, seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)]) # AUC after reward period
  
  
  
  # scale covariates to make contrasts more interpretable
  dat$trial <- dat$trial - trial_num / 2 # trial 50 is middle
  dat$session <- scale(dat$session, center = TRUE, scale = FALSE) # consider mroe interpertable version
  dat$licks <- scale(dat$licks, center = TRUE, scale = FALSE)
  dat$lick_time <- scale(dat$lick_time, center = TRUE, scale = FALSE)
  dat$iri <- scale(dat$iri, center = TRUE, scale = TRUE)
  dat$reward_num <- dat$trial + (dat$session - 1) * num_trials # number of rewards assuming 100 rewards/session
  dat$dat_AUC_pre <- dat_AUC_pre
  dat$dat_AUC_post <- dat_AUC_post
  dat$dat_AUC <- dat_AUC_post - dat_AUC_pre # described in "Experiment 1" of "Data Analysis" in Supplement
  dat <- dat[!is.infinite(dat$bout_rate),] # 1 trial has infinite bout rate
  dim(dat)
  dat <- dat[complete.cases(dat),] # remove one observation
  
}else{
  dat <- read.csv(paste0(data_path, "/science_lick_align.csv"))
}

# trial trim to reduce computational burden
trls_reduced <- as.integer(round( seq(min(dat$trial), max(dat$trial), by = trial_trim)   ))
dat <- dat[dat$trial %in% trls_reduced, ] # reduce number of trials analzed

# number of knots
nknots_min <- 40 # L = 100

# order data
dat_photo <- dat[order(dat$id, dat$session, dat$trial, decreasing = FALSE), ] %>%
                dplyr::as_tibble() %>%
                dplyr::rename(cue = iri) %>%
                as.data.frame()
n <- length(unique(dat_photo$id)) # number of subjects in (potentially) reduced sample
rm(dat)

# scale variables for stability
dat_photo$trial <- scale(dat_photo$trial)
dat_photo$cue <- scale(dat_photo$cue) # already was close to being scaled

############################################
# fit original model for ground truth
############################################
# fit model
timeStart <- Sys.time()

fit <- lfosr(formula = photometry ~ cue +
                  (cue | id), # thi sis equivalent to A + B
                   data = dat_photo, 
                   family = "gaussian", 
                   var = TRUE, 
                   parallel = TRUE, 
                   silent = FALSE,
                   analytic = TRUE,
                   nknots_min = nknots_min,
                   smooth_method = "GCV.Cp",
                   splines = "tp",
                   design_mat = TRUE,
                   G_return = TRUE,
                   residuals = TRUE,
                   wd = wd,
                   subj_ID = "id",
                   num_cores = num_cores,
                   multi_level = TRUE)

timeEnd <- Sys.time()

# pre_cue_idx <- seq(1, target_Hz * cue_analyze_start) # indices of photometry samples before cue (to remove) so we can derive summary statistics from trial and analyze
# colMeans(fit$aic[-pre_cue_idx,])

# number of covariates
p <- nrow(fit$betaHat)

##################
# simulate data
##################
set.seed(iter)

dat_sim <- photo_stimulate(X = NULL, # design matrix
                           Z = NULL,
                           include = dat_photo["trial"],
                           N = n_sim,
                           model = fit, # fLME model object
                           fixed_smooth = fixed_smooth, # whether to use smoothed betas for true fixed effects
                           IDs = dat_photo$id,  # col of IDs (N x 1 vector)
                           resid_var_indep = resid_var_indep, # assume residual variance is indepedent across fn. domain as in paper
                           resid_var_subj = resid_var_subj, # calculate error variance in subject-specific fashion as opposed to globally
                           fixed_effects_draw = FALSE,
                           fixed_effects_scale = 1,
                           random_effects_scale = 1,
                           rand_int = FALSE,
                           rand_int_sigma = FALSE, #rand_int_sigma,
                           ID_level_re = TRUE # indicates that random effects are subject-level only (do not vary within subject)
                           )
#-------------------------------------------------------------------------------------------------------
# individual observation plots
#-------------------------------------------------------------------------------------------------------
# filter data average over trials to get Session Avg.s for each animal
mean_data <- dat_sim$data %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>%
  dplyr::filter(cue == 1) %>%  # baseline
  #               group == "control", # control animals
  #               seshID == sess_max, # last session so well trained
  #               ids %in% ids_inc, # only animals with sufficient trials (only eliminates one animal)
  #               lickState == 1) %>%  # only look at lick+
  dplyr::mutate(time = as.numeric(time) / target_Hz - pre_min_tm) %>%
  dplyr::group_by(time, ID) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) )

# first plot average for animal 1 +/- s.e.m.
p6 <- mean_data %>%
  ggplot(aes(x = as.numeric(time), y = photometry, ymin=y_low, ymax=y_high, fill = as.factor(ID) )) + # ,  fill = "#ca0020"      , 
  #geom_ribbon(aes(fill = "#ca0020"), alpha=0.3, colour = NA) + 
  geom_line(size = 0.75, aes(colour = as.factor(ID))) +  #, colour = "#ca0020"
  scale_color_grey(start = 0.7, end = 0.1) +
  scale_fill_grey(start = 0.7, end = 0.1) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Simulated Long-Delay Animal Averages") +
  #coord_cartesian(ylim = c(-0.75, 3.25)) +
  # scale_fill_manual(values = c("#ca0020") ) + 
  theme_classic() + 
  geom_vline(xintercept= 3,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  # geom_vline(xintercept=300/ target_Hz - pre_min_tm, 
  #            linetype="dashed", 
  #            color = "gray", 
  #            size = rel(0.3),
  #            alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_text(size=12, face="bold"),
        axis.title.y = element_text(size=12, face="bold"),
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none")
# ------------------------------------------

#######################################
# individual traces
#######################################
# filter data
mean_data <- dat_sim$data %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(ID == 1 ) %>%  # arbitrarily choose the first control animal
  dplyr::mutate(time = as.numeric(time) / target_Hz - pre_min_tm) %>%
  dplyr::group_by(time, cue, trial) 

mean_data$cue <- ifelse(mean_data$cue == 0, "Short", "Long")
# # average over cue1 trials
# mean_dat1 <- mean_data  %>%
#   dplyr::filter(cue == 1 ) %>% 
#   dplyr::summarise(photometry = mean(value, na.rm = TRUE),
#                    y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
#                    y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) # +/- sem
#   )

# average over cue0 trials
mean_dat0 <- mean_data  %>%
  # dplyr::filter(cue == 0 ) %>% 
  group_by(cue, time) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) # +/- sem
  )

# filter and choose specific trials
mean_dat2 <- mean_data %>% 
  group_by(cue, time) %>%
  dplyr::slice(c(1, 25, 50)) %>% # subset of rows
  dplyr::rename(photometry = value)

# add trials
# mean_dat2$trial_num <- rep(rep(1:50, each = 50), 2)

# first plot average for animal 1 +/- s.e.m.
p4 <- 
  ggplot(NULL, aes(x = time, y = photometry)) + # ,  fill = "#ca0020"      , 
  # geom_ribbon(data = mean_dat1, aes(ymin=y_low, ymax=y_high, fill = "#ca0020"), alpha=0.5, colour = NA) + 
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 1) +  #, colour = "#ca0020"
  # scale_color_grey(start = 0.7, end = 0.1) +
  geom_line(data = mean_dat2, aes(colour = as.factor(cue), group = trial), size = 0.5, alpha = 0.5) +  #, colour = "#ca0020"
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Trial-Level Data") +
  # coord_cartesian(ylim = c(-2.5, 4)) +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  geom_vline(xintercept= 3,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_text(size=14, face="bold", hjust=0.5),# element_blank(),
        axis.title.y = element_text(size=12, face="bold")) + 
  # legend.position="none"
  guides(color= guide_legend(title="Delay"), fill = "none")




#----------------------------------------------------------
# A) functional LME
#----------------------------------------------------------
beta_idx <- 2 # this corresponds to cue effect --only evaluate performance of models on these to make comparable with t-tests

##############################
# fit model to simulated data
##############################

# fit model and get model CIs
timeStart <- Sys.time()
fit <- lfosr(photometry ~ cue * trial +
                 (trial| ID/session), # thi sis equivalent to A + B 
               data = dat_sim$data, 
               family = "gaussian", 
               var = TRUE, 
               parallel = FALSE, 
               silent = FALSE,
               analytic = TRUE,
               nknots_min = nknots_min,
               smooth_method = "GCV.Cp",
               splines = "tp",
               design_mat = TRUE,
               wd = wd,
               subj_ID = "ID",
               num_cores = num_cores)
timeEnd <- Sys.time()

