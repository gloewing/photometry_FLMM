# simulations using Science paper's data for data-driven sims
# makes figures showing individual observations as examples
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
# simParams <- expand.grid(TRUE, c(FALSE), FALSE, c(4:12), c(1) )
# colnames(simParams) <- c("fix_smooth", "resid_var_indep", "resid_var_subj", "n", "beta_mult")
# write.csv( simParams, "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photo_params.csv", row.names = FALSE)

cluserInd <- FALSE # whether running on computer or on cluster
num_cores <- 10 # parallelize

simParams <- read.csv( "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photo_params.csv" )

runNum <- 1
iter <- 2
data_path <-"/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper" # path to original Matlab files for data
source("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/photometry_sim_fLME_fn_multi.R") # simulations
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM.R')
wd <- "~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/"

# simulation parameters
fixed_smooth <- TRUE # whether to use smoothed coefficients as truth or not
reps <- 100 # replicates of simulations
boots <- 5000 # for naive bootstrap permutation test
knots_div <- 4 # number of functional samples divisor (so if it is 1, there is a knot for every obsveration on fn domain)
resid_var_indep <- FALSE # assume residual variance is not indepedent across fn. domain as in paper
resid_var_subj <- TRUE # residual variance subject specific or global
n_sim <- 7 
trial_trim <- 3 # take every other trial to reduce computational burden
target_Hz <- 15
cue_period_length <- 2
fixed_effects_scale <- 1 #simParams[runNum, 5]
random_effects_scale <- 1
resid_scale <- 5 # increase correlation across fn domain
comparison_methods <- TRUE # whether to fit other models besides fLME
rand_int <- FALSE
fixed_shift = 0
ids <- 7

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
pre_pro <- FALSE
#################################################
if( pre_pro ){
  
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
nknots_min <- round(L/2) # L = 100

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
                   subj_ID = "id",
                   num_cores = num_cores)

timeEnd <- Sys.time()

p <- nrow(fit$betaHat)

##################
# simulate data
##################
set.seed(iter)
beta_idx <- 2 # this corresponds to cue effect --only evaluate performance of models on these to make comparable with t-tests

dat_sim <- photo_stimulate(X = NULL, # design matrix
                           Z = NULL,
                           include = dat_photo["trial"],
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
  dplyr::mutate(time = as.numeric(time) / target_Hz - pre_min_tm) %>%
  dplyr::group_by(time, ID) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) )

# first plot average for animal 1 +/- s.e.m.
p6 <- mean_data %>%
  ggplot(aes(x = as.numeric(time), y = photometry, ymin=y_low, ymax=y_high, fill = as.factor(ID) )) + # ,  fill = "#ca0020"      , 
  geom_line(size = 0.75, aes(colour = as.factor(ID))) +  
  scale_color_grey(start = 0.7, end = 0.1) +
  scale_fill_grey(start = 0.7, end = 0.1) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Session Averages of 7 Simulated Animals") +
  theme_classic() + 
  geom_vline(xintercept= 0,
             linetype="dashed",
             color = "black") + #,
             # size = rel(0.3),
             # alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=16, face="bold"),
        axis.text=element_text(face="bold",color="black", size=rel(1.25)),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(face="bold", color="black", size = rel(2)), 
        legend.title = element_text(face="bold", color="black", size = rel(2)),
        strip.text.x = element_text(face="bold", color="black", size = rel(1)),
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none")


setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_IndividualAvgs_obs_delay.pdf",
        plot = p6,
        width = 8,
        height = 4)
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
  dplyr::filter(ID == ids ) %>%  # arbitrarily choose an example animal
  dplyr::mutate(time = as.numeric(time) / target_Hz - pre_min_tm) %>%
  dplyr::group_by(time, cue, trial) 

mean_data$cue <- ifelse(mean_data$cue == 0, "Short", "Long")

# average over cue0 trials
mean_dat0 <- mean_data  %>%
  group_by(cue, time) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) # +/- sem
  )

# filter and choose specific trials
set.seed(1)
trials_unique <- unique(mean_data$trial)
trials_unique <- trials_unique[sample.int(length(trials_unique), 8)]
mean_dat2 <- mean_data %>% 
  group_by(cue, time, trial) %>%
  dplyr::filter(trial %in% trials_unique,
                ID == ids) %>% # subset of rows # ,
  dplyr::rename(photometry = value)

mean_dat2$cue <- as.factor(mean_dat2$cue)
mean_dat2$trial[mean_dat2$cue == "Short"] <- mean_dat2$trial[mean_dat2$cue == "Short"] + 10
mean_dat2$trial <- as.factor(mean_dat2$trial)

# first plot average for animal 1 +/- s.e.m.
p4 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +  
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 1) +  
  geom_line(data = mean_dat2, aes(colour = as.factor(cue), group = as.factor(trial)), size = 0.5, alpha = 0.5) +  #, colour = "#ca0020"
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Simulated Trial-Level Data") +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  geom_vline(xintercept= 0,
             linetype="dashed",
             color = "black") + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=16, face="bold"),
        axis.text=element_text(face="bold",color="black", size=rel(1.25)),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(face="bold", color="black", size = rel(1.25)), 
        legend.title = element_text(face="bold", color="black", size = rel(1.25)),
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) + 
  guides(color= guide_legend(title="Delay"), fill = "none")

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_sims_science_IndividualTrials_obs_delay.pdf",
          plot = p4,
          width = 8,
          height = 4)

#####################################################
# original data
mean_data <- dat_photo %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::mutate(time = as.numeric(time) / target_Hz - pre_min_tm) %>%
  dplyr::group_by(time, cue, trial) 

mean_data$cue <- ifelse(mean_data$cue == 0, "Short", "Long")

# average over cue0 trials
mean_dat0 <- mean_data  %>%
  # average across trials within animal
  group_by(cue, time, id) %>%
  dplyr::summarise(photo_mean = mean(value, na.rm = TRUE)) %>% 
  # average across animals so se is uncertainty w.r.t animal averages
  group_by(cue, time) %>%
  dplyr::summarise(photometry = mean(photo_mean, na.rm = TRUE),
                   y_low = mean(photo_mean, na.rm = TRUE) - sd(photo_mean) / sqrt(n()), # +/- sem
                   y_high = mean(photo_mean, na.rm = TRUE) + sd(photo_mean) / sqrt(n()) )

# uncertainty if trials are independent
# mean_data  %>%
#   group_by(cue, time) %>%
#   dplyr::summarise(photometry = mean(value, na.rm = TRUE),
#                    y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
#                    y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) ) # +/- sem
#   


p1 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +      
  geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, fill = as.factor(cue)), alpha=0.5) + 
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 1) +  
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Original Data Session Averages") +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  # geom_vline(xintercept= 0,
  #            linetype="dashed",
  #            color = "black",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  # geom_vline(xintercept= 2,
  #            linetype="dashed",
  #            color = "grey 20",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  # geom_vline(xintercept= 2.5,
  #            linetype="dashed",
  #            color = "grey 40",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  # geom_vline(xintercept= 3,
  #            linetype="dashed",
  #            color = "grey 60",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=16, face="bold"),
        axis.text=element_text(face="bold",color="black", size=rel(1.25)),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(face="bold", color="black", size = rel(1.25)), 
        legend.title = element_text(face="bold", color="black", size = rel(1.25)),
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) + 
  guides(color= guide_legend(title="Delay"), fill = "none") + 
  coord_cartesian(ylim = c(-0.25, 7.5)) +
  # bars over top
  geom_segment(aes(x=0,xend=2,y=7,yend=7), color = "grey 20", size = 1) +
  geom_segment(aes(x=0,xend=2.5,y=7.1,yend=7.1), color = "grey 40", size = 1) +
  geom_segment(aes(x=0,xend=3,y=7.2,yend=7.2), color = "grey 60", size = 1) +
  # vertical dotted lines
  geom_segment(aes(x=0,xend=0,y=0,yend=7.2), color = "black", linetype="dashed") +
  #geom_segment(aes(x=2,xend=2,y=0,yend=6), color = "grey 20", linetype="dashed") +
  #geom_segment(aes(x=2.5,xend=2.5,y=0,yend=6.1), color = "grey 40", linetype="dashed") +
  #geom_segment(aes(x=3,xend=3,y=0,yend=6.2), color = "grey 60", linetype="dashed") + 
  annotate(geom="text", x=1.5, y=7.4, label="Extended Cue Periods", color = "black")
  
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_science_original_data_extended2.pdf",
        plot = p1,
        width = 8,
        height = 4)


##########################################################################################################

p1_2 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +   
  geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, fill = as.factor(cue)), alpha=0.5) + 
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 1) +  
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Original Data Session Averages") +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  # geom_vline(xintercept= 0,
  #            linetype="dashed",
  #            color = "black",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  # geom_vline(xintercept= 2,
  #            linetype="dashed",
  #            color = "grey 20",
  #            size = rel(0.3),
  #            alpha = 0.95) +
  # geom_vline(xintercept= 2.5,
#            linetype="dashed",
#            color = "grey 40",
#            size = rel(0.3),
#            alpha = 0.95) +
# geom_vline(xintercept= 3,
#            linetype="dashed",
#            color = "grey 60",
#            size = rel(0.3),
#            alpha = 0.95) +
geom_hline(yintercept=0, 
           linetype="dashed", 
           color = "black", 
           size = rel(0.5),
           alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=16, face="bold"),
        axis.text=element_text(face="bold",color="black", size=rel(1.25)),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(face="bold", color="black", size = rel(1.25)), 
        legend.title = element_text(face="bold", color="black", size = rel(1.25)),
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) + 
  guides(color= guide_legend(title="Delay"), fill = "none") + 
  coord_cartesian(ylim = c(-0.25, 7.75)) +
  # bars over top
  geom_segment(aes(x=0,xend=2,y=7,yend=7), color = "grey 20", size = 1) +
  geom_segment(aes(x=0,xend=2.5,y=7.4,yend=7.4), color = "grey 40", size = 1) +
  geom_segment(aes(x=0,xend=3,y=7.5,yend=7.5), color = "grey 60", size = 1) +
  # vertical dotted lines
  geom_segment(aes(x=0,xend=0,y=0,yend=7.5), color = "black", linetype="dashed") +
  # geom_segment(aes(x=2,xend=2,y=0,yend=7), color = "grey 20", linetype="dashed") +
  # geom_segment(aes(x=2.5,xend=2.5,y=0,yend=7.4), color = "grey 40", linetype="dashed") +
  # geom_segment(aes(x=3,xend=3,y=0,yend=7.5), color = "grey 60", linetype="dashed") + 
  annotate(geom="text", x=1.5, y=7.7, label="Extended Cue Periods", color = "black") +
  annotate(geom="text", x=1, y=7.2, label="Cue Period", color = "black")

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_science_original_data_extended.pdf",
        plot = p1_2,
        width = 8,
        height = 4)


##########################################################################################################

p2 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +      
  geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, fill = as.factor(cue)), alpha=0.5) + 
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 1) +  
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Original Data Session Averages") +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
geom_hline(yintercept=0, 
           linetype="dashed", 
           color = "black", 
           size = rel(0.5),
           alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=16, face="bold"),
        axis.text=element_text(face="bold",color="black", size=rel(1.25)),
        legend.key.size = unit(2, "line"), 
        legend.text = element_text(face="bold", color="black", size = rel(1.25)), 
        legend.title = element_text(face="bold", color="black", size = rel(1.25)),
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) + 
  guides(color= guide_legend(title="Delay"), fill = "none") + 
  coord_cartesian(ylim = c(-0.25, 7.25)) +
  # bars over top
  geom_segment(aes(x=0,xend=2,y=7,yend=7), color = "grey 20", size = 1) +
  # vertical dotted lines
  geom_segment(aes(x=0,xend=0,y=0,yend=7), color = "black", linetype="dashed") +
  # geom_segment(aes(x=2,xend=2,y=0,yend=6), color = "grey 20", linetype="dashed") +
  annotate(geom="text", x=1, y=7.2, label="Cue Period", color = "black")

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/sims/Final")
ggsave( "photo_science_original_data.pdf",
        plot = p2,
        width = 8,
        height = 4)


##########################################################################################################
# combine

# library(patchwork)
# set.seed(25)
# plot_combined <- p1 + p4 + p6 +
#   plot_layout(ncol = 2,
#               nrow = 2,
#               byrow = NULL,
#               widths = NULL,
#               heights = 10,
#               guides = NULL)
# 

