# Lick aligned Figure 4
library(data.table)
library(dplyr)

# NOTE: "sub-HJ-FP-M6_ses-Day22.nwb" file had an unequal number of reward cues and
#       rewards and so I removed the session from the folder.

# fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_new.R')
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/interval_label.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_adjust.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_RE.R') # plot random effects

path <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_4/data/"

# day when cue duration changes
mouseList <- c('HJ-FP-M2', 'HJ-FP-M3', 'HJ-FP-M4', 'HJ-FP-F1', 'HJ-FP-F2', 'HJ-FP-M6', 'HJ-FP-M7') # c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-M6','HJ-FP-M7', 'HJ-FP-F1','HJ-FP-F2')
session_max <- c(32, 27, 37, 22, 27, 22, 23) # last day of each condition
cue_change <- c(29, 24, 32, 19, 24, 19, 20) # c(29, 24, 32, 19, 20, 19, 24)
csPlus_vec <- c(15, 15, 16, 16, 15, 15, 16) # c(15, 15, 16, 15, 16, 16, 15)
cue_mat <- as.data.frame(cbind(mouseList, cue_change, csPlus_vec))
colnames(cue_mat) <- c("id", "day", "cs")

n <- length(mouseList)
preSessions <- 3 # number of days before change day
session_min <- cue_change - (preSessions + 1) # assumes this folder only has files we need
min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
pre_min_tm <- 2 # pre_lick period
post_min_tm <- 14 # post lick period
cue_period_length <- 12 # used for AIC
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
index_lick <- 5
index_bgdrw = 10 # reward
index_s1 <- 15 # cue1
index_s2 <- 16 # cue 2
index_sessionend <- 0
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
reward_period_idx <- seq(preSamps + 1, cue_period_length * target_Hz + preSamps) # for AIC

# remove NAs
dat <- dat[complete.cases(dat[,out_index]),] # remove trials with missing photometry data
dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
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
cue_period <- seq(preSamps + 1, preSamps + 3.5 * target_Hz)
cue_post <- seq(max(cue_period) + 1, preSamps + 10 * target_Hz)

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

############################
# flME Modeling
############################
# Comparable to random slope model but not selected because inference is less conservative
dat2 <- dat %>% as_tibble() %>%
            dplyr::filter(session %in% c(-1,0), # day 0 is when behavior changes
                          cs == 1) %>% # CS+ trials only
            as.data.frame()

dat2$trial <- scale(dat2$trial) # scale for interpretability

set.seed(1)
#########################################################################################
# random slope model
DA_delay <- fui(formula = photometry ~ delay + 
                    (delay | id), 
                  data = dat2,  
                  family = "gaussian", 
                  silent = FALSE,
                  var = TRUE, 
                  analytic = TRUE,
                  parallel = FALSE,
                  G_return = TRUE,
                  nknots_min = nknots_min,
                  smooth_method = "GCV.Cp",
                  splines = "tp",
                  subj_ID = "id",
                  REs = TRUE)

colMeans(DA_delay$aic[reward_period_idx,]) # do not use pre-reward period values for AIC

tm <- pre_min_tm * target_Hz # lick onset
fit_dat <- DA_delay

# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Mean Signal on Short Delay Trials", "Mean Signal Difference: Long-Short"))#,
  
  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1,0), c(0, 2)),
                                     x_interval_text = list(c(-1.25), c(1.75)), # 
                                     text_list = list("Baseline", "Cue Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "exp3_fig4g_delay_rand_slope.pdf",
        plot = fig,
        width = 8,
        height = 4)

#----------------------
# plot random effects
#----------------------
ylim = list(c(-0.25, 10.5), c(-2.25, 4.5))
# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI_RE(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               ylim = ylim[[prdtr]],
                               var_name = c("Random Intercept: Short Delay Mean", "Random Slope: Long-Short Difference"))#,
  
  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1,0), c(0, 2)),
                                     x_interval_text = list(c(-1.25), c(1.75)), # 
                                     text_list = list("Baseline", "Cue Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "exp3_fig4g_delay_rand_slope_REs.pdf",
        plot = fig,
        width = 8,
        height = 4)

rm(fig, plot.f4)

###############
# zoomed in 
###############
# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               ylim = list(NULL, c(-2,2))[[prdtr]],
                               xlim = c(-1, 3.5),
                               var_name = c("Mean Signal on Short Delay Trials", "Mean Signal Difference: Long-Short"))#,
  
  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1,0), c(0, 2)),
                                     text_list = list("Baseline", "Cue Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.001, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     ylim = list(NULL, c(-2,2.1))[[prdtr]],
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

plot.f4[[2]] <- interval_label(fig = plot.f4[[2]],
                               x_interval_list = list(c(0, 0.5), c(0.5, 2)),
                               text_list = list("(1)", "(2)"),
                               scl = 1.01, # percent for lines above original ylim values
                               x_scl = 0.001, # percent for gap in lines
                               txt_adjust = -0.278, # percent text is above line
                               txt_size = 3,
                               col_seq = c("#0868ac", "#0868ac"),
                               y_val = 1.75, 
                               alpha = 0.75)

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "exp3_fig4g_delay_rand_slope_zoom.pdf",
        plot = fig,
        width = 8,
        height = 4)

