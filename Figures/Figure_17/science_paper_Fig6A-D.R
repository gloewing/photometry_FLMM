# cue aligned Figure 6
library(data.table)
library(dplyr)

# NOTE: days for each animal taken from:
# https://github.com/namboodirilab/ANCCR/blob/master/analysis/fig6/backpropagation_pavlovian.m
# and https://docs.google.com/spreadsheets/d/1pmpQ5JFhg4Q7h18DQYifjNrW17HtaoPHBt2vJF4mxiU/edit#gid=0

# Note: sub-HJ-FP-F2_ses-Day13.nwb had different numbers of rewards and CS+s so I removed the file

# # fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_new.R')
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/interval_label.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_adjust.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_RE.R') # plot random effects


path <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_6/data/"

# day when cue duration changes
mouseList <- c('HJ-FP-M2', 'HJ-FP-M3', 'HJ-FP-M4', 'HJ-FP-F1', 'HJ-FP-F2', 'HJ-FP-M6', 'HJ-FP-M7') # c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-M6','HJ-FP-M7', 'HJ-FP-F1','HJ-FP-F2')
cue_change <- c(29, 24, 32, 19, 24, 19, 20) # c(29, 24, 32, 19, 20, 19, 24)
learned_day <- c(22,17,25,12,17,12,14) # when they learned to lick (only needed for name)
session_min <- c(16,12,13,7,13,7,6) # first day for each animal (from matlab script)
session_max <- cue_change - 1 # analyze all days before change in cue length
csPlus_vec <- c(15, 15, 16, 16, 15, 15, 16) # c(15, 15, 16, 15, 16, 16, 15)
cue_mat <- as.data.frame(cbind(mouseList, learned_day, csPlus_vec))
colnames(cue_mat) <- c("id", "day", "cs")

n <- length(mouseList)
min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
pre_min_tm <- 1 # pre_lick period
post_min_tm <- 6 # post lick period
cue_period_length <- 6 # used for AIC
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
index_bgdrw = 10 # reward
index_s1 <- 15 # cue1
index_s2 <- 16 # cue 2
index_sessionend = 0
index_trial_start <- 12 # for both cs+ and cs-

dat_lst <- vector(length = 5 * n, "list")
cnt <- 0

for(i in 1:n){
  
  sess <- as.integer(cue_mat$day[i]) # change of learning
  sess_vec <- seq(session_min[i], session_max[i])
  id <- cue_mat$id[i] # id name
  index_csPlus <- as.integer(cue_mat$cs[i] ) # index of cs+
  index_csMinus <- ifelse(index_csPlus == index_s1, index_s2, index_s1)
  cnt_i <- 0
  for(s in sess_vec){
    suff <- ifelse(s == sess, "-learned.nwb", ".nwb")  
    
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
      idx_ttls_cs <- sapply(ttls_cs$time, function(x) which(dat_photo$timestamps == x))  #
      idx_ttls_csMinus <- sapply(ttls_csMinus$time, function(x) which(dat_photo$timestamps == x)) # CS+ or reward
      
      
      idx_mat <- data.frame( idx = c(idx_ttls_cs, idx_ttls_csMinus), # indices
                             onset = c(ttls_cs$time, ttls_csMinus$time), # time when cue starts
                             reward_delay = c(reward_lengths, rep(NA, length(idx_ttls_csMinus) )), # delay to reward
                             csPlus = c(rep(1, length(reward_lengths)), rep(0, length(idx_ttls_csMinus) )) ) # cs+/cs-

      reward_times <- ttls_rew$time
      
      # lick times
      lick_times <- ttls$time[ which(ttls$event == index_lick) ] # reward times from behavior
      lick_idx <- sapply(lick_times, function(x) which(x == dat_photo$timestamps) ) # lick indices in photometry
      tot_licks <- sapply(reward_times, function(x)  sum(lick_times > x & lick_times <= (x + min_tm) ) ) # number of licks in post-reward period
      first_lick <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1] - x ) # time-to-first lick latency in post-reward period

      # consummatory licks
      first_lick_session <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1]) # session time of first lick in post-reward period in overall session time
      lick_bouts <- NA # 
      
      # first lick index
      first_lick_idx <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1] ) # number of licks in post-reward period
      first_lick_idx <- sapply(first_lick_idx, function(x) which(x == dat_photo$timestamps))
      first_lick_idx <- sapply(first_lick_idx, function(x) ifelse(length(x) > 0, x, NA)) # make sure to have NAs in empty places
      first_lick <- round(first_lick, 4)  
      
      # center around start of CS+
      idx <- t( sapply(idx_mat$idx, function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for trials based on CS start times
      dat <- as.data.frame( apply(idx, 2, function(x) dat_photo$data[x]) ) # photometry trials
      L <- ncol(dat)
      
      # downsample photometry
      by_num <- round( Hz / target_Hz )
      photo_idx_pre <- sort(-seq(-pre_samps, -1, by = by_num)) #
      photo_idx_post <- seq(pre_samps + by_num, L, by = by_num)
      photo_idx <- unique( c(photo_idx_pre, photo_idx_post) )
      dat <- dat[,photo_idx]
      L <- ncol(dat)
      colnames(dat) <- paste0("photometry.", 1:L)
      
      # combine data
      idx_mat[,3] <- seq(1, num_trials )
      sess_info <- as.data.frame( cbind(id, cnt_i, s, idx_mat[idx_mat$csPlus==1,3:4], tot_licks, first_lick ) ) # ,
      colnames(sess_info) <- c("id", "session", "day", "trial", "cs", "licks", "lick_time") # 
      
      dat_lst[[cnt]] <- as.data.frame( cbind(sess_info, dat[idx_mat$csPlus==1,]) )
      
      rm(dat, sess_info, ttls, dat_photo, idx)
        
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
dat$reward_num <- dat$trial + (dat$session - 1) * num_trials # number of rewards assuming 100 rewards/session

# remove NAs
dat <- dat[complete.cases(dat[,out_index]),] # remove trials with missing photometry data (169) in total
dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
dat2 <- dat
rm(dat_lst)
#########################################################################################
train_mat <- matrix(0, nrow = nrow(dat), ncol = 1)
colnames(train_mat) <- c("period") # periods

for(ii in 1:length(mouseList)){
  idx <- which(dat$id == mouseList[ii])
  
  max_sess <- max(dat$session[idx])
  
  # factor version of periods
  train_mat[idx,1] <- ifelse(dat$session[idx] <= 1 & 
                               dat$day[idx] < learned_day[ii], 1, train_mat[idx,1]) # if first 3 sessions and doesn't overlap with learned day
  
  train_mat[idx,1] <- ifelse(dat$session[idx] == 2, 2, train_mat[idx,1]) # 3-4
  train_mat[idx,1] <- ifelse(dat$session[idx] == 3, 3, train_mat[idx,1]) # 6 -7
  train_mat[idx,1] <- ifelse(dat$session[idx] > max_sess - 1, 4, train_mat[idx,1]) # max
}

# only save early/mid/later periods
train_mat <- as.data.frame(train_mat)
dat2 <- data.frame(period = as.factor(train_mat$period), dat)
dat2 <- dat2[dat2$period != 0,]
dat2$period <- as.factor(dat2$period)

# down-sample trials for computational reasons
dat2 <- dat2 %>% as_tibble() %>%
  dplyr::filter(cs == 1 # 
                ) %>% # CS+ trials only
  as.data.frame()

############################
# flME Modeling
############################
set.seed(1)
DA_early <- fui(formula = photometry ~ period +
                  (period | id), 
                  data = dat2,  
                  family = "gaussian", 
                  silent = FALSE,
                  var = TRUE, 
                  analytic = TRUE,
                  parallel = TRUE,
                  G_return = TRUE,
                  nknots_min = nknots_min,
                  smooth_method = "GCV.Cp",
                  splines = "tp",
                  subj_ID = "id",
                  REs = TRUE)
colMeans(DA_early$aic[reward_period_idx,]) # do not use pre-reward period values for AIC

tm <- pre_min_tm * target_Hz # lick onset
fit_dat <- DA_early

# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Mean Signal: 1st Session", 
                                            "Mean Change: 2nd-1st Session", 
                                            "Mean Change: 3rd-1st Session",
                                            "Mean Change: Last-1st Session"),
                               y_scal_orig = 0.01) #, 

  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(0,1), c(1, 2), c(2,3)),
                                     x_interval_text = NULL,#
                                     text_list = list("Early", "Mid", "Late"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )

  }

fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))
# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-fig6A")
ggsave( "fig6a_backprop_earlyMidLate.pdf",
        plot = fig,
        width = 8,
        height = 8)

#------------------------------------------------------------------
#----------------------
# plot random effects
#----------------------
ylim = list(c(-1, 16), c(-10, 5), c(-12, 7.5), c(-15, 12))
# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI_RE(r = prdtr, 
                                  Hz = target_Hz, 
                                  align = tm,
                                  ylim = ylim[[prdtr]],
                                  var_name = c("Mean Signal: 1st Session", 
                                               "Mean Change: 2nd-1st Session", 
                                               "Mean Change: 3rd-1st Session",
                                               "Mean Change: Last-1st Session"),
                                  y_scal_orig = 0.01) #, "Session", "Bout Rate x Session"
  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(0,1), c(1, 2), c(2,3)),
                                     x_interval_text = NULL,#
                                     text_list = list("Early", "Mid", "Late"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))
# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-fig6A")
ggsave( "fig6a_backprop_earlyMidLate_REs.pdf",
        plot = fig,
        width = 8,
        height = 8)
