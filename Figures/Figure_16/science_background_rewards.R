# Lick aligned Figure 4
library(data.table)
library(dplyr)
library(fastFMM)
library(lme4)

# fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_new.R')
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/interval_label.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_adjust.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_RE.R') # plot random effects

path <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/background_reward/"

# day when cue duration changes
mouseList <- c('HJ-FP-M2', 'HJ-FP-M3', 'HJ-FP-M4', 'HJ-FP-F1', 'HJ-FP-F2', 'HJ-FP-M6', 'HJ-FP-M7') # c('HJ-FP-M2','HJ-FP-M3','HJ-FP-M4','HJ-FP-M6','HJ-FP-M7', 'HJ-FP-F1','HJ-FP-F2')
session_max <- c(43, 38, 47, 27, 33, 30, 29) # last day of each condition
session_min <- c(33, 28, 38, 23, 28, 23, 24) # c(29, 24, 32, 19, 20, 19, 24)
csPlus_vec <- c(15, 15, 16, 16, 15, 15, 16) # c(15, 15, 16, 15, 16, 16, 15)
cue_mat <- as.data.frame(cbind(mouseList, NA, csPlus_vec))
colnames(cue_mat) <- c("id", "day", "cs")

n <- length(mouseList)
preSessions <- 3 # number of days before change day
min_tm <- 3 # 3 seconds before and after sucrose to be considered acceptable trial
iri_background <- 3 # remove background reward trials that are less than this many seconds apart
align <- c("cs", "b_rew", "b_rew_lick")[3] # what to align to: CS+, background reward delivery, or background reward delivery first lick
if(align == "cs"){
  pre_min_tm <- 2 # pre_delivert period
  post_min_tm <- 14 # post delivery period
  cue_period_length <- 12 # used for AIC
}else{
  pre_min_tm <- 2 # pre_lick period
  post_min_tm <- 5 # post lick period
  cue_period_length <- 2 # used for AIC
}

total_trial_length <- pre_min_tm + post_min_tm # total trial length
trial_num <- 100 # trials per session
ibi_length <- 1 # 1 second inter-lick-interval for bout calculator
trial_length <- post_min_tm + pre_min_tm # total length of session in seconds
cs_rew_diff <- c(8.95, 9.05) # CS timestamp - Reward timestamp can be within this range to be considered (should be 9seconds)

# samples before after
Hz <- 1 / 0.008
target_Hz <- 25
pre_samps <- round(pre_min_tm * Hz) # samples before lick
post_samps <- round(post_min_tm * Hz) # samples after lick

# event indices
index_lick <- 5
index_fix <- 10 # reward
back_bckrew <- 7
index_s1 <- 15 # cue1
index_s2 <- 16 # cue 2
index_sessionend <- 0
index_trial_start <- 12 # for both cs+ and cs-
# 7, 6, 14
dat_lst <- vector(length = 5 * n, "list")
cnt <- 0

for(i in 1:n){
  
  sess_vec <- seq(session_min[i], session_max[i]) # session range
  id <- cue_mat$id[i] # id name
  index_csPlus <- as.integer(cue_mat$cs[i] ) # index of cs+
  index_csMinus <- ifelse(index_csPlus == index_s1, index_s2, index_s1)
  cnt_i <- 0
  for(s in sess_vec){
    #suff <- ifelse(s == sess, "-8s.nwb", ".nwb")  
    suff <- ".nwb" # altered file *names* to remove the "-8s-bgd6s" suffix on some because was inconsistent across animals -- file content was kept identical
    
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
      ttls_rew <- ttls %>% as_tibble() %>% dplyr::filter(event == index_fix ) # fixed reward
      ttls_csMinus <- ttls %>% as_tibble() %>% dplyr::filter(event == index_csMinus ) # CS-
      ttls_bckRew <- ttls %>% as_tibble() %>% dplyr::filter(event == back_bckrew ) # CS-

      # find rewards that are background rewards
      rewards_idx_ttl <- sapply(ttls_rew$time, function(x)  min( (x - ttls_cs$time)[(x - ttls_cs$time) > 0] )   ) # shortest time between reward and a CS+ -- if any rewards occur before first CS+, R throws an Inf error but this results in correct indices (no issues)
      rew_idx <- which(rewards_idx_ttl >= cs_rew_diff[1] & rewards_idx_ttl <= cs_rew_diff[2]) # which one are the right distance from CS+ (roughly 9s)
      back_idx_ttl <- ttls_bckRew$time # background rewards
      rew_idx_ttl <- ttls_rew$time[rew_idx] # CS+ rewards
      back_rew_trials <- length(back_idx_ttl) # number of background reward trials
      
      reward_lengths <- round(rew_idx_ttl - ttls_cs$time, 2) # find time between CS+ and reward associated with CS+ (should be ~9s)
      
      # lick times
      reward_times <- back_idx_ttl # background reward times (not indices)
      bgrwd_iri <- diff(reward_times) # time between successive background rewards
      rm_idx_lick <- rm_idx <- NULL
      trial_rm <- c() # indicator of trial numbers (indices of rows) to remove #1:back_rew_trials # inclusion indicator
      trial_orig <- 1:back_rew_trials
      
      # ***** this assumes the first lick has to be before reward delivery of a subsequent reward (either background or cue-associated)
      lick_times <- ttls$time[ which(ttls$event == index_lick) ] # reward times from behavior
      lick_idx <- sapply(lick_times, function(x) which(x == dat_photo$timestamps) ) # lick indices in photometry
      
      # *** inclusion criteria in first_lick 
      first_lick_idx <- first_lick_raw <- sapply(back_idx_ttl, function(x)  lick_times[lick_times > x & 
                                                                             lick_times < min( back_idx_ttl[back_idx_ttl > x])][1] ) # time-to-first lick latency in post-reward period
      first_lick_latency <- first_lick_idx - back_idx_ttl # time-to-first lick latency in post-reward period
      trial_rm <- c(trial_rm, which(is.na(first_lick_idx))) # remove licks that do not meed criteria above
      
      # remove licks with less than 3 s apart per Huijeong email
      if(length(trial_rm) > 0){
        fl_idx <- first_lick_idx[-trial_rm]
      }else{
        fl_idx <- first_lick_idx
      }
      if(any(diff(fl_idx) < 3)){
        rm_idx_lick <- which(diff(fl_idx) < 3)
        rm_idx_lick <- (trial_orig[-trial_rm])[unique(c(rm_idx_lick, rm_idx_lick+1))] # **** removes BOTH lick before AND after lick that happens within 3s
        trial_rm <- c(trial_rm, rm_idx_lick[!is.na(rm_idx_lick)]) # update removal list
      }
      first_lick_idx <- sapply(first_lick_idx, function(x) which(x == dat_photo$timestamps))
      first_lick_idx <- sapply(first_lick_idx, function(x) ifelse(length(x) > 0, x, NA)) # make sure to have NAs in empty places
      first_lick_latency <- round(first_lick_latency, 4)  
      
      # indices of CS+/CS-
      idx_ttls_cs <- sapply(ttls_cs$time, function(x) which(dat_photo$timestamps == x))  #which(ttls$event == index_csPlus ) # CS+ or reward
      idx_ttls_csMinus <- sapply(ttls_csMinus$time, function(x) which(dat_photo$timestamps == x)) # CS+ or reward
      idx_ttls_backrwd <- sapply(back_idx_ttl, function(x) which(dat_photo$timestamps == x)) # photometry data aligned to background reward onset

      # combine data
      if(align == "cs"){
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
        sess_info <- as.data.frame( cbind(id, cnt_i - (preSessions + 1), c(1:length(reward_lengths), rep(NA, length(ttls_csMinus$time))), 
                                          idx_mat[,3:4] ) )
        colnames(sess_info) <- c("id", "session", "trial", "delay", "cs")
      }else{
        # align to background reward or first lick after background reward delivery
        ######
        trials <- seq(1, back_rew_trials )
        # remove trials before photometry starts
        min_photo_time <- first_lick_idx - pre_samps
        if(any(min_photo_time[!is.na(min_photo_time)] < 0)){
          rm_pre_photometry <- which(first_lick_idx - pre_samps < 0)
          trial_rm <- c(trial_rm, rm_pre_photometry)
        }
        if(length(trial_rm) > 0){
          # remove bad trials
          trial_rm <- sort(unique(trial_rm))
          trials <- trials[-trial_rm]
          first_lick_idx <- first_lick_idx[-trial_rm] # remove NAs
          first_lick_latency <- first_lick_latency[-trial_rm]
          first_lick_raw <- first_lick_raw[-trial_rm]
        }  
        
        # center around licks
        if(align == "b_rew_lick"){
          idx <- t( sapply(first_lick_idx, function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for first licks
        }else if (align == "b_rew"){
          idx <- t( sapply(idx_ttls_backrwd, function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for reward deliveries--note these do not remove deliveries that are not followed by "valid lick"
        }
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
        
        # new samples before/after
        pre_samps_target <- length(photo_idx_pre) # samples before lick with new sampling rate
        post_samps_target <- length(photo_idx_post) # samples after lick with new sampling rate
        
        # combine data
        sess_info <- as.data.frame( cbind(id, s, s - min(sess_vec) + 1, trials, first_lick_latency, first_lick_raw) ) # , lick_bouts[,3:4]
        colnames(sess_info) <- c("id", "session_raw", "session", "trial", "lick_latency", "lick_session_time") # , "bout_events", "bout_rate"
        
      }

      dat_lst[[cnt]] <- as.data.frame( cbind(sess_info, dat) )
      
      rm(dat, trials, sess_info, ttls, dat_photo, idx)
        
    }
  }
}

# cocnatenate results
dat <- do.call(rbind, dat_lst)

#####################################################################################

# photometry sample length
out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
L <- length(out_index)
nknots_min <- round(L/4)
cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates
preSamps <- pre_min_tm * target_Hz
reward_period_idx <- seq(preSamps + 1, cue_period_length * target_Hz + preSamps) # for AIC

# remove NAs
dat <- dat[complete.cases(dat[,out_index]),] # remove trials with missing photometry data
photo_dat <- dat[, out_index] # photometry data
dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
rm(dat_lst)
#########################################################################################
# creat AUC variables for model comparison

# AUC calculation
if(align == "b_rew_lick"){
  pre_reward_period_length <- 1.5 # length of pre-reward period (baseline)
  pre_lick_reward_period <- 0.5
  reward_period_length <- 1.5
  pre_min_tm <- pre_reward_period_length + pre_lick_reward_period # pre_lick period
  post_lick_reward_period <- 1 # 1 second after 
  post_reward_length <- 1.5 # 2 seconds after reward period ends
  post_min_tm <- post_reward_length + post_lick_reward_period # post lick period
  
  pre_lick_reward_period_samps <- round(target_Hz * pre_reward_period_length) # samples of pre_delivery period
  reward_period_max_idx <- round( (pre_reward_period_length + reward_period_length) * target_Hz ) # max_idx of post_cue period
  dat_AUC_pre <- rowMeans( photo_dat[, 1:pre_lick_reward_period_samps] ) # AUC before reward period
  dat_AUC_post <- rowMeans(photo_dat[, seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)]) # AUC after reward period
  
  reward_period_idx=seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)

  # baseline adjusted AUCs
  dat$dat_AUC_post <- dat_AUC_post
  dat$dat_AUC <- dat_AUC_post - dat_AUC_pre # baseline adjusted reward period
}

#########################################################################################
# session-wise figures
#######################################################
# OLS pooled
#######################################################
sess_range <- 1:10 # take sessions where most animals have data
sim_dat <- data.frame(id = dat$id, session = dat$session, trial = dat$trial, 
                      DA = dat$dat_AUC, 
                      x = dat$trial )

simPoint <- sim_dat[sim_dat$session %in% sess_range,] # remove other sessions
simPoint$point <- 0 # initialize vector

# save coefficients
coef_mat <- data.frame(matrix(NA, nrow = length(unique(simPoint$id)) * length(sess_range), ncol = 4))
colnames(coef_mat) <- c("session", "id", "beta0", "beta1")
cnt <- 1

for(s in sess_range){
  for(i in unique(simPoint$id)){
    idx <- which(simPoint$id == i & simPoint$session == s) # rows for animal/session
    if(length(idx) > 0){
      mm <- lm(DA ~ trial, data = simPoint[idx, ]) # fit OLS
      coef_mat[cnt,] <- c(s, i, coef(mm)) # save coefficients
      cnt <- cnt + 1

      # minimum trial in each session
      min_trial <- which(simPoint$trial[idx] == min(simPoint$trial[idx]))
      min_trial_removed <- idx[-min_trial] # all except minimum trial
      simPoint[ idx[min_trial], ]$point <- predict(mm, data.frame(trial = min(simPoint$trial[idx]))) # predict on first trial
      simPoint <- simPoint[-min_trial_removed,] # remove everything except first point
      rm(mm)
    }
  }
}

# calculate average intercept (fitted value of first trial)
simPoint_mean <- simPoint %>%
  as_tibble() %>%
  group_by(id) %>%
  summarise(m = mean(point))
simPoint <- left_join(simPoint, simPoint_mean, by = "id")

# # animal/session names
sim_dat$session_name = paste('Session', as.numeric(as.factor(sim_dat$session)))
simPoint$session_name = paste('Session', as.numeric(as.factor(simPoint$session)))

library(latex2exp)
fig <-
  sim_dat %>%
  dplyr::filter(session %in% sess_range) %>%
  ggplot(aes(x=x, y=DA), fill="white") +
  geom_point(aes(color = id),  fill="white", size=1, alpha=0.9, shape = 21) +
  facet_grid(id ~ session_name, scales = "free_y") +
  theme_classic() +
  # mixed model individual curves
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            color = "black",
            alpha = 0.75) +
  geom_point(data = simPoint, aes(x = trial, y = point),
             fill = "black", alpha = 0.75, pch = 21, size = 2) + 
  geom_hline(data = simPoint,
             aes(yintercept=m), linetype="dashed", color = "black") +
  theme(axis.text = element_text(size=10, face="bold"),
        axis.title = element_text(size=14, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),
        legend.title = element_text(face = "bold", size = 12)) +
  guides(fill = "none", alpha = "none",
         color = "none") +
  ggtitle("Reward Number OLS") +
  xlab("Trial") +
  ylab("Reward Period DA")

# # save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science Background Rewards/Figures")
ggsave( "Background_rew_AUC_per_animal_session_OLS.pdf",
        plot = fig,
        width = 8,
        height = 8)
############################
# flME Modeling
############################
dat2 <- dat 
dat2$trial_scl <- (dat2$trial - 1) / 180 # scale for interpretability
dat2$session_scl <- scale(dat2$session)
set.seed(1)
#########################################################################################
# session pooled analysis
DA_delay2 <- fastFMM::fui(formula = photometry ~ trial_scl + session + 
                            (1 | id/session), 
                          data = dat2[dat2$session <= 6,],  
                          family = "gaussian", 
                          silent = FALSE,
                          parallel = TRUE,
                          nknots_min = nknots_min,
                          subj_ID = "id")

tm <- pre_min_tm * target_Hz # lick onset
fit_dat <- DA_delay2

# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Intercept: Mean on First Trial", "Trial Number", "Session Number"))#,
  
  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1.5,0), c(0, 1.5)),
                                     x_interval_text = list(c(-0.75), c(0.75)), # 
                                     text_list = list("Baseline", "Cue Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 3))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science Background Rewards/Figures")
ggsave( "Background_Rewards_Pooled.pdf",
        plot = fig,
        width = 4,
        height = 12)

# -------------------------
# session-specific analysis
# -------------------------
for(s in 1:6){
  DA_delay2 <- fastFMM::fui(formula = photometry ~ trial_scl + 
                              (1 | id), 
                            data = dat2[dat2$session == s,],  # very interesting: [dat2$session >= 4 & dat2$session <= 6,],  
                            family = "gaussian", 
                            silent = FALSE,
                            parallel = TRUE,
                            nknots_min = nknots_min,
                            subj_ID = "id")
  
  
  tm <- pre_min_tm * target_Hz # lick onset
  fit_dat <- DA_delay2
  
  # plot
  plot.f4 <- list()
  for(prdtr in 1:nrow(fit_dat$betaHat)){
    plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = tm,
                                 var_name = c(paste0("Intercept: First Trial (Session ",s,")"), paste0("Trial Number (Session ", s, ")")))#,
    
    # add interval names
    plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                       x_interval_list = list(c(-1.5,0), c(0, 1.5)),
                                       x_interval_text = list(c(-0.75), c(0.75)), # 
                                       text_list = list("Baseline", "Cue Period"),
                                       scl = 1.01, # percent for lines above original ylim values
                                       x_scl = 0.0025, # percent for gap in lines
                                       txt_adjust = 0.03, # percent text is above line
                                       txt_size = 3,
                                       col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
  }
  
  # align plots
  fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))
  
  # save
  setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science Background Rewards/Figures")
  ggsave( paste0("Background_Rewards_Session-",s, ".pdf"),
          plot = fig,
          width = 4,
          height = 8)
  
}
