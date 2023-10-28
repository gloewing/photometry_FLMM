# Lick aligned Experiment 1 
# Simpsons Paradox
library(data.table)
library(dplyr)
library(parallel)
library(lme4)
library(mgcv)

# fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_new.R')
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/interval_label.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_adjust.R') # plot

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
short_iri_adjust <- TRUE # if true remove both trials before and after. if false remove only after short iris

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
      iri <- round(iri, 4)
      
      # variance of mean of 100 trials is (100 / (1/12)^2 ) / 100^2 = 1.44, so 8 seconds is > 3 * sd less than mean
      if(iri_length > 5){
        
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
        if(short_iri_adjust){
          trial_iri <- unique(c(trial_iri, trial_iri + 1)) # remove trials before and after short iris because their pre-peroid is not long enough
        }else{
          trial_iri <-  which(diff(reward_times) < min_tm) # which trials to remove
          trial_iri <- unique(trial_iri + 1) # remove trials only after short iris because their pre-peroid is not long enough
        }
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
        if(short_iri_adjust){
          lick_bouts <- bout_finder_trial(x = lick_times, Hz = Hz, ibi = ibi_length, bout_start_time = first_lick_session, include = trial_iri)
        }else{
          lick_bouts <- matrix(NA, nrow = length(trial_iri), ncol = 4)
        }
        
        # first lick index
        first_lick_idx <- sapply(reward_times, function(x)  lick_times[I(lick_times > x & lick_times <= (x + min_tm) )][1] ) # number of licks in post-reward period
        first_lick_idx <- sapply(first_lick_idx, function(x) which(x == dat_photo$timestamps))
        first_lick_idx <- sapply(first_lick_idx, function(x) ifelse(length(x) > 0, x, NA)) # make sure to have NAs in empty places
        first_lick <- round(first_lick, 4)  
        
        # center around licks
        idx <- t( sapply(first_lick_idx[trial_iri], function(x) seq(x - pre_samps, x + post_samps) ) ) # matrix of indices for trials
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
if(short_iri_adjust)    dat <- dat[!is.na(dat$bout_rate),] # remove trials with problematic bouts (1.4% of trials)
dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
dat2 <- dat
rm(dat_lst)

#####################################################################################
# photometry sample length
out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
L <- length(out_index)
nknots_min <- round(L/2)
cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates

# create AUC variables 
photo_dat <- dat[, -cov_idx] # photometry data
pre_lick_reward_period_samps <- round(target_Hz * pre_reward_period_length) # samples of pre_cue period
reward_period_max_idx <- round( (pre_reward_period_length + reward_period_length) * target_Hz ) # max_idx of post_cue period
dat_AUC_pre <- rowMeans( photo_dat[, 1:pre_lick_reward_period_samps] ) # AUC before reward period
dat_AUC_post <- rowMeans(photo_dat[, seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)]) # AUC after reward period

reward_period_idx=seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)

# scale covariates to make contrasts more interpretable
dat$reward_num <- dat$trial + (dat$session - 1) * num_trials # number of rewards assuming 100 rewards/session
dat$session_orig <- dat$session
dat$session <- scale(dat$session, center = TRUE, scale = FALSE) # consider mroe interpertable version
dat$trial_orig <- dat$trial 
dat$trial <- dat$trial - trial_num / 2 # trial 50 is middle
# dat$licks <- scale(dat$licks, center = TRUE, scale = FALSE)
# dat$lick_time <- scale(dat$lick_time, center = TRUE, scale = FALSE)
dat$iri <- scale(dat$iri, center = TRUE, scale = FALSE)
dat$dat_AUC_pre <- dat_AUC_pre
dat$dat_AUC_post <- dat_AUC_post
dat$dat_AUC <- dat_AUC_post - dat_AUC_pre # described in "Experiment 1" of "Data Analysis" in Supplement
tm <- pre_min_tm * target_Hz # lick onset
dat$dat_AUC <- photo_dat[, reward_period_idx[10]] - dat_AUC_pre

if(short_iri_adjust){
  dat <- dat[!is.infinite(dat$bout_rate),] # 1 trial has infinite bout rate
  dat$bout_rate <- scale(dat$bout_rate, center = TRUE, scale = FALSE)
  dat <- dat[complete.cases(dat),] # remove one observation
}else{
  dat <- dat[complete.cases(subset(dat, select = -c(bout_events,bout_rate) )),] # remove one observation
  
}    
dim(dat)

#######################################################
# OLS pooled
#######################################################
sess_range <- 1:6 # take sessions where most animals have data
sim_dat <- data.frame(id = dat$id, session = dat$session_orig, trial = dat$trial_orig, 
                      DA = dat$dat_AUC, iri = dat$iri, rn = dat$reward_num,
                      x = dat$trial_orig )

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
sim_dat$id = paste('Animal', as.numeric(as.factor(sim_dat$id)))
sim_dat$session_name = paste('Session', as.numeric(as.factor(sim_dat$session)))
simPoint$id = paste('Animal', as.numeric(as.factor(simPoint$id)))
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
            #linetype = 5,
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
  ylab("Reward Period AUC") 

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Simpsons_Fig")
ggsave( "Simpsons_OLS.pdf",
        plot = fig,
        width = 8,
        height = 8)
#######################################################
# simpsons OLS dual lines
sim_dat2 <- sim_dat
sim_dat2$x <- sim_dat2$rn # use reward number as x

simPoint2 <- simPoint
simPoint2$rn <- (simPoint2$session - 1) * trial_num + simPoint2$trial
library(latex2exp)
fig <- 
  sim_dat2 %>%
  dplyr::filter(session %in% sess_range) %>%
  ggplot(aes(x=x, y=DA), fill="white") +
  geom_point(aes(color = session_name),  fill="white", size=1, alpha=0.75, shape = 21) +
  facet_grid(id ~., scales = "free_y") +
  theme_classic() + 
  # mixed model individual curves
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1.25,
            #linetype = 5,
            color = "black",
            alpha = 1) +
  geom_line(data = sim_dat2[sim_dat2$session %in% sess_range,], 
            stat="smooth", method = "lm", formula = y ~ x,
            size = 0.95,
            aes(color = session_name, group = session_name)) +
  geom_point(data = simPoint2, aes(x = rn, y = point),
             fill = "black", alpha = 0.85, pch = 21, size = 2) + 
  scale_fill_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252") ) +
  scale_color_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252")) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(1,600), breaks = c(1,200,400,600)) +
  geom_hline(data = simPoint, 
             aes(yintercept=m), linetype="dashed", color = "black") +
  theme(axis.text = element_text(size=10, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"), 
        legend.position="bottom", 
        legend.title = element_blank()) +
  guides(fill = "none", alpha = "none") +
  ggtitle("Dopamine - Reward Number Association by Session") +
  xlab("Reward Number") +
  ylab("Reward Period AUC") 

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Simpsons_Fig")
ggsave( "Simpsons_Dual.pdf",
        plot = fig,
        width = 6,
        height = 8)

########################
# Simpsons Dual subset
########################
# subset first 3 animals
id_subset <- 1:3
id_subset <- paste("Animal",id_subset)
sim_dat2 <- sim_dat2 %>%
  dplyr::filter(session %in% sess_range,
                id %in% id_subset)

simPoint2 <- simPoint2 %>%
  dplyr::filter(session %in% sess_range,
                id %in% id_subset)

simPoint <- simPoint %>%
  dplyr::filter(session %in% sess_range,
                id %in% id_subset)

fig <- 
  sim_dat2 %>%
  ggplot(aes(x=x, y=DA), fill="white") +
  geom_point(aes(color = session_name),  fill="white", size=1, alpha=0.75, shape = 21) +
  facet_grid(id ~., scales = "free_y") +
  theme_classic() + 
  # mixed model individual curves
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1.25,
            color = "black",
            alpha = 1) +
  geom_line(data = sim_dat2, 
            stat="smooth", method = "lm", formula = y ~ x,
            size = 0.95,
            aes(color = session_name, group = session_name)) +
  geom_point(data = simPoint2, aes(x = rn, y = point),
             fill = "black", alpha = 0.85, pch = 21, size = 2) + 
  scale_fill_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252") ) +
  scale_color_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252")) +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(1,600), breaks = c(1,200,400,600)) +
  geom_hline(data = simPoint, 
             aes(yintercept=m), linetype="dashed", color = "black") +
  theme(axis.text = element_text(size=10, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"), 
        legend.position="bottom", 
        legend.title = element_blank()) +
  guides(fill = "none", alpha = "none") +
  ggtitle("Dopamine - Reward Number Association by Session") +
  xlab("Reward Number") +
  ylab("Reward Period AUC") 

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Simpsons_Fig")
ggsave( "Simpsons_Dual_subset.pdf",
        plot = fig,
        width = 6,
        height = 6)

###################
# coefficients
###################
coef_mat <- coef_mat[complete.cases(coef_mat),] # remove NAs
coef_mat$beta0 <- as.numeric(coef_mat$beta0)
coef_mat$beta1 <- as.numeric(coef_mat$beta1)
coef_mat$session <- as.numeric(coef_mat$session)
coef_mat$id = paste('Animal', as.numeric(as.factor(coef_mat$id)))

####################
# beta 0 vs. session
####################
# inference
mod <- nlme::lme(beta0 ~ session, 
                 random = ~session | id, 
                 data = coef_mat,
                 control = lmeControl(maxIter = 5000, msMaxIter = 5000, niterEM = 5000) )
mod_coef <- summary(mod)$tTable[,1]
mod <- summary(mod)$tTable[2,c(1,5)]
coef_mat$Animal <- as.factor(as.integer(as.factor(coef_mat$id)))

fig1 <- 
  coef_mat %>%
  dplyr::filter(session %in% sess_range) %>%
  ggplot(aes(x=session, y=beta0, fill = Animal)) +
  theme_classic() + 
  # mixed model individual curves
  geom_abline(slope = mod_coef[2], intercept = mod_coef[1], # mixed model fit
              color = "black", alpha = 0.75, size = 2) +
  geom_point(aes(fill = Animal), size=2, alpha=0.9, shape = 21) + #  fill="white", 
  scale_fill_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252", "lightgray", "yellow") ) +
  scale_color_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252", "lightgray", "yellow")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(axis.text = element_text(size=10, face="bold"), 
       axis.title = element_text(size=14, face="bold"),
       plot.title = element_text(hjust = 0.5, face="bold"), 
       legend.position="bottom", 
       legend.title = element_text(size=14, face="bold")) +
  ggtitle("Intercepts Increase Across Sessions") +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(1,6), breaks = 1:6) +
  guides(alpha = "none") + 
  xlab("Session") +
  ylab(TeX('Animal/Session Intercepts:  $\\beta_0$')) +
  annotate(geom = "text", x = 5, y = 6, 
           label = TeX("p \\approx 0.0243"),
           color = "black", size = 4.25)

rm(mod)

####################
# beta 1 vs. session
####################
# inference
mod <- nlme::lme(beta1 ~ session, 
                 random = ~session | id, data = coef_mat,
                 control = nlmeControl(maxIter = 5000, niterEM = 5000))
mod_coef <- summary(mod)$tTable[,1]
mod <- summary(mod)$tTable[2,c(1,5)]
fig2 <- 
  coef_mat %>%
  dplyr::filter(session %in% sess_range) %>%
  ggplot(aes(x=session, y=beta1), fill="white") +
  theme_classic() + 
  # mixed model individual curves
  geom_abline(slope = mod_coef[2], intercept = mod_coef[1], # mixed model fit
              color = "black", alpha = 0.75, size = 2) +
  geom_point(aes(fill = Animal), size=2, alpha=0.9, shape = 21) + 
  scale_fill_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252", "lightgray", "yellow") ) +
  scale_color_manual(values = c("#ca0020", "darkgreen", "#0868ac",  "darkgray", "#E69F00", "#525252", "lightgray", "yellow")) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  theme(axis.text = element_text(size=10, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"), 
        legend.position="bottom", 
        legend.title = element_text(size=14, face="bold")) +
  ggtitle("Few Slopes are Significanty Positive") +
  scale_x_continuous(expand = c(0.02, 0.02), limits = c(1,6), breaks = 1:6) +
  guides(alpha = "none") + 
  xlab("Session") +
  ylab(TeX('Animal/Session Slopes:  $\\beta_1$')) +
  annotate(geom = "text", x = 5, y = -0.03, 
           label = TeX("p \\approx 0.804"),
           color = "black", size = 4.25)

rm(mod)

# combine beta0 and beta1 and plot
library(patchwork)
set.seed(25)
plot_combined <- fig1 + fig2 + 
  plot_layout(ncol = 2,
              nrow = 1,
              byrow = NULL,
              widths = NULL,
              heights = 6,
              guides = 'collect') &
  theme(legend.position='bottom') 

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Simpsons_Fig")
ggsave( "Simpsons_OLS_beta_fit_combined.pdf",
        plot = plot_combined,
        width = 8,
        height = 4)
