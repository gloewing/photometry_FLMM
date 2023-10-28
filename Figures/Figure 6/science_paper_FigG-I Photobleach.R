# Lick aligned Figure 4--check whether activity falls with trial-- Photobleahing?
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
cue_len <- 3 # 3 seconds between cue and reward delivery: for anticipatory lick counting (Vijay request)
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
      CS_ttls <- c(ttls_cs$time, ttls_csMinus$time) # must be in order: c( CS+, CS-  ) for ordering below
      
      # save "reward" times that are closest to cs+
      rew_idx_ttl <- lapply(ttls_cs$time, function(x) which( (ttls_rew$time - x) > 0  )[1] )
      rew_idx_ttl <- do.call(c, rew_idx_ttl )
      ttls_rew <- ttls_rew[unique(rew_idx_ttl), ]
      
      # lick times
      lick_times <- ttls$time[ which(ttls$event == index_lick) ] # reward times from behavior
      lick_idx <- sapply(lick_times, function(x) which(x == dat_photo$timestamps) ) # lick indices in photometry
      tot_licks <- sapply(CS_ttls, function(x)  sum(lick_times > x & lick_times < (x + cue_len) ) ) # number of anticipatory licks
      first_lick <- sapply(CS_ttls, function(x)  lick_times[I(lick_times > x & lick_times < (x + cue_len) )][1] - x ) # time-to-first lick latency in post-cue period
      
      
      # ensure number of rewards == number of CS+
      if(length(ttls_rew$time) != length(ttls_cs$time)){
        print(paste("ERROR: i", i, "cnt", cnt_i))
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
      trials <- order(order(c(idx_ttls_cs, idx_ttls_csMinus))) #seq(1, num_trials )
      sess_info <- as.data.frame( cbind(id, 
                                        cnt_i - (preSessions + 1), 
                                        trials,#c(1:length(reward_lengths), rep(NA, length(ttls_csMinus$time))), 
                                        idx_mat[,3:4],
                                        tot_licks, 
                                        first_lick
                                        ) )
      colnames(sess_info) <- c("id", "session", "trial", "delay", "cs", "total_licks", "first_lick")
      
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

############################
# flME Modeling
############################
# Comparable to random slope model but not selected because inference is less conservative
# use only 
dat$delay <- ifelse(dat$session >= 0, 1, 0) # 1 means it is long delay session
dat2 <- dat %>% as_tibble() %>%
            dplyr::filter(delay == 0, # use only short sessions
                          cs == 1) %>% # CS+ trials only
            as.data.frame()
dat2$session <- as.factor(dat2$session)

# average number of trials per session for each animal and session
num_trial <- dat2 %>% as_tibble() %>%
            dplyr::group_by(session, id) %>%
            dplyr::mutate(maxTrial = max(trial)) %>%
            summarise(m = mean(maxTrial))

dat2 <- left_join(dat2, num_trial, by = c("id", "session"))
dat2$trial <- (dat2$trial - 1) / (dat2$m-1) # standardize trial to be [0,1]

set.seed(1)

#########################################################################################
# best model
DA_delay <- fastFMM::fui(formula = photometry ~ trial + 
                           (1 | id), 
                         data = dat2,  
                         parallel = FALSE,
                         G_return = TRUE,
                         nknots_min = nknots_min,
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
                               var_name = c("Intercept: Mean Signal on 1st Trial", "Trial Number"))#,
  
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
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "Short_Delay_Trial_Effect_LastSesson.pdf",
        plot = fig,
        width = 4,
        height = 8)

#########################################################################################
# best model
DA_delay <- fastFMM::fui(formula = photometry ~ trial + total_licks +
                           (1 | id), 
                         data = dat2,  
                         parallel = FALSE,
                         G_return = TRUE,
                         nknots_min = nknots_min,
                         subj_ID = "id",
                         REs = TRUE)

colMeans(DA_delay$aic[reward_period_idx,]) # do not use pre-reward period values for AIC

tm <- pre_min_tm * target_Hz # lick onset
fit_dat <- DA_delay
# AIC      BIC     cAIC 
# 870.4217 889.6539       NA 

# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Intercept: Mean Signal on 1st Trial", "Trial Number", "Total Licks"))#,
  
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
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "Short_Delay_Trial_Effect_LastSesson_withLicks.pdf",
        plot = fig,
        width = 8,
        height = 8)


#######################################################
# OLS pooled
#######################################################
dat2 <- dat %>% as_tibble() %>%
  dplyr::filter(delay == 0, # use only short sessions
                cs == 1 # CS+ trials only
  ) %>% 
  as.data.frame()
dat2$session <- as.factor(dat2$session)

# average number of trials per session for each animal and session
num_trial <- dat2 %>% as_tibble() %>%
  dplyr::group_by(session, id) %>%
  dplyr::mutate(maxTrial = max(trial)) %>%
  summarise(m = mean(maxTrial))

dat2 <- left_join(dat2, num_trial, by = c("id", "session"))
#dat2$trial <- (dat2$trial - 1) / (dat2$m-1) # standardize trial to be [0,1]

set.seed(1)

sess_range <- seq(-3,-1) # take sessions where most animals have data
sim_dat <- data.frame(id = dat2$id, session = dat2$session, trial = dat2$trial, 
                      DA = dat2$total_licks, 
                      x = dat2$trial )

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
             fill = "black", alpha = 0.75, pch = 21, size = 2) + # aes(fill = id), color = "black",
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
  ggtitle("Short Delay Anticipatory Licks") +
  xlab("Trial") +
  ylab("Total Anticipatory Licks")

# # save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "Delay_Length_Licks_session_OLS.pdf",
        plot = fig,
        width = 8,
        height = 8)


num_subj <- length(unique(dat2$id))

# last session only
library(latex2exp)
fig <-
  sim_dat %>%
  dplyr::filter(session %in% c(-1)) %>%
  ggplot(aes(x=x, y=DA), fill="white") +
  geom_point(aes(color = id),  fill="white", size=1, alpha=0.9, shape = 21) +
  facet_wrap(id ~., nrow = num_subj, scales = "free_y") +
  theme_classic() +
  # mixed model individual curves
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            color = "black",
            alpha = 0.75) +
  geom_point(data = simPoint[simPoint$session == -1,], aes(x = trial, y = point),
             fill = "black", alpha = 0.75, pch = 21, size = 2) + # aes(fill = id), color = "black",
  geom_hline(data = simPoint[simPoint$session == -1,],
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
  ggtitle("Short Delay Anticipatory Licks") +
  xlab("Trial") +
  ylab("Total Anticipatory Licks")

# # save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "Delay_Length_Licks_Lastsession_OLS.pdf",
        plot = fig,
        width = 6,
        height = 12)


############################################################################################################
# plot average signal
############################################################################################################
# average signal
col_highlight <- "darkgreen"

# original data
mean_dat0 <- dat2 %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::mutate(time = round(as.numeric(time) / target_Hz - pre_min_tm, 2)) %>%
  group_by(time, id) %>%
  dplyr::summarise(photo_mean = mean(value, na.rm = TRUE)) %>% 
  # average across animals so se is uncertainty w.r.t animal averages
  group_by(time) %>%
  dplyr::summarise(photometry = mean(photo_mean, na.rm = TRUE),
                   y_low = mean(photo_mean, na.rm = TRUE) - sd(photo_mean) / sqrt(n()), # +/- sem
                   y_high = mean(photo_mean, na.rm = TRUE) + sd(photo_mean) / sqrt(n()),
                   sem = sd(photo_mean) / sqrt(n()) )

# trick variables to use plotting function to just plot the mean
mean_signal <- DA_delay
mean_signal$betaHat[2,] <- mean_dat0$photometry
diag(mean_signal$betaHat.var[,,2]) <- mean_dat0$sem
mean_signal$qn <- rep(0,3)

plot.f4 <- list()
prdtr <- 2
fit_dat <- mean_signal
plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               x_lab = "Time from Cue Onset (s)",
                               ylim = c(-0.5,7),
                               var_name = c("Intercept", "Average Signal")) + 
  ylab(latex2exp::TeX("Photometry Signal  $(\\Delta F / F)$" ))

# add interval names
plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                     x_interval_list = list(c(-1,0), c(0, 2)),
                                     x_interval_text = list(c(-1.25), c(1.75)), # 
                                     text_list = list("Baseline", "Cue Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     #ylim = c(-1, 6 * 1.01),
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G")
ggsave( "Photobleach Average Signal.pdf",
        plot = fig,
        width = 4,
        height = 4)

