library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(refund)
library(dplyr)

# fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fReg_new.R')
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/interval_label.R') # plot
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_adjust.R') # plot
path <- "/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/data/science_paper/fig_6/data/"

# read data and extract variable names
setwd("~/Desktop/NIMH Research/Causal/DA_adapts_learn_rate")
matdat <- R.matlab::readMat("seshMerge.mat")
matdat <- matdat$seshMerge


# indices in dataset of variables
photo_idx <- 7
opt_idx <- 4
session_idx <- 2
trial_idx <- 1
sess_idx <- c(1:6)
iti_idx <- 15
Hz <- 100
target_Hz <- 20

# indices in photometry data of time points (e.g., cue onset/offset)
# cued trials (0 and 2) -- not sure about cued trials==1
cue_onset <- 151
cue_offset <- 200
post_cue_end <- 250 # analyze window from cue_onset:post_cue_end according to Luke
reward_onset <- 301
reward_off <- 500

# animal group ids
control <-  c(1, 2, 4, 6, 9, 11, 15, 16, 19)
slPlus <- c(7, 8, 12, 13, 17) 
slMinus <- c(3, 5, 10, 14, 18, 20) 
sPlusLPlus <- c(21, 22, 23, 24, 25)
n <- 24 # sample size


id_mat <- data.frame( ids = 1:n, group = rep(NA, n) )
id_mat[control,2] <- "control"
id_mat[slPlus,2] <- "slPlus"
id_mat[slMinus,2] <- "slMinus"
id_mat[sPlusLPlus,2] <- "sPlusLPlus"

# session data
nms <- rownames(matdat[,1,])
cue_nms <- nms[sess_idx] # names with same dimension about trial type
dat <- do.call(cbind, lapply(sess_idx, function(x) do.call(c, matdat[x,,])  ) ) # seq_along(cue_nms)
colnames(dat) <- cue_nms
ids <- do.call(c, lapply(seq_along(1:n), function(x) rep(x, length( matdat[1,,][[x]])  ) ))
trial_num <- do.call(c, sapply(seq_along(1:n), function(x)  1:length( ids[ids == x] ) )   ) # trial number (ignores session structure) # assumes in order within animal ID
ids <- data.frame(ids = ids, trial_num = trial_num)

# treatment groups
ids <- left_join(ids, id_mat, by = "ids")

# ITIs (not provided in all sessions)
ITIs <- vector(length = nrow(ids) )
iti_vec <- do.call(c, matdat[iti_idx,,])  # ITIs 
ITIs[1:length(iti_vec)] <- iti_vec # assumes (correctly) ITIs happen for first part of vector and then ends

# photometry
photo <- do.call(rbind, matdat[photo_idx,,] )
colnames(photo) <- paste0("photometry.", 1:ncol(photo) )

# summary of photometry for outcomes
### In cued trials (trialID==0 or 2), cues start at data point 151 in the trials and end at data point 200 (i.e., 0.5 sec cue that ends 1 s before reward delivery)
## "Analysis windows were chosen to capture the extent of mean phasic activations following each kind of stimulus. For NAc–DA and VTA–DA, reward responses were quantified from 0 to 2 s after reward delivery and 
## cue responses were quantified from 0 to 1 s after cue delivery. DS–DA exhibited much faster kinetics, and thus reward and cue responses were quantified from 0 to 0.75 s after delivery."
# 2) Quantify DA (subtract out mean baseline): 
#       -mean baseline DA[1:150]
#       -mean cue DA[151:250] # includes
#       -mean reward DA[301:500]

base_DA <- rowMeans( photo[ , 1:(cue_onset-1)] ) # baseline DA before cue onset
cue_DA <- rowMeans( photo[ , cue_onset:post_cue_end] ) # cue DA 
reward_DA <- rowMeans( photo[ , reward_onset:reward_off] ) # reward DA 
cue_DA_max <- apply( photo[ , cue_onset:post_cue_end], 1, max ) # cue DA 
reward_DA_max <- apply( photo[ , reward_onset:reward_off], 1, max ) # reward DA 

# DA during periods relative to pre-cue baseline
DA_mat <- data.frame( cue_DA =    cue_DA - base_DA,
                      reward_DA = reward_DA - base_DA,
                      cue_DA_max = cue_DA_max - base_DA, # t**akes max during period and subtracts off MEAN during baseline**
                      reward_DA_max = reward_DA_max - base_DA,
                      cue_DA_noBase =    cue_DA,
                      reward_DA_noBase = reward_DA,
                      cue_DA_max_noBase = cue_DA_max, # t**akes max during period and subtracts off MEAN during baseline**
                      reward_DA_max_noBase = reward_DA_max)

# join data together
dat <- as.data.frame( cbind(ids, ITIs, dat, DA_mat) )
dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(iti_vec, ids, ITIs, DA_mat, photo, trial_num)


# downsample photometry
pre_samps <- 150
out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
cov_index <- seq(1, ncol(dat))[-out_index]
L <- 701
by_num <- round( Hz / target_Hz )
photo_idx_pre <- sort(-seq(-pre_samps, -1, by = by_num)) # seq(1, pre_samps, by = by_num)
photo_idx_post <- seq(pre_samps + by_num, L, by = by_num)
photo_idx <- unique( c(photo_idx_pre, photo_idx_post) )
dat <- dat[ , c(cov_index, out_index[photo_idx]) ]
L <- length(photo_idx)
nknots_min <- round(L/2)
pre_min_tm <- 1.5
colnames(dat[,-cov_index]) <- paste0("photometry.", 1:L)

########################################
# correlation b/w prep behavior and DA
########################################
# Unexpectedly, initial NAc–DA reward signals were neg-
# atively correlated with the amount of preparatory behaviour at the end of training (Fig. 2c and Extended
# Data Fig. 3f; NAc–DA rewardtrials 1–100 versus preparatory indextrials 700–800, r=0.85, P=0.004)

prep700_800 <- dat %>% 
  as_tibble() %>%
  dplyr::filter(trial_num >= 700,
                trial_num <= 800,
                stimState == 0, # no stim
                latency <= 1000, # remove crazy outliers (< 1%)
                trialID == 0 # CUED reward
                ) %>%  
  dplyr::group_by(ids) %>%
  dplyr::summarise(prepAvg = mean(lickState, na.rm = TRUE),
                   latencyAvg = mean(latency, na.rm = TRUE)) # average prep-behavior

dat <- left_join(dat, prep700_800, by = "ids")

dat2 <- dat %>% dplyr::filter(group == "control",
                              trialID == 0, # CUED reward
                              stimState == 0,
                              trial_num <= 100)
dat2$seshID <- as.factor(dat2$seshID)
dat2$ids <- as.factor(dat2$ids)
dat2$trial_num <- as.numeric(dat2$trial_num)
# functional LME

# reproduce figure 2e
dat2 %>% dplyr::filter(group == "control",
                       trialID == 0,
                       stimState == 0,
                       trial_num <= 100) %>%
  group_by(ids) %>%
  summarise(mean_reward_DA = mean(reward_DA, na.rm = TRUE), # reward_DA_max also works
            latencyAvg = latencyAvg  # latency is in miliseconds
            ) %>% # this is already averaged so no need to take any functions of it
  unique() %>%
  ggplot(aes(x = mean_reward_DA, y = latencyAvg, color = ids, group = ids)) +
            geom_point()

##########################################################
# average fosr -- average latency
dat3 <- dat2 %>%
  dplyr::filter(latency <= 1000,
                group == "control",
                trialID == 0,
                stimState == 0,
                trial_num <= 100) %>%
  dplyr::group_by(ids) %>%
  dplyr::select(-group, -seshID) %>%
  dplyr::summarise_all(mean) # throws errors for IDs and non-numeric that can be ignored safely

L <- sum( grepl( "photometry" , names( dat3 )) )
photo_dat <- as.matrix(dat3[, grepl( "photometry" , names( dat3 ) )])

dat3_x <- as.matrix(scale(dat3$latencyAvg / 1000, center = TRUE, scale = FALSE) ) # for interpretability -- dividing by 1000 converts units into seconds
tm <- 1.5 * target_Hz # lick onset

# Penalized OLS with smoothing parameter chosen by grid search
olsmod = refund::fosr(Y = photo_dat, X = cbind(1, dat3_x), method="OLS", nbasis = round(L/4)) # 
fit_dat <- olsmod

# plot
plot.f4 <- list()
for(prdtr in 1:ncol(olsmod$se.func)){
  plot.f4[[prdtr]] <- plot.freg(olsmod,
                                r = prdtr, 
                                Hz = target_Hz, 
                                align = tm,
                                var_name = c("Mean Signal: Avg. Latency Subject", 
                                             "Final Latency"),
                                y_scal_orig = 0.01) 
  
  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1.5,0), c(0, 1), c(1.5,3.5)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Cue", "Reward"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
  
  
  # centered to 1.5 sec so above all have 1.5 subtracted off
  #       -mean baseline DA[1:150]
  #       -mean cue DA[151:250] # includes
  #       -mean reward DA[301:500] 
}

fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/da_adapts_learn/Final")
ggsave( "da_causal_100trialsDA_predictsBehavior_fosr.pdf",
        plot = fig,
        width = 4,
        height = 8)

##################
# prepatory lick probability
dat3_x <- as.matrix(scale(dat3$prepAvg, center = TRUE, scale = FALSE) ) # for interpretability 

# Penalized OLS with smoothing parameter chosen by grid search
olsmod = refund::fosr(Y = photo_dat, X = cbind(1, dat3_x), method="OLS", nbasis = round(L/4)) # 
fit_dat <- olsmod

tm <- 1.5 * target_Hz # lick onset

# plot
plot.f4 <- list()
for(prdtr in 1:ncol(olsmod$se.func)){
  plot.f4[[prdtr]] <- plot.freg(olsmod,
                               r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Mean Signal: Avg. Lick Prob. Subject", 
                                            "Lick Probability"),
                               y_scal_orig = 0.01) 
  
  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1.5,0), c(0, 1), c(1.5,3.5)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Cue", "Reward"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )

}



fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))


# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/da_adapts_learn/Final")
ggsave( "da_causal_100trialsDA_predictsPrepLick_fosr.pdf",
        plot = fig,
        width = 4,
        height = 8)

##########################################################
# Figure 2e 
##########################################################
# Notably, individual differences in initial NAc–DA reward signals were not correlated with 
# the learning of NAc–DA cue signals (rewardtrials 1–100 versus cuetrials 700–800, P = 0.5; Fig. 2e). 
dat2 <- dat %>% dplyr::filter(group == "control",
                              trialID == 0, # CUED reward
                              stimState == 0,
                              latency <= 1000)
dat2$seshID <- as.factor(dat2$seshID)
dat2$ids <- as.factor(dat2$ids)
dat2$trial_num <- as.numeric(dat2$trial_num)

# average signal over first 100 trials
dat_100 <- dat2 %>%
  dplyr::filter(trial_num <= 100) %>%
  dplyr::group_by(ids) %>%
  dplyr::select(-group, -seshID) %>%
  dplyr::summarise_all(mean) # throws errors for IDs and non-numeric that can be ignored safely

dat100_reward <- as.numeric(dat_100$reward_DA) # average reward period DA (mean across trials)
ids_100 <- as.factor(dat_100$ids)

# average signal over last 100 trials
dat_700 <- dat2 %>%
  dplyr::filter(trial_num <= 800,
                trial_num >= 700) %>%
  dplyr::group_by(ids) %>%
  dplyr::select(-group, -seshID) %>%
  dplyr::summarise_all(mean) # throws errors for IDs and non-numeric that can be ignored safely

dat700_cue <- as.numeric(dat_700$cue_DA) # average cue across trials

# combine
photo_dat100 <- as.matrix(dat_100[, grepl( "photometry" , names( dat_100 ) )])
photo_dat700 <- as.matrix(dat_700[, grepl( "photometry" , names( dat_700 ) )])


# fosr: predict entire 700-800 trials (function) from reward period 100 OLS with smoothing parameter chosen by grid search
dat2_fd <- tibble(photo700 = photo_dat700,
                  photo100 = dat100_reward,
                  ids = as.factor(ids_100)) %>%
  as.data.frame()

# fosr model
olsmod <- refund::fosr(photo700 ~ photo100, data = dat2_fd, method="OLS", nbasis = round(L/4)) # 
fit_dat <- olsmod

tm <- 1.5 * target_Hz # lick onset

# plot
plot.f4 <- list()
for(prdtr in 1:ncol(olsmod$se.func)){
  plot.f4[[prdtr]] <- plot.freg(olsmod,
                                r = prdtr, 
                                Hz = target_Hz, 
                                align = tm,
                                var_name = c("Mean Signal on Avg. Reward DA Subject", 
                                             "Reward DA"),
                                y_scal_orig = 0.01) 
  
  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1.5,0), c(0, 1), c(1.5,3.5)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Cue", "Reward"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
  
}


# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/da_adapts_learn/Final")
ggsave( "da_causal_100trialsRewardDA_Avg_predictsCueDA700_fosr.pdf",
        plot = fig,
        width = 4,
        height = 8)

# scalar-on-function regression: predict entire 700-800 trials (function) from reward period 100 OLS with smoothing parameter chosen by grid search
dat2_fd <- tibble(photo700 = dat700_cue,
                  photo100 = photo_dat100,
                  ids = as.factor(ids_100)) %>%
  as.data.frame()

# predict functional trials 1-100 based on scalar of trials 700-800
olsmod = refund::fosr(photo100 ~ photo700, data = dat2_fd, method="OLS", nbasis = round(L/4)) # 
fit_dat <- olsmod

# plot
plot.f4 <- list()
for(prdtr in 1:ncol(olsmod$se.func)){
  plot.f4[[prdtr]] <- plot.freg(olsmod,
                                r = prdtr, 
                                Hz = target_Hz, 
                                align = tm,
                                var_name = c("Mean Signal on Avg. Cue DA Subject", 
                                             "Cue DA"),
                                y_scal_orig = 0.01) 
  
  
  # # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1.5,0), c(0, 1), c(1.5,3.5)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Cue", "Reward"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
  
}

fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))


# save
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/da_adapts_learn/Final")
ggsave( "da_causal_predicts_100trialsRewardDA_CueDA700Avg_fosr.pdf",
        plot = fig,
        width = 4,
        height = 8)
