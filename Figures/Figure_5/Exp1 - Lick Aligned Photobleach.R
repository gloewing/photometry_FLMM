# Lick aligned Experiment 1 
# Final Version
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

#########################################################################################
# create AUC variables for EDA

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
dat$trial <- (dat$trial - 1) / 100 # trial 50 is middle
dat$licks <- scale(dat$licks, center = TRUE, scale = FALSE)
dat$lick_time <- scale(dat$lick_time, center = TRUE, scale = FALSE)
dat$iri_orig <- dat$iri
dat$iri <- scale(dat$iri, center = TRUE, scale = FALSE)
dat$dat_AUC_pre <- dat_AUC_pre
dat$dat_AUC_post <- dat_AUC_post
dat$dat_AUC <- dat_AUC_post - dat_AUC_pre # described in "Experiment 1" of "Data Analysis" in Supplement
tm <- pre_min_tm * target_Hz # lick onset

if(short_iri_adjust){
  dat <- dat[!is.infinite(dat$bout_rate),] # 1 trial has infinite bout rate
  dat$bout_rate <- scale(dat$bout_rate, center = TRUE, scale = FALSE)
  dat <- dat[complete.cases(dat),] # remove one observation
}else{
  dat <- dat[complete.cases(subset(dat, select = -c(bout_events,bout_rate) )),] # remove one observation
  
}    
dim(dat)
iri_95 <- quantile(dat$iri, probs = 0.95)

########################################################
# flMM Modeling
########################################################
# ---------------------------------------------------
# Trial Number
# ---------------------------------------------------
# compare different random effect structures
# best model
DA_trial <- fui(photometry ~ session + trial +
                   (lick_time | id/session), 
                   data = dat,  
                   parallel = TRUE,
                   nknots_min = nknots_min,
                   subj_ID = "id")
#AIC/BIC: 13360.91 13421.66

# model 2
DA_trial2 <- fui(photometry ~ session + trial +
                  (1 | id/session), 
                data = dat,  
                parallel = TRUE,
                nknots_min = nknots_min,
                subj_ID = "id")
#AIC/BIC: 13471.74 13508.19

# model 3
DA_trial3 <- fui(photometry ~ session + trial +
                   (trial | id/session), 
                 data = dat,  
                 parallel = TRUE,
                 nknots_min = nknots_min,
                 subj_ID = "id")
#AIC/BIC: 13373.19 13433.93

# compare model fits
rbind(colMeans(DA_trial$aic[reward_period_idx,]),
      colMeans(DA_trial2$aic[reward_period_idx,]),
      colMeans(DA_trial3$aic[reward_period_idx,]))

# best model is model 1
fit_dat <- DA_trial

# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("Intercept", "Session Number Association", "Trial Number Association"))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                     x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Reward Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
plot.f4 <- plot_adjust(plot.f4)
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned Trial Effect.pdf",
        plot = fig,
        width = 4,
        height = 8)


############################################################################################################
# plot average signal
############################################################################################################
# average signal
col_highlight <- "darkgreen"

# average all sessions
sess_vec <- 1:8

# original data
mean_dat0 <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(session_orig %in% sess_vec) %>%
  dplyr::mutate(time = round(as.numeric(time) / target_Hz - pre_min_tm, 2)) %>%
  group_by(time, id) %>%
  dplyr::summarise(photo_mean = mean(value, na.rm = TRUE)) %>% 
  # average across animals so se is uncertainty w.r.t animal averages
  group_by(time) %>%
  dplyr::summarise(photometry = mean(photo_mean, na.rm = TRUE),
                   y_low = mean(photo_mean, na.rm = TRUE) - sd(photo_mean) / sqrt(n()), # +/- sem
                   y_high = mean(photo_mean, na.rm = TRUE) + sd(photo_mean) / sqrt(n()),
                   sem = sd(photo_mean) / sqrt(n()) )


p2 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +      
  geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, alpha=0.4)) + # fill = as.factor(session_orig), 
  geom_line(data = mean_dat0, colour = "black", size = 0.75) +  # colour = "black",
  xlab("Time from Lick (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Photometry Signal") +
  # scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.25),
             alpha = 0.7) +
  geom_vline(xintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.25),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.text.y = element_text(size=5, face="bold"), 
        axis.text.x = element_text(size=6, face="bold"), 
        axis.title = element_text(size=7, face="bold"),
        legend.position = "none",
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) 


############################################################################################################
# trick variables to use plotting function to just plot the mean
mean_signal <- DA_trial
mean_signal$betaHat[2,] <- mean_dat0$photometry
diag(mean_signal$betaHat.var[,,2]) <- mean_dat0$sem
mean_signal$qn <- rep(0,3)



plot.f4 <- list()
prdtr <- 2
fit_dat <- mean_signal
plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               x_lab = "Time from Lick (s)",
                               ylim = c(-0.5,6.25),
                               var_name = c("Intercept", "Average Signal", "Trial Number Association")) + 
                                ylab(latex2exp::TeX("Photometry Signal  $(\\Delta F / F)$" ))

# add interval names
plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                     x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Reward Period"),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     #ylim = c(-1, 6 * 1.01),
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )

# pre/post lick
plot.f4[[1]] <- interval_label(fig = plot.f4[[1]],
                               x_interval_list = list(c(-0.5, 0), c(0, 1)),
                               text_list = list("Pre-lick   ", "Post-lick"),
                               scl = 1.01, # percent for lines above original ylim values
                               x_scl = 0.005, # percent for gap in lines
                               txt_adjust = -0.0625, # percent text is above line
                               txt_size = 2.5,
                               col_seq = c("#0868ac", "#0868ac"),
                               y_val = 0.535, 
                               alpha = 0.75) +
  ggplot2::coord_cartesian(ylim = c(-0.2, layer_scales(plot.f4[[1]])$y$range$range[2]), #c(-0.2, 0.5642675 * 1.05),
                           xlim = c(-2, 2.5))

plot.f4[[2]] <- interval_label(fig = plot.f4[[2]],
                               x_interval_list = list(c(-0.5, 0), c(0, 1)),
                               text_list = list("Pre-lick   ", "Post-lick"),
                               scl = 1.01, # percent for lines above original ylim values
                               x_scl = 0.005, # percent for gap in lines
                               txt_adjust = -0.0625, # percent text is above line
                               txt_size = 2.5,
                               col_seq = c("#0868ac", "#0868ac"),
                               y_val = 0.6625, 
                               alpha = 0.75) +
  ggplot2::coord_cartesian(ylim = c(-1.6, layer_scales(plot.f4[[2]])$y$range$range[2]),
                           xlim = c(-2, 2.5))

# align plots
plot.f4 <- plot_adjust(plot.f4)
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned Average Signal.pdf",
        plot = fig,
        width = 4,
        height = 8)

############################################################################################################
# session stratified
sess_vec <- seq(1,8)

# original data
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(session_orig %in% sess_vec) %>%
  dplyr::mutate(time = round(as.numeric(time) / target_Hz - pre_min_tm, 2)) %>%
  dplyr::group_by(time, trial, session_orig) 


# average signal
mean_dat0 <- mean_data  %>%
  # average across trials within animal
  group_by(time, id, session_orig) %>%
  dplyr::summarise(photo_mean = mean(value, na.rm = TRUE)) %>% 
  # average across animals so se is uncertainty w.r.t animal averages
  group_by(time, session_orig) %>%
  dplyr::summarise(photometry = mean(photo_mean, na.rm = TRUE),
                   y_low = mean(photo_mean, na.rm = TRUE) - sd(photo_mean) / sqrt(n()), # +/- sem
                   y_high = mean(photo_mean, na.rm = TRUE) + sd(photo_mean) / sqrt(n()))





p2 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +      
  #geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, fill = as.factor(session_orig), alpha=0.4)) + # 
  geom_line(data = mean_dat0, aes(colour = as.factor(session_orig)),  size = 0.75) +  # colour = "black",
  xlab("Time from Lick (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Photometry Signal") +
  #scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.25),
             alpha = 0.7) +
  geom_vline(xintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.25),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.text.y = element_text(size=5, face="bold"), 
        axis.text.x = element_text(size=6, face="bold"), 
        axis.title = element_text(size=7, face="bold"),
        legend.position = "none",
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) 
# guides(color= guide_legend(title="Delay"), fill = "none") + 
# coord_cartesian(ylim = c(-0.25, 2.75),
#                 xlim = c(7.5,11.5)) + # , clip = "off"
# bars over top
# geom_segment(aes(x=0,xend=2,y=7,yend=7), color = "grey 20", size = 1) +
# # vertical dotted lines
# geom_segment(aes(x=example_time,xend=example_time,
#                  yend=mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Long"],#mean_data$value[mean_data$time == example_time],
#                  y=mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Short"]), 
#              color = col_highlight, linetype="dashed") +
# # # geom_segment(aes(x=2,xend=2,y=0,yend=6), color = "grey 20", linetype="dashed") +
# annotate(geom="text", x=8, y=1.1, label="Long", color = "#ca0020", fontface =2, size = 2.5) +
# annotate(geom="text", x=8, y=-0.2, label="Short", color = "#0868ac", fontface =2, size = 2.5) +
# geom_point(aes(x=example_time, y=peaks[1]), 
#            size = 1.5, colour="#ca0020") +
# geom_point(aes(x=example_time, y=peaks[2]), 
#            size = 1.5, colour="#0868ac") +
# scale_x_continuous(breaks = c(8, 9.4, 11),
#                    labels = c(8, 9.4, 11)) 
#scale_y_continuous(position = "right") + # put y-axis on right
# brace
# geom_brace(aes(c(9.55, 9.85), c(peaks[2], peaks[1]), label=" "),  # TeX("$\\beta_1(2.3)$ = Peak \n Difference")
#            labelsize=3.5, rotate = 90) + # inherit.data=FALSE, 
# annotate(geom="text", x=10.5, y=0.85, 
#          #label="Long - Short \nat 2.33s", 
#          label=latex2exp::TeX("$\\beta_1 (9.4) = 1.6", bold = TRUE),
#          color = col_highlight, fontface =2, size = 2.5) +
# annotate(geom="text", x=10.5, y=1.3, 
#          #label="Long - Short \nat 2.33s", 
#          label="Long - Short:",
#          color = col_highlight, fontface =1, size = 2.5) 


DA_trial <- fui(photometry ~ session + trial +
                  (lick_time | id/session), 
                data = dat[dat$session_orig %in% sess_vec,],  
                parallel = TRUE,
                nknots_min = nknots_min,
                subj_ID = "id")


# best model is model 1
fit_dat <- DA_trial

# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = tm,
                                 var_name = c("Intercept", "Session Number Association", "Trial Number Association"))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                       x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                       x_interval_text = NULL,
                                       text_list = list("Baseline", "Reward Period"),
                                       scl = 1.01, # percent for lines above original ylim values
                                       x_scl = 0.0025, # percent for gap in lines
                                       txt_adjust = 0.03, # percent text is above line
                                       txt_size = 3,
                                       col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
plot.f4 <- plot_adjust(plot.f4)
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned Trial Effect Photobleach.pdf",
        plot = fig,
        width = 4,
        height = 8)

# ---------------------------------------------------
# IRI
# ---------------------------------------------------
# compare different random effect structures
# best model
DA_iri <- fui(photometry ~ iri + 
                  (lick_time | id/session), 
                data = dat,  
                parallel = TRUE,
                nknots_min = nknots_min,
                subj_ID = "id")
#AIC/BIC: 13362.52 13417.19

# model 2
DA_iri2 <- fui(photometry ~ iri + 
                   (1 | id/session), 
                 data = dat,  
                 parallel = TRUE,
                 nknots_min = nknots_min,
                 subj_ID = "id")
#AIC/BIC: 13472.49 13502.86

# model 3
DA_iri3 <- fui(photometry ~ iri + 
                   (iri | id/session), 
                 data = dat,  
                 parallel = TRUE,
                 nknots_min = nknots_min,
                 subj_ID = "id")
#AIC/BIC: 13367.15 13421.83

# compare model fits
rbind(colMeans(DA_iri$aic[reward_period_idx,]),
      colMeans(DA_iri2$aic[reward_period_idx,]),
      colMeans(DA_iri3$aic[reward_period_idx,]))

# best model is model 1
fit_dat <- DA_iri

# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = tm,
                                 var_name = c("Intercept", "Inter-reward Interval Association"))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                       x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                       x_interval_text = NULL,
                                       text_list = list("Baseline", "Reward Period"),
                                       scl = 1.01, # percent for lines above original ylim values
                                       x_scl = 0.0025, # percent for gap in lines
                                       txt_adjust = 0.03, # percent text is above line
                                       txt_size = 3,
                                       col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252") )
}

# align plots
figs <- list()
figs[[1]] <- plot.f4[[1]]
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned IRI Effect.pdf",
        plot = fig,
        width = 4,
        height = 4)


# example animal

#######################################################
# LME
#######################################################
# lme model pooled
mod <- lme4::lmer(dat_AUC ~ iri + 
                    ( scale(lick_time) | id/session_orig), #
                  data = dat)
m <- summary(mod)
sim_dat <- data.frame(id = dat$id, session = dat$session_orig, trial = dat$trial_orig, 
                      DA = dat$dat_AUC, iri = dat$iri, rn = dat$reward_num,
                      x = dat$trial_orig )

sim_dat$fitted <- fitted(mod) # fitted values with BLUPs
sim_dat$fixed <- predict(mod,re.form=NA,newdata=dat) # fixed effect only fitted values
sim_dat <- sim_dat[sim_dat$session %in% 1:8,] # easier to visualize

# fit this to get p-values
mod <- nlme::lme(dat_AUC ~ reward_num,
                 random = ~scale(lick_time)| id/session_orig, 
                 data = dat) # [dat$session_orig <= 2,]
m <- summary(mod)
beta1 <- round(mod$coefficients$fixed[2], 4)
beta1_p <- round(m$tTable[2,5], 4)
rm(mod)

##################
# example animal
##################
iri_95_orig <- quantile(dat$iri_orig, probs = 0.95)

dd <- dat %>%
  dplyr::filter(id == unique(dat$id)[1], # first animal arbitrarily
                iri_orig <= iri_95_orig) 

cor.test(dd$iri_orig, dd$dat_AUC)
# p-value = 2.511e-10
# correlation coefficient 0.2533079

# Reward Number LMM pooled across sessions and animals
corPlot3 <- 
  dd %>%  # remove IRI outliers for visualization
  ggplot(aes(x=iri_orig, y=dat_AUC), fill="white") +
  geom_point(fill="white", size=1, alpha=0.9, shape = 21) + # aes(color = id),  
  theme_classic() + 
  # mixed model individual curves
  # geom_line(aes(color = id, x = rn, y = fixed), color = "black",#aes(group = id, color = id, x = x, y = fixed),
  #           size = 1) +
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            color = "black",
            alpha = 0.75) +
  theme(axis.text = element_text(size=12, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        plot.title = element_text(hjust = 0.5, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        legend.title = element_text(face = "bold", size = 12)) +
  guides(fill = "none", alpha = "none", 
         color = "none") + 
  coord_cartesian(ylim = c(-5, 7)) +
  ggtitle("IRI-DA Correlation in Example Animal") +
  xlab("Inter-reward Interval") +
  ylab("Reward Period DA") +
  annotate(geom = "text", x = 15, y = 5, 
           label = TeX("$p < 2.5e-10$, $\\hat{\\rho} = 0.253$"),
           color = "black", size = 3.75)

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "example_IRI_correlation.pdf",
        plot = corPlot3,
        width = 4,
        height = 4)
# ---------------------------------------------------
# Lick Rate
# ---------------------------------------------------
# compare different random effect structures
# best model
DA_br <- fui(photometry ~ bout_rate + 
                  (lick_time | id/session), 
                data = dat,  
                parallel = TRUE,
                nknots_min = nknots_min,
                subj_ID = "id")
# AIC/BIC: 13409.71 13464.38

# model 2
DA_br2 <- fui(photometry ~ bout_rate + 
                   (1 | id/session), 
                 data = dat,  
                 parallel = TRUE,
                 nknots_min = nknots_min,
                 subj_ID = "id")

# AIC/BIC: 13522.94 13553.31

# model 3
DA_br3 <- fui(photometry ~ bout_rate + 
                   (trial | id/session), 
                 data = dat,  
                 parallel = TRUE,
                 nknots_min = nknots_min,
                 subj_ID = "id")

#AIC/BIC: 13380.74 13435.42

# compare model fits
rbind(colMeans(DA_br$aic[reward_period_idx,]),
      colMeans(DA_br2$aic[reward_period_idx,]),
      colMeans(DA_br3$aic[reward_period_idx,]))

# model 3 is best

ylimList <- list(NULL, c(-0.2, 0.2))

fit_dat <- DA_br3
# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               fig = fit_dat,
                               Hz = target_Hz, 
                               align = tm,
                               ylim = ylimList[[prdtr]],
                               var_name = c("Mean Signal on Mean Lick Rate Trial", "Lick Rate Association")) #, "Session", "Bout Rate x Session"

  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                     x_interval_text = NULL,
                                     text_list = list("Baseline", "Reward Period"),
                                     scl = 1.02, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("#ca0020", "#0868ac", "#E69F00", "#525252"),
                                     ylim = ylimList[[prdtr]])
  
}

# align y-axes
figs[[2]] <- plot.f4[[2]]
fig <- do.call("grid.arrange", c(plot.f4, nrow = 2))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned Lick Rate Effect.pdf.pdf",
        plot = fig,
        width = 4,
        height = 8)





# remove intercepts and join across trial ranges
figs_join <- plot_adjust(figs)
fig <- do.call("grid.arrange", c(figs_join, nrow = 2))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Exp1/Final/Lick Aligned")
ggsave( "Lick-Aligned Lick Rate_IRI Effect.pdf",
        plot = fig,
        width = 4,
        height = 8)
