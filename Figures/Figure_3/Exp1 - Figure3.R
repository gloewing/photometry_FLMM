# show temporal correlation figures
# Lick aligned Experiment 1 
# Final Version
library(data.table)
library(dplyr)
library(ggplot2)

# fGLMM code
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_smooth.R") # function code
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_multi.R") # function code
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/lfosr3s_wrap.R") # function code

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
# consider centering iri, trial, etc, so intercept is interpretable and more likely to be significant

# photometry sample length
out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
L <- length(out_index)
nknots_min <- round(L/2)
cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates

#########################################################################################
# creat AUC variables for EDA

photo_dat <- dat[, -cov_idx] # photometry data
pre_lick_reward_period_samps <- round(target_Hz * pre_reward_period_length) # samples of pre_cue period
reward_period_max_idx <- round( (pre_reward_period_length + reward_period_length) * target_Hz ) # max_idx of post_cue period
dat_AUC_pre <- rowMeans( photo_dat[, 1:pre_lick_reward_period_samps] ) # AUC before reward period
dat_AUC_post <- rowMeans(photo_dat[, seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)]) # AUC after reward period

reward_period_idx=seq(pre_lick_reward_period_samps + 1, reward_period_max_idx)

# scale covariates to make contrasts more interpretable
dat$reward_num <- dat$trial + (dat$session - 1) * num_trials # number of rewards assuming 100 rewards/session
dat$licks <- scale(dat$licks, center = TRUE, scale = FALSE)
dat$lick_time <- scale(dat$lick_time, center = TRUE, scale = FALSE)
dat$dat_AUC_pre <- dat_AUC_pre
dat$dat_AUC_post <- dat_AUC_post
dat$dat_AUC <- dat_AUC_post - dat_AUC_pre # described in "Experiment 1" of "Data Analysis" in Supplement
dat <- dat[!is.infinite(dat$bout_rate),] # 1 trial has infinite bout rate
dat$bout_rate <- scale(dat$bout_rate, center = TRUE, scale = FALSE)
dim(dat)
dat <- dat[complete.cases(dat),] # remove one observation

dat %>% as_tibble() %>% group_by(session) %>% summarise(m = n())
dat$id <- as.factor(as.numeric(as.factor(dat$id)))
# -----------------------------------------------------------------
# plots

# plot session boxplots (with trial level jitter points) of reward-DA
f1 <- dat %>% 
  as_tibble() %>% 
  group_by(session, id) %>% 
  dplyr::filter(session <= 8) %>%
  ggplot(aes(y = dat_AUC, x = as.factor(session), group = as.factor(session))) +
  geom_boxplot(lwd = 1,
               fatten = 0.5) + 
  geom_jitter(aes(colour = id),
              size = 2, shape = 21, alpha = 0.65, width = 0.3) +
  ylab("Reward Period DA" ) + 
  xlab( "Session" ) + 
  scale_fill_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black") ) +
  scale_color_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black")) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(-2, 8) ) + 
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(1, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5)) ) + 
  guides(colour= guide_legend(title="Animal"))  
# 
#--------------------------------------------
# plot trial lines of reward-DA
f3 <- dat %>% 
  as_tibble() %>% 
  group_by(session, id) %>% 
  dplyr::filter(session == 3,
                trial <= 50) %>%
  ggplot(aes(y = dat_AUC, x = as.factor(trial), group = as.factor(trial))) +
  geom_point(aes(colour = id),
             size = 2, shape = 21) +
  geom_smooth(aes(group = id, colour = id),
              method = "loess",
              span = 0.3,
              alpha = 0.65,
              se = FALSE) +
  ylab("Reward Period DA" ) + 
  xlab( "Trial" ) + 
  scale_fill_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black") ) +
  scale_color_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black")) +
  theme_classic(base_size = 12) +
  coord_cartesian(ylim = c(-2, 8) ) + 
  scale_x_discrete(labels=c(1, seq(10, 100, by = 10)),
                   breaks=c(1, seq(10, 100, by = 10))) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(1, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5))
         ) + 
  guides(colour= guide_legend(title="Animal"))  

#-----------------------------------------------------------
iri_95 <- quantile(dat$iri, probs = 0.95) # remove outliers

f2 <- 
  dat %>% as_tibble() %>% 
  dplyr::filter(session <= 6,#,
                session >= 3,
                iri <= iri_95) %>% # some high leverage points
  ggplot(aes(y = dat_AUC, x = iri)) +
  geom_point(aes(colour = id),
             size = 2, shape = 21, alpha = 0.65) + 
  facet_wrap( ~ session, nrow = 1) +
  geom_hline(yintercept=0.95,
             linetype="dashed",
             color = "black",
             size = rel(0.5),
             alpha = 0.7) + #
  ylab("Reward Period DA" ) + 
  xlab( "Inter-Reward Interval (IRI)" ) + 
  geom_smooth(aes(group = id, colour = id),
              method = "lm",
              alpha = 0.65,
              se = FALSE
  ) +
  scale_x_continuous( #
                     breaks=c(0,10,20,30)) +
  scale_fill_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black") ) +
  scale_color_manual(values = c("#ca0020", "lightgrey", "#0868ac", "darkgray", "#E69F00", "#525252", "orange", "black")) +
  theme_classic(base_size = 12) +
  theme( plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
         axis.text=element_text(face="bold",color="black", size=rel(2)),
         axis.title = element_text(face="bold", color="black", size=rel(2)),
         legend.key.size = unit(1, "line"), # added in to increase size
         legend.text = element_text(face="bold", color="black", size = rel(2)), # 
         legend.title = element_text(face="bold", color="black", size = rel(2)),
         strip.text.x = element_text(face="bold", color="black", size = rel(2.5))
  ) + guides(colour= guide_legend(title="Animal"))  

#----------------------------------------------------
# sample size figure

max_trials <- max(dat$trial)
# data on sessions
plot <- dat %>% 
  as_tibble() %>% 
  mutate(IDs = as.numeric(as.factor(id)),
         Session = session) %>% 
  group_by(IDs, Session) %>% 
  summarise(Trials = n()  ) %>%
  # heatmap plot
  ggplot(aes(y = IDs, x = Session, fill = Trials)) +
  geom_tile() + 
  theme_classic() +
  scale_fill_gradient(low="white", high="red") +
  theme(plot.title = element_text(hjust = 0.5, color="black", size=rel(2.5), face="bold"),
        axis.text=element_text(color="black", size=rel(2), face="bold"),
        axis.title = element_text(face="bold", color="black", size=rel(2)),
        legend.key.size = unit(2, "line"), # added in to increase size
        legend.text = element_text(face="bold", color="black", size = rel(2)), # 
        legend.title = element_text(face="bold", color="black", size = rel(2)),
        strip.text.x = element_text(face="bold", color="black", size = rel(2.5))
  ) + guides(fill= guide_legend(title="Trials")) +
  ylab('Animal ID') +
  xlab("Session Number") +
  coord_cartesian(xlim = c(1, 11) )  + 
  scale_x_continuous(breaks= scales::pretty_breaks(n=11)) + 
  scale_y_continuous(breaks= scales::pretty_breaks(n=n))
#----------------------------------------------------


# combine plots
library(patchwork)
set.seed(25)
plot_combined <- f1 + f3 + f2 + plot + 
  plot_layout(ncol = 2,
              nrow = 2,
              byrow = NULL,
              widths = NULL,
              heights = 10,
              guides = NULL)

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Figure 3")
ggsave( "Exp1_grid.pdf",
        plot = plot_combined,
        width = 20,
        height = 12)
