# FLMM Explanation Inset for explanation
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
iri_95 <- quantile(dat$iri_orig, probs = 0.95)
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

# ---------------------------------------------------
# lick_time
# ---------------------------------------------------
# simple example model
DA_iri <- fui(photometry ~ lick_time + 
                (1 | id/session), 
              data = dat,  
              parallel = TRUE,
              nknots_min = nknots_min,
              subj_ID = "id")
# ---------------------------------------------------
#----------------------
# plot random effects
#----------------------
# simple example model
DA_re <- fui(photometry ~ lick_time + 
                (lick_time | id), 
              data = dat,  
              parallel = TRUE,
              nknots_min = nknots_min,
              REs = TRUE,
              subj_ID = "id")
# ---------------------------------------------------

# custom plot.FUI
plot.FUI <- function(r, 
                     fig = NULL, 
                     align, 
                     Hz, 
                     var_name = NULL, 
                     title = NULL,
                     y_val_lim = 1.1,
                     ylim = NULL,
                     xlim = NULL,
                     col = "black",
                     y_scal_orig = 0.05){
  library(gridExtra)
  name = NULL
  if(!is.null(fig))     fit_dat <- fig
  
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaHat[r,],
                             lower = fit_dat$betaHat[r,] - 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaHat[r,] + 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaHat[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaHat[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          ) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), 
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"=col)) + # "blue3"
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  
  p.CI <- p.CI + 
    labs(x = "Trial Time-point (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
         title = var_name[r]) +
    theme(legend.position = "none",
          axis.title.y = element_text(size=16, face="bold")
          ) #+
  
  # make x and y intercepts
  if(!is.null(ylim)){
    p.CI <- p.CI + coord_cartesian(ylim = ylim)
  }else{  
    ylim <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint)) #layer_scales(p.CI)$y$range$range
    y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint)) #layer_scales(p.CI)$y$range$range
    ylim[1] <- ylim[1] - y_adjust # just scale bottom because top is scaled below
  }     
  
  if(is.null(xlim)){
    xlim = layer_scales(p.CI)$x$range$range
  }
  
  x_range <- diff(xlim) * 0.1
  y_range <- diff(ylim) * 0.1
  y_range_up <- diff(ylim) * 0.02
  
  # extend upper limit
  y_val_lim <- c(1, y_val_lim)
  y_top <- (0.975) * diff(ylim*y_val_lim) + ylim[1]*y_val_lim[1]
  
  p.CI <- p.CI + 
    coord_cartesian(ylim = ylim*y_val_lim,
                    xlim = xlim) + 
    geom_segment(aes(x=xlim[1], xend=xlim[2],
                     y=0,yend=0),
                 color = "black", lwd=0.5, alpha = 0.5, linetype = "dashed") #+ # x - intercept
  geom_segment(aes(y=ylim[1], yend=y_top,#ylim[2]+y_range_up,
                   x=0,xend=0), # don't extend up
               color = "white", lwd=0.5, alpha = 0, linetype = "dashed")
  
  return(p.CI)
}

# best model is model 1
fit_dat <- DA_iri
b_iri <- DA_iri$betaHat[2,]
col_highlight <- "red"

# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = 1,
                                 col = col_highlight,
                                 var_name = c("Intercept", " "))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                       x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                       x_interval_text = NULL,
                                       text_list = list("Baseline", "Reward Period "),
                                       scl = 1.01, # percent for lines above original ylim values
                                       x_scl = 0.0025, # percent for gap in lines
                                       txt_adjust = 0.03, # percent text is above line
                                       txt_size = 3,
                                       col_seq = c("white", "white") #c("#ca0020", "#0868ac", "#E69F00", "#525252")
  )
}

# align plots
figs <- list()
figs[[1]] <- plot.f4[[1]]
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

###########################
# make inset figure for p1
###########################
n_samp <- trial_size <- 10 # 10 trials per animal
nn <- 4
tm <- 1
time_point <- 1.72 # time after lick
time_idx <- tm + round(time_point * target_Hz) 
beta1 <- b_val <- fit_dat$betaHat[2,time_idx+1]
fit_dat$betaHat[2,time_idx-1] # b-value shown
beta0 <- fit_dat$betaHat[1,time_idx]
dd <- dat[1:(nn*trial_size),]
dd$Y <- dd[paste0("photometry.",time_idx)]
dd$Y <- as.numeric(dd$Y[,1])
dd$id <- rep(1:nn, each = trial_size)
dot_color <- "darkgrey"

# simulate data according to fit for easy visualization
set.seed(1)
dd$sim <- dd$lick_time * (beta1-0.3) + rnorm(nrow(dd), 0, 0.15) 

###########################################
# example animal
beta1_2 <- b_val_2 <- DA_re$betaHat[2,time_idx+1] + DA_re$randEff[animal_highlight, 2,time_idx+1]
DA_re$betaHat[2,time_idx-1] # b-value shown
beta0 <- DA_re$betaHat[1,time_idx]

# simulate data according to fit for easy visualization
set.seed(1)
id_num <- unique(dd$id)[animal_highlight]
dd$sim[dd$id == 1] <- dd$lick_time[dd$id == 1] * beta1_2 + rnorm(nrow(dd[dd$id == 1,]), 0, 0.15) 
#############################################

mod_coef <- coef(lm(sim ~ lick_time, data = dd))

corPlot3 <- 
  dd %>% 
  as_tibble() %>%
  dplyr::mutate(Y = sim, 
                IRI = lick_time) %>%
  dplyr::select(Y, id, IRI) %>%  # 
  # remove IRI outliers for visualization
  ggplot(aes(x=IRI, y=Y), fill="white") +
  geom_point(fill="white", size=1, alpha=0.75, 
             color = dot_color) + #
  theme_classic() + 
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            color = col_highlight) +
  theme(axis.text = element_blank(),#element_text(size=4, face="bold"), 
        axis.title = element_text(size=8.5, face="bold"),
        plot.title = element_blank(),#element_text(hjust = 0.5, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        legend.title = element_text(face = "bold", size = 12)) +
  guides(fill = "none", alpha = "none", 
         color = "none") + 
  xlab("Latency") +
  ylab(latex2exp::TeX("\\textbf{Y(1.7):  $\\Delta F/F$ at 1.7sec}")) + #TeX("\\textbf{Y(0.1):}  \\Delta \\textbf{F / F Value at 0.1s}")) + #("Y(0.1):  Value at 0.1s") +
  annotate(geom = "text", x = 0.8, y = 0.85, 
           label = latex2exp::TeX("\\textbf{Slope:}  $\\beta_1 (1.7) = 0.82", bold = TRUE),
           angle = 39.5, # atan(beta1) * ( 180.0 / pi )
           color = col_highlight, size = 3.5) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5 ,5.5), "pt")) 
  
  

set.seed(25) # for jitter (position of highlighted dot)
p1_combined <- plot.f4[[1]] + 
  geom_point(aes(x=time_point, y=b_val), size = 3, colour=col_highlight) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5 ,5.5), "pt")) +
  geom_segment(aes(x =time_point, y = b_val, xend = time_point, yend = -0.71), #yend = -0.825),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75)  + 
  geom_segment(aes(x =-0.2, y = b_val, xend = time_point, yend = b_val),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75) +
  coord_cartesian(xlim = c(0, 4.5), 
                  ylim = c(-0.6, 1.85),  
                  clip = 'off') +
  scale_x_continuous(breaks = c(0, 4),
                     labels = c(0, 4)) +
  scale_y_continuous(breaks = c(0,  1.5),
                     labels = c(0, 1.5)) +
  annotate(geom = "text", x = -0.5, y = beta1 * 1.01, 
           label =  "0.82",
           color = col_highlight, size = 4) +
  annotate(geom = "text", x = time_point, y = -0.8, 
           label =  "1.7",
           color = col_highlight, size = 4) +
  patchwork::inset_element(corPlot3, 0.55, 0.51, 1, 1, align_to = 'full')

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration/Final")
ggsave( "FLMM Explanation Inset.pdf",
        plot = p1_combined,
        width = 4,
        height = 4)

l_combined <- list()
l_combined[[1]] <- p1_combined
######################################################
# random slope
######################################################
source('~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/plot_fGLMM_RE.R') # plot random effects
# ---------------------------------------------------


## Figure 4
plot.FUI_RE <- function(r, 
                        fig = NULL, 
                        animal = NULL, # highlighted animal
                        align, 
                        col = "black",
                        col_fixed = "red",
                        Hz, 
                        var_name = NULL, 
                        title = NULL,
                        y_val_lim = 1.1,
                        ylim = NULL,
                        y_scal_orig = 0.05){
  library(gridExtra)
  name = NULL
  if(!is.null(fig))     fit_dat <- fig
  
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaHat[r,],
                             lower = fit_dat$betaHat[r,] - 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaHat[r,] + 2*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaHat[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaHat[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    #theme_bw() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = col_fixed), #",#ca0020" "#0868ac", "darkgray", "#E69F00", "#525252",
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  
  p.CI <- p.CI + 
    labs(x = "Trial Time-point (s)", y = "",#bquote(paste("Animal 2 Slope: ", beta[.(r-1)], "(s) + ", gamma[2][","][.(r-1)], "(s)")), 
         title = var_name[r]) +
    theme(legend.position = "none") #+
  
  if(is.null(animal)){
    RE <- data.frame(ID = rownames(fit_dat$randEff[,r,]), 
                     t(apply( as.matrix(fit_dat$randEff[,r,]), 1, function(x) x + as.numeric(fit_dat$betaHat[r,]) )))
    RE_g <- tidyr::pivot_longer(data = RE, cols = starts_with("X"), names_prefix = "X")
    p.CI <- p.CI +
      geom_line(data = RE_g, size = 0.5, 
                aes(colour = as.factor(ID), x = as.numeric(name) / Hz - align/Hz - 1/Hz, y = value)) +  #, colour = "#ca0020"
      scale_color_grey(start = 0.7, end = 0.1)
    
    p.CI <- p.CI +
      geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta), #",#ca0020" "#0868ac", "darkgray", "#E69F00", "#525252",
                data = beta.hat.plt, alpha = 1, size = 1, colour = col_fixed)
    
    }else{
      # one animal to highlight
    RE <- data.frame(ID = rownames(fit_dat$randEff[,r,]), 
                     t(apply( as.matrix(fit_dat$randEff[,r,]), 1, function(x) x + as.numeric(fit_dat$betaHat[r,]) )))
    RE_g <- tidyr::pivot_longer(data = RE, cols = starts_with("X"), names_prefix = "X")
    p.CI <- p.CI +
      geom_line(data = RE_g[RE_g$ID == unique(RE_g$ID)[animal],], size = 1, 
                aes(x = as.numeric(name) / Hz - align/Hz - 1/Hz, y = value), colour = col)
      
    p.CI <- p.CI +
      geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta), #",#ca0020" "#0868ac", "darkgray", "#E69F00", "#525252",
                data = beta.hat.plt, #alpha = 0.5, 
                size = 1, colour = col_fixed)
    }
  
  
  # plot random effects

  # make x and y intercepts
  if(!is.null(ylim)){
    p.CI <- p.CI + coord_cartesian(ylim = ylim)
  }else{  
    ylim <- c(min(beta.hat.plt$lower.joint), max(beta.hat.plt$upper.joint)) #layer_scales(p.CI)$y$range$range
    y_adjust <- y_scal_orig * (max(beta.hat.plt$upper.joint) - min(beta.hat.plt$lower.joint)) #layer_scales(p.CI)$y$range$range
    ylim[1] <- ylim[1] - y_adjust # just scale bottom because top is scaled below
  }     
  
  xlim = layer_scales(p.CI)$x$range$range
  
  x_range <- diff(xlim) * 0.1
  y_range <- diff(ylim) * 0.1
  y_range_up <- diff(ylim) * 0.02
  
  # extend upper limit
  y_val_lim <- c(1, y_val_lim)
  y_top <- (0.975) * diff(ylim*y_val_lim) + ylim[1]*y_val_lim[1]
  
  p.CI <- p.CI + 
    coord_cartesian(ylim = ylim*y_val_lim,
                    xlim = xlim) + 
    geom_segment(aes(x=0, xend=xlim[2] + x_range,
                     y=0,yend=0),
                 color = col_fixed, lwd=0.5, alpha = 0.5, linetype = "dashed")  # x - intercept

  return(p.CI)
}


# plot
fit_dat <- DA_re
fit_dat$betaHat[2,] <- b_iri # use fit from above for consistency
animal_highlight <- 8 # animal to highlight
col_highlight_fixed <- col_highlight #"#0868ac" 
col_highlight <- "#00A64F" 

# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI_RE(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = 1,
                                 animal = animal_highlight,
                                 col = col_highlight,
                                 col_fixed = col_highlight_fixed,
                                 var_name = c("Intercept", " "))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
                                       x_interval_list = list(c(-2,-0.5), c(-0.5, 1)),
                                       x_interval_text = NULL,
                                       text_list = list("Baseline", "Reward Period "),
                                       scl = 1.01, # percent for lines above original ylim values
                                       x_scl = 0.0025, # percent for gap in lines
                                       txt_adjust = 0.03, # percent text is above line
                                       txt_size = 3,
                                       col_seq = c("white", "white") #c("#ca0020", "#0868ac", "#E69F00", "#525252")
  )
}

# align plots
figs <- list()
figs[[1]] <- plot.f4[[1]]
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))

###########################
# make inset figure for p2
###########################

dd1 <- dd %>% 
  dplyr::filter(id == 1) %>% # highlighted animal
  as_tibble() %>%
  dplyr::mutate(Y = sim, 
                IRI = lick_time) %>%
  dplyr::select(Y, id, IRI) #%>%  # 

dd2 <- dd %>% 
  dplyr::filter(id != 1) %>% # highlighted animal
  as_tibble() %>%
  dplyr::mutate(Y = sim, 
                IRI = lick_time) %>%
  dplyr::select(Y, id, IRI)

corPlot3 <- 
  # remove IRI outliers for visualization
  ggplot(data = dd1, aes(x=IRI, y=Y), fill="white") +
  geom_abline(intercept = mod_coef[1], 
              slope = mod_coef[2], 
              color = col_highlight_fixed, size = 1#, 
              ) +
  geom_point(data = dd2,
             fill="white", size=1, 
             alpha=0.5, 
             color = dot_color) + #
  geom_point(fill="white", size=1, 
             alpha=0.75, 
             color = col_highlight) + #
  theme_classic() + 
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 0.75,
            color = col_highlight) +
  theme(axis.text = element_blank(),#element_text(size=4, face="bold"), 
        axis.title = element_text(size=8.5, face="bold"),
        plot.title = element_blank(),#element_text(hjust = 0.5, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        legend.title = element_text(face = "bold", size = 12)) +
  guides(fill = "none", alpha = "none", 
         color = "none") + 
  xlab("Latency") +
  ylab(latex2exp::TeX("\\textbf{Y(1.7):  $\\Delta F/F$ at 1.7sec}")) + #TeX("\\textbf{Y(0.1):}  \\Delta \\textbf{F / F Value at 0.1s}")) + #("Y(0.1):  Value at 0.1s") +
  annotate(geom = "text", x = 0.35, y = 1.3,
           label = latex2exp::TeX("\\textbf{Animal 2 Slope:}", bold = TRUE),
           angle = 46, # atan(beta1) * ( 180.0 / pi )
           color = col_highlight, size = 3) +
  annotate(geom = "text", x = 0.65, y = 1,
           label = latex2exp::TeX("$\\beta_1 (1.7) + \\gamma_{2,1} (1.7) = 0.95", bold = TRUE),
           angle = 46, # atan(beta1) * ( 180.0 / pi )
           color = col_highlight, size = 3) +
  annotate(geom = "text", x = 1.8, y = 0.875,
           label = latex2exp::TeX("$\\beta_1 (1.7)", bold = TRUE),
           angle = 46, # atan(beta1) * ( 180.0 / pi )
           color = col_highlight_fixed, size = 3#, 
           ) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5 ,5.5), "pt")) 



set.seed(25) # for jitter (position of highlighted dot)
p1_combined <- plot.f4[[1]] + geom_point(aes(x=time_point, y=b_val_2), size = 3, colour=col_highlight) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5 ,5.5), "pt")) +
  geom_segment(aes(x =time_point, y = b_val_2, xend = time_point, yend = -0.71),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75)  + 
  geom_segment(aes(x =-0.2, y = b_val_2, xend = time_point, yend = b_val_2),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75) +
  coord_cartesian(xlim = c(0, 4.5), 
                  ylim = c(-0.6, 1.85),  
                  clip = 'off') +
  scale_y_continuous(breaks = c(0, 1.5)) +#,
  scale_x_continuous(breaks = c(0, 4),
                     labels = c(0, 4)) +
  annotate(geom = "text", x = -0.5, y = b_val_2 * 1.01, 
           label =  "0.95",
           color = col_highlight, size = 4) +
  annotate(geom = "text", x = time_point, y = -0.8, 
           label =  "1.7",
           color = col_highlight, size = 4) +
  annotate(geom = "text", x = 3, y = -0.35,
           label = latex2exp::TeX("$\\beta_1 (s) + \\gamma_{2,1} (s)$"),
           color = col_highlight, size = 4) +
  annotate(geom = "text", x = 2.6, y = 0.125,
           label = latex2exp::TeX("$\\beta_1 (s)"),
           color = "black", size = 4, alpha = 0.75) +
  patchwork::inset_element(corPlot3, 0.57, 0.51, 1, 1, align_to = 'full')


# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration/Final")
ggsave( "FLMM Explanation RE Inset.pdf",
        plot = p1_combined,
        width = 4,
        height = 4)

# combine and save
l_combined[[2]] <- p1_combined
l_combined[[1]] + l_combined[[2]] 

setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration/Final")
ggsave( "FLMM Explanation RE Inset Combined.pdf",
        width = 8,
        height = 4)
