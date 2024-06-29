# Lick aligned Figure 4 - Include Inset to explain figure
library(data.table)
library(dplyr)
library(mgcv)
library(lme4)
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
data_path <-"/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/simulations/science_paper" # path to original Matlab files for data

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
target_Hz_orig <- target_Hz <- 25
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
dat <- dat[complete.cases(dat$delay),]
dat[,-1] <- apply(dat[,-1], 2, as.numeric) # all except id
rm(dat_lst)
#########################################################################################
# creat AUC variables for model comparison
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



# trial trim to reduce computational burden
dat$id <- as.numeric(as.factor(dat$id)) # make 1:n

# number of knots
L <- ncol(dat[, (grep(paste0("^", "photometry"), colnames(dat)))])
nknots_min <- round(L/4) # L = 100

# order data
dat_photo <- dat[order(dat$id, dat$session, dat$trial, decreasing = FALSE), ] #%>%
n <- length(unique(dat_photo$id)) # number of subjects in (potentially) reduced sample

dat_photo$trial <- scale(dat_photo$trial) # important for stability of model

#########################################################################################

######################################################

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
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), # , color = "Estimate"
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"="black")) + # "blue3"
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  
  p.CI <- p.CI + 
    labs(x = "Time from Cue Onset (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
         title = var_name[r]) +
    theme(legend.position = "none") #+
  
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
    geom_segment(aes(x=xlim[1] - x_range, xend=xlim[2] + x_range,
                     y=0,yend=0),
                 color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed") + # x - intercept
    geom_segment(aes(y=ylim[1] - y_range, yend=y_top,#ylim[2]+y_range_up,
                     x=0,xend=0), # don't extend up
                 color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed") 
  
  #ggplot2::layer_scales(p.CI)$y$range$range <- ylim*y_val_lim
  
  return(p.CI)
}

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


# plot
plot.f4 <- list()
for(prdtr in 2:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr-1]] <- plot.FUI(r = prdtr, 
                                 Hz = target_Hz, 
                                 align = tm,
                                 var_name = c("Mean Signal on Short Delay Trials", "Mean Signal Difference: Long-Short"))#,
  
  # add interval names
  plot.f4[[prdtr-1]] <- interval_label(fig = plot.f4[[prdtr-1]],
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
plot.f <- plot.f4
fig2 <- do.call("grid.arrange", c(plot.f4, nrow = 1))


# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G/Final/Figure 5 Final")
ggsave( "exp3_fig4g_delay_rand_slope_original.pdf",
        plot = fig2,
        width = 4,
        height = 4)


#----------------------
# plot Inset
#----------------------
ylim = list(c(-0.25, 10.5), c(-2.25, 4.5))
# plot
#####################################################
library(ggbrace)

dat <- dat2
example_time <- 9.4

# trial trim to reduce computational burden
dat$id <- as.numeric(as.factor(dat$id)) # make 1:n

# number of knots
L <- ncol(dat[, (grep(paste0("^", "photometry"), colnames(dat)))])
nknots_min <- round(L/2) # L = 100

# order data
dat_photo <- dat[order(dat$id, dat$session, dat$trial, decreasing = FALSE), ] #%>%
n <- length(unique(dat_photo$id)) # number of subjects in (potentially) reduced sample
dat_photo$trial <- scale(dat_photo$trial) # important for stability of model

###############
# inset
###############
col_highlight <- "darkgreen"

# original data
mean_data <- dat_photo %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::mutate(time = round(as.numeric(time) / target_Hz - pre_min_tm,2),
                cue = delay) %>%
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

# peaks
peaks <- c(mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Long"],
           mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Short"])

p2 <- 
  ggplot(NULL, aes(x = time, y = photometry)) +      
  geom_ribbon(data = mean_dat0, aes(ymin=y_low, ymax=y_high, fill = as.factor(cue)), alpha=0.4) + 
  geom_line(data = mean_dat0, aes(colour = as.factor(cue)), size = 0.75) +  
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Original Data Session Averages") +
  scale_colour_manual(values = c("#ca0020", "#0868ac") ) + 
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.25),
             alpha = 0.7) +
  theme(plot.title = element_blank(),#element_text(size=14, face="bold", hjust=0.5),
        axis.text.y = element_text(size=5, face="bold"), 
        axis.text.x = element_text(size=6, face="bold"), 
        axis.title = element_text(size=7, face="bold"),
        legend.position = "none",
        strip.text.x = element_text(face="bold", color="black", size = rel(1))) + 
  guides(color= guide_legend(title="Delay"), fill = "none") + 
  coord_cartesian(ylim = c(-0.25, 2.75),
                  xlim = c(7.5,11.5)) + # , clip = "off"
  # bars over top
  # # vertical dotted lines
  geom_segment(aes(x=example_time,xend=example_time,
                   yend=mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Long"],#mean_data$value[mean_data$time == example_time],
                   y=mean_dat0$photometry[mean_dat0$time == example_time & mean_dat0$cue == "Short"]), 
               color = col_highlight, linetype="dashed") +
  # # geom_segment(aes(x=2,xend=2,y=0,yend=6), color = "grey 20", linetype="dashed") +
  annotate(geom="text", x=8.1, y=1.1, label="Long Average", color = "#ca0020", fontface =2, size = 2.5) +
  annotate(geom="text", x=8.1, y=-0.2, label="Short Average", color = "#0868ac", fontface =2, size = 2.5) +
  geom_point(aes(x=example_time, y=peaks[1]), 
             size = 1.5, colour="#ca0020") +
  geom_point(aes(x=example_time, y=peaks[2]), 
             size = 1.5, colour="#0868ac") +
  scale_x_continuous(breaks = c(8, 9.4, 11),
                     labels = c(8, 9.4, 11)) +
  # brace
  geom_brace(aes(c(9.55, 9.85), c(peaks[2], peaks[1]), label=" "), 
             labelsize=3.5, rotate = 90) + # inherit.data=FALSE, 
annotate(geom="text", x=10.5, y=0.85, 
         label=latex2exp::TeX("$\\beta_1 (9.4) = 1.6", bold = TRUE),
         color = col_highlight, fontface =2, size = 2.5) +
  annotate(geom="text", x=10.5, y=1.3, 
           label="Long - Short:",
           color = col_highlight, fontface =1, size = 2.5) 
  

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
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank()) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), # , color = "Estimate"
              data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"="black")) + # "blue3"
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+
  
  p.CI <- p.CI + 
    labs(x = "Time from Cue Onset (s)", 
         title = var_name[r]) +
    theme(legend.position = "none",
          axis.title.y = element_blank()) #+
  
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
                 color = "black", lwd=0.5, alpha = 0.75, linetype = "dashed") #+ # x - intercept
  geom_segment(aes(y=ylim[1], yend=y_top,#ylim[2]+y_range_up,
                   x=0,xend=0), # don't extend up
               color = "white", lwd=0.5, alpha = 0, linetype = "dashed")
  
  return(p.CI)
}

###############
# zoomed in 
###############
# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, 
                               Hz = target_Hz, 
                               align = tm,
                               var_name = c("", ""))#,
  
  # add interval names
  plot.f4[[prdtr]] <- interval_label(fig = plot.f4[[prdtr]],
                                     x_interval_list = list(c(-1,0), c(0, 2)),
                                     x_interval_text = list(c(-1.25), c(1.75)), # 
                                     text_list = list("", ""),
                                     scl = 1.01, # percent for lines above original ylim values
                                     x_scl = 0.0025, # percent for gap in lines
                                     txt_adjust = 0.03, # percent text is above line
                                     txt_size = 3,
                                     col_seq = c("white", "white") )
}

# align plots
fig <- do.call("grid.arrange", c(plot.f4, nrow = 1))


b_val <- fit_dat$betaHat[2, (example_time + pre_min_tm) * target_Hz]
set.seed(25) # for jitter (position of highlighted dot)
p1_combined <- 
  plot.f4[[2]] + 
  # fig +
  geom_point(aes(x=example_time, y=b_val), size = 3, colour=col_highlight) +
  theme(plot.margin = unit(c(5.5, 5.5, 5.5 ,25), "pt")) +
  geom_segment(aes(x =example_time, y = b_val, xend = example_time, yend = -2.4),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75)  + 
  geom_segment(aes(x =-2.6, y = b_val, xend = example_time, yend = b_val),
               colour = col_highlight,
               linetype = "dashed",
               size = 0.75) +
  coord_cartesian(xlim = plot.f4[[2]]$coordinates$limits$x, 
                  ylim = plot.f4[[2]]$coordinates$limits$y,  clip = 'off') + # 
  scale_x_continuous(breaks = c(0, 5),
                     labels = c(0, 5)) +
  annotate(geom = "text", x = -3.5, y = b_val * 0.96, 
           label =  round(b_val,1),
           color = col_highlight, size = 4) +
  annotate(geom = "text", x = example_time-0.2, y = -2.5, # y= -2.7
           label =  example_time,
           color = col_highlight, size = 4) +
  geom_segment(aes(x =example_time, y = b_val, xend = 8.5, yend = b_val* 1.5),
               lineend = "round", linejoin = "bevel", colour = col_highlight, # grey35
               size = 0.5, 
               arrow = arrow(length = unit(0.3, "cm")) ) + #, ends = "both"
  patchwork::inset_element(p2, 0, 0.61, 0.65, 1, align_to = 'full')


library(patchwork)
plot.f[[1]] + p1_combined + 
  plot_layout(ncol = 2,
              nrow = 1,
              byrow = NULL,
              widths = NULL,
              guides = NULL)

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/Science-Fig4G/Final/Figure 5 Final")
ggsave( "Delay Effect Inset Combined.pdf",
        width = 8,
        height = 4)


