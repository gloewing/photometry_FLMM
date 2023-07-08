library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(refund)
library(dplyr)

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

# join data together
dat <- as.data.frame( cbind(ids, ITIs, dat) )
dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(iti_vec, ids, ITIs, photo, trial_num)

###################
# plot mean
###################
Hz <- 100
min_time <- 1.5 # time when cue starts from start of trials

lickPlus_color <- "#ca0020"
lickMinus_color <- "darkgray"

mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control",
                seshID == 8,
                ids == control[1]) %>%  # only look at trals before
  dplyr::group_by(time, lickState) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()), # +/- sem
                   lickState = ifelse(lickState == 1, "Lick+", "Lick-") )

mean_data$time <- as.numeric(mean_data$time)
mean_data$lickState <- as.factor(mean_data$lickState)
mean_data$lickState <- relevel(mean_data$lickState, ref = c("Lick+")) # relevel for proper colors

t_vec <- c(2.5, 4.5)
#####################################################
# method explanation figure
#####################################################
# raw photometry signal
p1 <- mean_data %>% 
  dplyr::filter(lickState == "Lick-") %>%
  ggplot(aes(x = time / Hz , y = photometry,  fill = lickState)) + # 
  geom_line(size = 1, aes(colour = lickState)) + 
  xlab("Time") +
  ylab("Photometry") +
  coord_cartesian(ylim = c(-0.75, 4)) +
  scale_fill_manual(values = c(lickPlus_color, lickMinus_color) ) + 
  scale_color_manual(values = c(lickPlus_color, lickMinus_color) ) +
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size=12, face="bold"), 
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none") #+


########################################################################################

p1 <- mean_data %>% 
  dplyr::filter(lickState == "Lick-") %>%
  ggplot(aes(x = time / Hz , y = photometry,  fill = lickState)) + # 
  geom_line(size = 1, aes(colour = lickState)) + 
  xlab("Time") +
  ylab("Photometry") +
  coord_cartesian(ylim = c(-0.75, 4)) +
  scale_fill_manual(values = c(lickPlus_color, lickMinus_color) ) + # 
  scale_color_manual(values = c(lickPlus_color, lickMinus_color) ) +
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_blank(),
        axis.title = element_text(size=12, face="bold"), #element_blank(),
        axis.text = element_blank(),
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none") 

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "method_photo_trace.pdf",
        plot = p1,
        width = 4,
        height = 2)


# methods
source("/Users/loewingergc/Desktop/NIMH Research/Photometry/code/existing_funs.R") # for time align functions
source("~/Desktop/NIMH Research/git_repos/fast_univ_fda/code/fui.R") # function code

dat2 <- dat[dat$seshID == 8,]

# photometry sample length
out_index <- grep(paste0("^", "photometry"), names(dat)) # indices that start with the outcome name
L <- length(out_index)
nknots_min <- round(L/2)
cov_idx <- seq(1:ncol(dat))[-out_index] # indices of covariates
out_idx <- seq(150, 600, length = 150) # down-sample
idx <- c(cov_idx, out_idx)
dat2$cueLatency <- scale(dat2$cueLatency)

DA <- fui(formula = photometry ~ lickState + latency +
                  (1| ids), 
                  data = dat2[,idx],  
                  family = "gaussian", 
                  silent = FALSE,
                  var = TRUE, 
                  analytic = TRUE,
                  parallel = TRUE,
                  G_return = TRUE,
                  smooth_method = "GCV.Cp",
                  splines = "tp",
                  subj_ID = "id")

Hz <- 30
fit_dat <- DA
set.seed(1)

# ylim of figures in original version (no noise)
ylim_mat <- matrix(NA, ncol = 2, nrow = 3)
ylim_mat[1,] <- c(-0.5120882,  3.4939680)
ylim_mat[2,] <- c(-0.7217547,  0.5873330)
ylim_mat[3,] <- c(-5.885771e-05,  1.604620e-05)

##### add noise to estimates for figures to accentuate effect of smoothing for explanatory purposes
for(r in 1:nrow(fit_dat$betaTilde)){
  # make "realistic" (correlated) noise using variance of beta hats
  
  Sig <- diag(0.01 * diag(fit_dat$betaHat.var[,,r]))
  noise_vec <- MASS::mvrnorm(n = 1, 
                             mu = rep(0, length(fit_dat$betaTilde[r,])),
                             Sigma = Sig)
  fit_dat$betaTilde[r,] <- fit_dat$betaTilde[r,] + noise_vec

}


## just estimates not smoothed
plot.FUI <- function(r, align = NULL, Hz, var_name = NULL, title = NULL){
  library(gridExtra)
  name = NULL

  if(is.null(var_name))    var_name <- rownames(fit_dat$betaTilde)
  if(nrow(fit_dat$betaTilde) != length(var_name) )  var_name <- rownames(fit_dat$betaTilde)
  if(is.null(align)) align <- 0
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaTilde), length.out = ncol(fit_dat$betaTilde)), 
                             beta = fit_dat$betaTilde[r,],
                             lower = fit_dat$betaTilde[r,] - 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaTilde[r,] + 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaTilde[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaTilde[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "white", alpha = 0) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "white", alpha = 1) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, size = 1) + # , lty = 5
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+

  if(r == 1){
    p.CI <- p.CI + labs(x = "", y = expression(beta[0](s)), title = var_name[r]) + #  
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") + # 
      theme(legend.title=element_blank(),
            legend.position = "None", # 
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            axis.text.x = element_blank(),
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "", y = bquote(paste(beta[.(r-1)], "(s)")), # 
                        title = var_name[r]) +
      theme(legend.position = "none",
            axis.text.x = element_blank(),) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") #+

  }
  
  return(p.CI)
}


# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, Hz = Hz, 
                               var_name = c("Intercept", "Covariate 1", "Covariate 2")) + 
                               coord_cartesian(ylim = ylim_mat[prdtr,])
}

fig1 <- do.call("grid.arrange", c(plot.f4, 
                                 # Hz, var_name = NULL,
                                 nrow = 1))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "method_fig1.pdf",
        plot = fig1,
        width = 8,
        height = 2)
#################################################
Hz <- 30
# fit_dat <- DA
# non-smoothed estimates with pointwise CIs
plot.FUI <- function(r, align = NULL, Hz, var_name = NULL, title = NULL){
  library(gridExtra)
  name = NULL
  
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  if(is.null(align)) align <- 0
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaTilde[r,],
                             lower = fit_dat$betaTilde[r,] - 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaTilde[r,] + 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaTilde[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaTilde[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "white", alpha = 0) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, size = 1) + # 
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+

  if(r == 1){
    p.CI <- p.CI + labs(x = "", y = expression(beta[0](s)), title = var_name[r]) + # 
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") + # 
      theme(legend.title=element_blank(),
            legend.position = "None", # 
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            axis.text.x = element_blank(),
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "", y = bquote(paste(beta[.(r-1)], "(s)")), #
                        title = var_name[r]) +
      theme(legend.position = "none",
            axis.text.x = element_blank()) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") #+

  }
  
  return(p.CI)
}


# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, Hz = Hz, 
                               var_name = c("", "", "")) + 
    coord_cartesian(ylim = ylim_mat[prdtr,])
}

fig1 <- do.call("grid.arrange", c(plot.f4, 
                                  nrow = 1))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "method_fig2a.pdf",
        plot = fig1,
        width = 8,
        height = 2)

#################################################

# smoothed estimates with pointwise CIs
fit_dat <- DA
plot.FUI <- function(r, align = NULL, Hz, var_name = NULL, title = NULL){
  library(gridExtra)
  name = NULL
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  if(is.null(align)) align <- 0
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaHat[r,],
                             lower = fit_dat$betaHat[r,] - 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaHat[r,] + 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaHat[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaHat[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "white", alpha = 0) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, size = 1) + #
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+

  if(r == 1){
    p.CI <- p.CI + labs(x = "", y = expression(beta[0](s)), title = var_name[r]) + 
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") + # 
      theme(legend.title=element_blank(),
            legend.position = "None", # 
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            axis.text.x = element_blank(),
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "", y = bquote(paste(beta[.(r-1)], "(s)")), #
                        title = var_name[r]) +
      theme(legend.position = "none",
            axis.text.x = element_blank()) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") #+

  }
  
  return(p.CI)
}


# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, Hz = Hz, 
                               var_name = c("", "", "")) + 
    coord_cartesian(ylim = ylim_mat[prdtr,])
}

fig2 <- do.call("grid.arrange", c(plot.f4, 
                                  nrow = 1))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "method_fig2.pdf",
        plot = fig2,
        width = 8,
        height = 2)


# smoothed estimates with joint CIs
plot.FUI <- function(r, align = NULL, Hz, var_name = NULL, title = NULL){
  library(gridExtra)
  name = NULL
  
  if(is.null(var_name))    var_name <- rownames(fit_dat$betaHat)
  if(nrow(fit_dat$betaHat) != length(var_name) )  var_name <- rownames(fit_dat$betaHat)
  if(is.null(align)) align <- 0
  decimal <- c(2,2,2,2,3)
  beta.hat.plt <- data.frame(s = seq(1, ncol(fit_dat$betaHat), length.out = ncol(fit_dat$betaHat)), 
                             beta = fit_dat$betaHat[r,],
                             lower = fit_dat$betaHat[r,] - 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper = fit_dat$betaHat[r,] + 1.5*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             lower.joint = fit_dat$betaHat[r,] - fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])),
                             upper.joint = fit_dat$betaHat[r,] + fit_dat$qn[r]*sqrt(diag(fit_dat$betaHat.var[,,r])))
  p.CI <- ggplot() +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper.joint, ymin = lower.joint), data = beta.hat.plt, fill = "gray20", alpha = 0.2) +
    geom_ribbon(aes(x = s / Hz - align/Hz - 1/Hz, ymax = upper, ymin = lower), data = beta.hat.plt, fill = "gray10", alpha = 0.4) +
    geom_line(aes(x = s / Hz - align/Hz - 1/Hz, y = beta, color = "Estimate"), data = beta.hat.plt, alpha = 1, size = 1) + #
    scale_colour_manual(name="", values=c("Estimate"="blue3")) +
    scale_y_continuous(labels=function(x) sprintf(paste0("%.", decimal[r], "f"), x)) #+

  if(r == 1){
    p.CI <- p.CI + labs(x = "Time (s)", y = expression(beta[0](s)), title = var_name[r]) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") + # 
      theme(legend.title=element_blank(),
            legend.position = "None", # 
            legend.justification = c("left", "top"),
            legend.box.just = "right",
            axis.text.x = element_blank(),
            legend.margin = margin(0, 0, 0, 0),
            legend.background = element_rect(fill=alpha('white', 0)))
  }else{
    p.CI <- p.CI + labs(x = "Time (s)", y = bquote(paste(beta[.(r-1)], "(s)")), 
                        title = var_name[r]) +
      theme(legend.position = "none",
            axis.text.x = element_blank()) +
      geom_hline(yintercept=0, color = "black", lwd=0.5, linetype = "dotted") #+

  }
  
  return(p.CI)
}


# plot
plot.f4 <- list()
for(prdtr in 1:nrow(fit_dat$betaHat)){
  plot.f4[[prdtr]] <- plot.FUI(r = prdtr, Hz = Hz, 
                               var_name = c("", "", "")) + 
    coord_cartesian(ylim = ylim_mat[prdtr,])
}

fig3 <- do.call("grid.arrange", c(plot.f4, 
                                  nrow = 1))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "method_fig3.pdf",
        plot = fig3,
        width = 8,
        height = 2)

##############################################
m <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control",
                lickState == 0,
                ids == control[1]) %>%  # only look at trals before
  dplyr::group_by(time, lickState, seshID) %>% 
  dplyr::summarise(photometry = mean(value, na.rm = TRUE))

m$seshID <- as.factor(m$seshID)
m$time <- as.numeric(m$time) / 100
m$photometry <- as.numeric(m$photometry)
m <- m[order(m$time),]

# add in a dot
t = 251 / 100 # time when put dot
photo_vals <- m %>%
  dplyr::filter(time == t) %>%
  dplyr::mutate(p1 = photometry) %>%
  dplyr::select(time, seshID, p1)

# join original dataset and dot dataset
m <- left_join(m, photo_vals, by = c("time", "seshID"))
line_dat <- data.frame(y2=t, x1=0) # for horizontal line
m <- cbind(m, line_dat)

################
# 3D - plot
################
library(plotly)
fig <- plot_ly(as.data.frame(m), x = ~photometry, y = ~time, z = ~seshID, 
               type = 'scatter3d', mode = 'lines', color = ~seshID, line = list(width = 4)) 


fig <- fig %>% 
        layout(scene = list(xaxis = list(title = 'Photometry', showgrid = F, showline = T, autotick = F, tickmode = "array", tickvals = c(0,2,4)),
                            yaxis = list(title = 'Time', showgrid = F, showline = T, autotick = F, tickmode = "array", tickvals = c(0,2,4,6)),
                            zaxis = list(title = 'Trial', showgrid = F, showline = T, showticklabels = FALSE)),
               showlegend = FALSE) %>%
        add_trace(x = ~p1, y = ~y2, z = ~seshID, mode = "markers",
                  marker = list(size = 3,
                  line = list(color = 'rgba(152, 0, 0, .8)',
                              width = 1)),
                  symbols = c('o')) %>% 
        add_lines(x = ~x1, y = ~y2, z = ~seshID, 
                  data = as.data.frame(m),
                  color = "grey",
                  line = list(color= "grey", widthh=0.5, dash="dash"))


# Remove the confidence interval
f1 <- ggplot(mtcars, aes(x=wt, y=scale(mpg))) + 
  geom_point(size = 3)+
  theme_classic() +
  theme(axis.title = element_text(size=22, face="bold"),
        axis.text = element_blank()) +
  geom_smooth(method=lm, se=FALSE, color = "gray50", size = 3) +
  xlab("X") +
  ylab(expression(Y(s)))

# save
setwd("/Users/loewingergc/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/method illustration")
ggsave( "regression_example.pdf",
        plot = f1,
        width = 3,
        height = 3)

######################
# table - methods
######################
library(lubridate)
library(xtable)

# functional x Longitudinal
tab_df <- data.frame(matrix( c("Longitudinal Functional (fGLMM)", "Generalized Linear Mixed-Effects (GLMM)",
                       "Functional Regression", "Generalized Linear Models (GLM)"), nrow = 2))

rownames(tab_df) <- c("Functional", "Scalar")
colnames(tab_df) <- c("Longitudinal", "Cross-Sectional")

xtable::xtable(tab_df)

# functional x functional
tab_df <- data.frame(matrix( c("Function-on-Function Regression (FoFR)", "Function-on-Scalar Regression (FoSR)",
                               "Scalar-on-Function Regression (SoFR)", "Generalized Linear Models (GLM)"), nrow = 2))

rownames(tab_df) <- c("Functional Predictors", "Scalar Predictors")
colnames(tab_df) <- c("Functional Outcome", "Scalar Outcome")

xtable::xtable(tab_df)
