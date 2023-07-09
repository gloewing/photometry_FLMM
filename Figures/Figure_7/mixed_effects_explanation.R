library(R.matlab)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(ggbrace)
library(latex2exp)

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

# photometry
photo <- do.call(rbind, matdat[photo_idx,,] )
colnames(photo) <- paste0("photometry.", 1:ncol(photo) )

base_DA <- rowMeans( photo[ , 1:(cue_onset-1)] ) # baseline DA before cue onset
cue_DA <- rowMeans( photo[ , cue_onset:post_cue_end] ) # cue DA 
reward_DA <- rowMeans( photo[ , reward_onset:reward_off] ) # reward DA 

# DA during periods relative to pre-cue baseline
DA_mat <- data.frame( cue_DA =    cue_DA - base_DA,
                      reward_DA = cue_DA - base_DA)

# join data together
dat <- as.data.frame( cbind(ids, dat, DA_mat) )
dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(ids, DA_mat, photo, trial_num)


###################
# plot 1
###################
Hz <- 100
min_time <- 1.5 # time when cue starts from start of trials

lickPlus_color <- "#ca0020"
lickMinus_color <- "#0868ac"
id_num <- 3 # 2  1 not bad

mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control",
                seshID == 8,
                ids == control[id_num]) %>%  # only look at trals before
  dplyr::group_by(time, lickState) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - 2*sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + 2*sd(value) / sqrt(n()), # +/- sem
                   lickState = ifelse(lickState == 1, "Lick+", "Lick-") )

mean_data$time <- as.numeric(mean_data$time)
mean_data$lickState <- as.factor(mean_data$lickState)
mean_data$lickState <- relevel(mean_data$lickState, ref = c("Lick+")) # relevel for proper colors

ls <- 0 # lick state # 0

dat_filt <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control",
                seshID == 7, # 7
                ids == control[id_num], # 2
                lickState == ls # "Lick+"
                )

# arbitrarily take first 20 trials
trial_samp <- unique(dat_filt$trial_num)[1:10] # arbitrarily take sample of 20 trials
dat_filt <- dat_filt %>% dplyr::filter(trial_num %in% trial_samp)

# mean of each condition
mean_data <- dat_filt %>%  # only look at trials before
  dplyr::group_by(time) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()), # +/- sem
                   sem = sd(value) / sqrt(n)) %>%
  unique()

dat_filt <- dat_filt %>% 
  dplyr::mutate(photometry = value) %>%
  dplyr::select(time, photometry, trial_num)

# change variable types
mean_data$time <- as.numeric(mean_data$time)
dat_filt$time <- as.numeric(dat_filt$time)
dat_filt$lickState <- as.factor(1) # use as indicator of mean
mean_data$lickState <- as.factor(0)

p1 <- #%>% 
  ggplot(NULL, aes(x = as.numeric(time) / Hz - min_time, y = photometry)) +
  geom_ribbon(data = mean_data, aes(ymin=y_low, ymax=y_high,  fill = lickPlus_color), 
              alpha=0.3, colour = NA) + 
  geom_line(data = mean_data, size = 1, aes(colour = lickPlus_color)) + 
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  geom_line(data = dat_filt, size = 0.75, alpha = 0.55, 
            aes(colour = as.factor(lickPlus_color), group = as.factor(trial_num))) + 
  labs(title = "Session Average: Photometry vs. Behavior") +
  coord_cartesian(ylim = c(-2, 9)) +
  scale_fill_manual(values = c(lickPlus_color, lickMinus_color) ) + # 
  scale_color_manual(values = c(lickPlus_color, lickMinus_color) ) +
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=14, face="bold"), 
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none")  +
  # lick + label
  annotate(geom = "text", x = 0.16, y = 5.475, 
           label = "paste(bold(Animal~1:~Session~Average))", 
           parse = TRUE,
           alpha = 1,
           color = lickPlus_color, size = 4) +
  # lick- label
  geom_segment(aes(x = 1.9, y = 5.1, xend = 2.9, yend = 5.475), 
             color = lickPlus_color, size = 0.4, alpha = 1, 
             arrow = arrow(length = unit(0.03, "npc")) )
###########
# inset
###########
# simualte behavioral data

# find peaks of data
# random effects parameters
beta <- 0.5 # mean beta
re_sd <- 0.1 # random effects sd 0.05
int_mean <- 0 # intercept # no random intercept
int_sd <- 0.1 # random intercept sd
noise_sd <- 0.25 # error sd

# cue and reward period 
cue_idx <- 100:400
peaks <- dat_filt %>% group_by(trial_num) %>% summarise(peak = max(photometry[cue_idx]))

# simulate other animals
set.seed(26) # 3, 7
n <- 10 # trials per subject
K <- 5 # subjects
beta_vec <- rnorm(K, mean = beta, sd = re_sd)
rand_int <- rnorm(K, mean = int_mean, sd = int_sd)
n_vec <- c(10, sample(5:15, K - 1, replace = TRUE)) # vector of sample sizes
X <- lapply(1:(K-1), function(x) runif(n_vec[x+1], min = (x-0.25), max = (x+1.25) ) ) 
X <- append( list(peaks$peak), X) # add in subject photometry data
x_mean <- sapply(X, mean) # mean of X for each subject
y <- lapply(1:K, function(j) X[[j]] * beta_vec[j] + 
              rand_int[j] + # intercept
              rnorm(n_vec[j], mean = 0, sd = noise_sd) ) # error term
y_mean <- sapply(y, mean) # mean of y for each subject
sim_dat <- data.frame(x = do.call(c, X), 
                      y = do.call(c, y), 
                      id = do.call(c, sapply(1:K, function(x) rep(x, each = n_vec[x]) ) ) 
                      )
sim_dat$id <- as.factor(sim_dat$id)
sim_dat$lickState <- as.factor(ifelse(sim_dat$id == 1, 0, 1)) # use this for coloring

# dataframe of averages for inset
mean_df <- data.frame(y = y_mean, x = x_mean, id = as.factor(1:K))
mean_df$lickState <- as.factor(ifelse(mean_df$id == 1, 0, 1)) # use this for coloring

# make boxplot for inset
p_inst <- mean_df %>%
  ggplot(aes(x=x, y=y), fill="white") +
  geom_point(aes(fill = id, colour = id), size=3, alpha=0.9, shape = 21) +  
  scale_fill_manual(values = c(lickPlus_color, rep("white", 4)) ) + 
  scale_color_manual(values = c(lickPlus_color, "#525252","#E69F00","darkgreen", lickMinus_color) ) +
  theme_classic() + 
  theme(axis.text = element_text(size=6, face="bold"), 
        axis.title = element_text(size=8, face="bold"),
        legend.position = c(.4, .95), 
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),  
        legend.title = element_text(face = "bold", size = 8),
        legend.key.height = unit(0.001, 'cm')
        ) +
  guides(fill = "none", alpha = "none", color = guide_legend(title="Animal")) +
  ylab("Session Avg. Behavior") +
  xlab("Session Avg. Photometry") + 
  geom_smooth(method=lm, se=FALSE, color = "black", size = 1, linetype = 5)


set.seed(25) # for jitter (position of highlighted dot)
p1_combined <- p1 + 
  patchwork::inset_element(p_inst, 0.65, 0.5, 1, 0.9, align_to = 'full') #+

p1_combined
#######################################################
# trial by trial
#######################################################
# find minimum value that is easy to point to on inset
peaks2 <- dat_filt %>% 
  group_by(trial_num) %>% 
  summarise(peak = max(photometry[1:250]))

min_idx <- peaks2$trial_num[which.max(peaks2$peak)]
dat_filt$min_value <- 0.5
dat_filt$min_value[dat_filt$trial_num == min_idx] <- 1
dat_filt$sizeL <- 0.75
dat_filt$sizeL[dat_filt$trial_num == min_idx] <- 1
x_val <- peaks$peak[peaks$trial_num == min_idx]


p2 <- 
  ggplot(NULL, aes(x = as.numeric(time) / Hz - min_time, y = photometry)) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  geom_line(data = dat_filt, 
            aes(colour = as.factor(lickState), group = as.factor(trial_num), alpha = min_value), size = 0.75) + 
  labs(title = "Trial-by-Trial: Photometry vs. Behavior") +
  geom_line(data = dat_filt[dat_filt$trial_num == min_idx,], 
            aes(colour = as.factor(lickState)), alpha = 1, size = 1) + 
  coord_cartesian(ylim = c(-2, 9)) +
  scale_fill_manual(values = c(lickPlus_color, lickMinus_color) ) + 
  scale_color_manual(values = c(lickPlus_color, lickMinus_color) ) +
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title = element_text(size=14, face="bold"), 
        legend.position="none")   +
  # lick- label
  annotate(geom = "text", x = 0.3, y = 5.5, 
           label = "paste(bold(Animal~1:~Trial~j))", 
           parse = TRUE,
           color = lickPlus_color, size = 4) +
  geom_segment(aes(x = 1.9, y = 5.82, xend = 2.9, yend = 6.2), 
             color = lickPlus_color, size = 0.4, alpha = 1, 
             arrow = arrow(length = unit(0.03, "npc")) ) +
  scale_alpha_continuous(range=range(dat_filt$min_value)) 


min_x <- x_val 
sim_dat$alphaVec <- ifelse(sim_dat$id == 1, 1, 0.65)
sim_dat$highlight <- as.factor(ifelse(sim_dat$id == 1 & sim_dat$x == min_x, 0, 1))

# make inset
p_inst2 <- sim_dat %>%
  ggplot(aes(x=x, y=y)) + # 
  geom_point(aes(alpha = alphaVec, fill = highlight, colour = id), size=3, shape = 21) + #  
  scale_fill_manual(values = c(lickPlus_color, rep("white", 4)) ) + 
  scale_color_manual(values = c(lickPlus_color, "#525252","#E69F00","darkgreen", lickMinus_color) ) +
  theme_classic() + 
  theme(axis.text = element_text(size=6, face="bold"), 
        axis.title = element_text(size=8, face="bold"),
        legend.position = c(1, .47),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6),  
        legend.title = element_text(face = "bold", size = 8),
        legend.key.height = unit(0.001, 'cm')
        ) + 
  guides(fill = "none", alpha = "none", color = guide_legend(title="Animal")) +
  ylab("Trial Behavior") +
  xlab("Trial Peak Photometry") + 
  geom_smooth(method=lm, se=FALSE, color = "black", size = 1, linetype = 5) + 
  scale_alpha_continuous(range=range(sim_dat$alphaVec))  +
  coord_cartesian(ylim = c(0, 2.75),
                  xlim = c(1, 8))

p2_combined <- p2 + 
  patchwork::inset_element(p_inst2, 0.65, 0.5, 1, 0.9, align_to = 'full') #+

#######################################################

library(latex2exp)
corPlot2 <- 
  sim_dat %>%
  ggplot(aes(x=x, y=y), fill="white") +
  geom_point(aes(color = id),  fill="white", size=3, alpha=0.9, shape = 21) +
  scale_fill_manual(values = c(lickPlus_color, rep("white", 4)) ) + 
  scale_color_manual(values = c(lickPlus_color, "#525252","#E69F00","darkgreen", lickMinus_color) ) +
  theme_classic() + 
  geom_line(stat="smooth",method = "lm", formula = y ~ x,
            size = 1,
            linetype = 5,
            color = "black",
            alpha = 1) +
  geom_smooth(aes(color = id), method=lm, se=FALSE, size = 1.25) +
  theme(axis.text = element_text(size=12, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        legend.title = element_text(face = "bold", size = 12),
        ) +
  guides(fill = "none", alpha = "none", 
         color = guide_legend(title="Animal", 
                              override.aes = list(linetype = 0) # removes lines from dots in legend 
                              )) +
  xlab("Trial Peak Photometry") +
  ylab("Trial Behavior")  +
  coord_cartesian(ylim = c(0.3, 2.75),
                  xlim = c(1, 7.75)) +
  geom_segment(aes(x = 2.71, y = 2.45, xend = 2.35, yend = 1.33),
             color = "#E69F00", size = 0.2, alpha = 1, 
             arrow = arrow(length = unit(0.03, "npc")) ) +
  annotate(geom = "text", x = 1.8, y = 2.5, label = "paste(bold(Animal~3~OLS~Fit))", parse = TRUE,
           color = "#E69F00", size = 4.25) +
  annotate(geom = "text", x = 1.75, y = 2.37, label = TeX("Parameters: $\\hat{\\beta}^{(OLS)}_3$"),
           color = "#E69F00", size = 3.75) +
# ols population
  annotate(geom = "text", x = 6.59, y = 2.8, label = "paste(bold(All~Animals~OLS~Fit))", parse = TRUE,
           color = "black", size = 4.25) +
  annotate(geom = "text", x = 6.365, y = 2.68, label = TeX("Parameters: $\\hat{\\beta}^{(OLS)}$"),
           color = "black", size = 3.75) #+

#######################################################
# LME
#######################################################
# lme model
mod <- lme4::lmer(y~x + (x | id),
                  data = sim_dat)
m=summary(mod)
beta_fixed <- m$coefficients[,1] # fixed effects
sim_dat$fixed <- beta_fixed[1] + sim_dat$x * beta_fixed[2]
sim_dat$fitted <- fitted(mod) # fitted values with BLUPs

# animal 3 fitted values
m3 <- lm(y~x, data = sim_dat[sim_dat$id == 3,])
m3 <- m3$fitted.values 
sim_dat$fitted3 <- NA # NAs for all other animals
sim_dat$fitted3[sim_dat$id == 3] <- m3

corPlot3 <- 
  sim_dat %>%
  ggplot(aes(x=x, y=y), fill="white") +
  geom_point(aes(color = id),  fill="white", size=3, alpha=0.9, shape = 21) +
  scale_fill_manual(values = c(lickPlus_color, rep("white", 4)) ) + 
  scale_color_manual(values = c(lickPlus_color, "#525252","#E69F00","darkgreen", lickMinus_color) ) +
  theme_classic() + 
  # main OLS line
  geom_line(stat="smooth", method = "lm", formula = y ~ x,
            size = 1,
            linetype = 5,
            color = "black",
            alpha = 0.5) +
  # fixed effects of lme
  geom_line(aes(x = x, y = fixed),
            size = 1, 
            color = "black",
            linetype = 5) +
  # mixed model individual curves
  geom_line(aes(group = id, color = id, x = x, y = fitted),
            size = 1.5) +
  # individual ols
  geom_line(aes(color = id, x = x, y = fitted3),
            size = 1.5,
            alpha = 0.5) +
  theme(axis.text = element_text(size=12, face="bold"), 
        axis.title = element_text(size=14, face="bold"),
        legend.position = c(1, .5),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6), 
        legend.title = element_text(face = "bold", size = 12),
  ) +
  guides(fill = "none", alpha = "none", 
         color = guide_legend(title="Animal", 
                              override.aes = list(linetype = 0) # removes lines from dots in legend 
         )) +
  xlab("Trial Peak Photometry") +
  ylab("Trial Behavior") +
  coord_cartesian(ylim = c(0.3, 2.75),
                  xlim = c(1, 7.75)) +
  geom_curve(aes(x = 2.5, y = 2.3, xend = 2.7, yend = 1.4),
             color = "#E69F00", size = 0.2, alpha = 1, curvature = 0.32, # 0.5
             arrow = arrow(length = unit(0.03, "npc")) ) +
  annotate(geom = "text", x = 1.85, y = 2.5, label = "paste(bold(Animal~3~LMM~Fit))", parse = TRUE,
           color = "#E69F00", size = 4.25) +
  annotate(geom = "text",  x = 2, y = 2.37, label = "paste(Parameters:~~~~~~~~~~~+~hat(gamma)[3])", parse=TRUE,
           color = "#E69F00", size = 3.75) +
  annotate(geom = "text", x = 2.43, y = 2.382, label = TeX("$\\hat{\\beta}^{(LMM)}$"),
           color = "black", size = 3.75) +
  # LME
  annotate(geom = "text", x = 6.69, y = 2.8, label = "paste(bold(LMM~Fixed-Effects~Fit))", parse = TRUE,
           color = "black", size = 4.25) +
  annotate(geom = "text", x = 6.365, y = 2.68, label = TeX("Parameters: $\\hat{\\beta}^{(LMM)}$"),
           color = "black", size = 3.75)

#######################################################
library(patchwork)
set.seed(25)
plot_combined <- p1_combined + p2_combined + corPlot2 + corPlot3 + 
  plot_layout(ncol = 2,
              nrow = 2,
              byrow = NULL,
              widths = NULL,
              heights = 6,
              guides = NULL)

set.seed(25)
setwd("~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/mixed_effects_explanation/Final")
ggsave( "mixed_effects_explan.pdf",
        plot = plot_combined,
        width = 12,
        height = 12)
