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
# cued trials (0 and 2) 
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

base_DA <- rowMeans( photo[ , 1:(cue_onset-1)] ) # baseline DA before cue onset
cue_DA <- rowMeans( photo[ , cue_onset:post_cue_end] ) # cue DA 
reward_DA <- rowMeans( photo[ , reward_onset:reward_off] ) # reward DA 

# DA during periods relative to pre-cue baseline
DA_mat <- data.frame( cue_DA =    cue_DA - base_DA,
                      reward_DA = cue_DA - base_DA)

# join data together
dat <- as.data.frame( cbind(ids, ITIs, dat, DA_mat) )
 dat <- cbind(dat, photo) # join full photometry matrix to trial/animal design matrix
rm(iti_vec, ids, ITIs, DA_mat, photo, trial_num)


# figure lick+/lick-
###################
# plot mean
###################
Hz <- 100
min_time <- 1.5 # time when cue starts from start of trials

lickPlus_color <- "#ca0020"
lickMinus_color <- "grey40"
avg_color <- "#0868ac"
# session 8 is the latest session in which all animals are recorded
dat %>% 
  as_tibble() %>%
  dplyr::filter(stimState == 0, # no opto
                group == "control",
                seshID >= 7) %>%  # only look at trals before
  dplyr::group_by(seshID, ids) %>%
  dplyr::summarise(m = n() ) %>% print(n = Inf)

# average data
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

# amplitude difference
library(ggbrace)
peaks <- mean_data %>% group_by(lickState) %>% summarise(peak = max(photometry))


# prep lick +/- trials 300 - 600
p1 <- mean_data %>% 
  ggplot(aes(x = time / Hz - min_time, y = photometry, ymin=y_low, ymax=y_high,  fill = lickState)) +
  geom_ribbon(aes(), alpha=0.3, colour = NA) + 
  geom_line(size = 1, aes(colour = lickState)) + 
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Peak Amplitude (Lick+ vs. Lick-)") +
  coord_cartesian(ylim = c(-0.75, 4)) +
  scale_fill_manual(values = c(lickPlus_color, lickMinus_color) ) + 
  scale_color_manual(values = c(lickPlus_color, lickMinus_color) ) +
  theme_classic() + 
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_text(size=14, face="bold"), #element_blank(),
        axis.title.y = element_blank(),
        legend.position="none") + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none")  +
  # geom_brace(aes( c(1.5,1.7), c(peaks$peak[2], peaks$peak[1]), label="Peak \nDifference"), 
  #            inherit.data=FALSE, labelsize=4, rotate = 270) +
  # lick + label
  annotate(geom = "text", x = -0.5, y = 1.2, label = "Lick+", 
           color = lickPlus_color, size = 5) +
  # lick- label
  annotate(geom = "text", x = -0.5, y = 0.35, label = "Lick-", 
           color = lickMinus_color, size = 5) 

###########################
# make inset figure for p1
###########################
# find difference in average curve peaks between each animal
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                seshID == 8) %>%  # only look at trals before
  dplyr::group_by(time, lickState, ids) %>% # time
  dplyr::summarise(photometry = mean(value, na.rm = TRUE) ) %>%
  dplyr::group_by(lickState, ids) %>%
  dplyr::summarise(photometry = max(photometry, na.rm = TRUE),
                   lickState = ifelse(lickState == 1, "Lick+", "Lick-") ) %>%
  unique() # save only unique rows (since all others are the same)

mean_data$lickState <- as.factor(mean_data$lickState)
mean_data$lickState <- relevel(mean_data$lickState, ref = c("Lick-")) # relevel for proper colors
mean_data$fill_col <- mean_data$lickState
mean_data$fill_col <- ifelse(mean_data$ids == control[1], mean_data$lickState, "Lick0")
mean_data$fill_col <- as.factor(mean_data$fill_col)

# make boxplot for inset
p_inst <- mean_data %>%
  ggplot(aes(x=lickState, y=photometry, color = lickState), fill="white") +
  geom_boxplot() +
  scale_fill_manual(values = c(lickMinus_color, lickPlus_color, "white") ) + 
  scale_color_manual(values = c(lickMinus_color, lickPlus_color, "white") ) +
  geom_jitter(aes(color = lickState, fill = fill_col), size=3, alpha=0.9, shape = 21) + 
  theme_classic() + 
  theme(axis.text.x = element_text(size=14, face="bold"), 
        axis.text.y = element_text(size=12, face="bold"),
        axis.title.y = element_blank(),
        legend.position="none") +
  ylab("Photometry Signal") +
  xlab("")

set.seed(25) # for jitter (position of highlighted dot)
p1_combined <- p1 + 
                  geom_curve(aes(x = 2, y = peaks$peak[1] * 0.94, xend = 4.05, yend = 1.45),
                            color = lickPlus_color, size = 0.2, alpha = 0.75, curvature = 0.55, 
                            arrow = arrow(length = unit(0.03, "npc")) ) + 
                  geom_curve(aes(x = 1.98, y = peaks$peak[2], xend = 3, yend = peaks$peak[2]* 0.98),
                             color = lickMinus_color, size = 0.2, alpha = 0.75, curvature = 0, 
                             arrow = arrow(length = unit(0.03, "npc")) ) + 
                  patchwork::inset_element(p_inst, 0.565, 0.5, 0.86, 0.9, align_to = 'full') # 0.65, 0.5, 1, 0.9

#######################
# area under the curve
#######################
time_min <- 151/Hz - min_time # start shading at this photometry index
time_max <- 500/Hz - min_time# 320  # stop shading
time_cue <- 300/Hz - min_time# 320  # stop shading

peak_color <- "darkblue"
time_color <- "darkgreen" 
mean_color <- "darkred" 

mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                seshID == 8,
                ids == control[1],
                group == "control",
                lickState == 0) %>%  
  dplyr::group_by(time, lickState) %>%
  dplyr::mutate(time = as.numeric(time) / Hz - min_time) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE) )

plot.f <- list()
mean_data$time <- as.numeric(mean_data$time)
mean_data$lickState <- as.factor(mean_data$lickState)

# prep lick +/- trials 300 - 600
p2 <-  mean_data %>% 
  ggplot(aes(x = time, y = photometry)) + #, / Hz - min_time
  geom_line(size = 1, colour = "black") + 
  coord_cartesian(ylim = c(-0.75, 4)) +
  geom_area(mapping = aes(x = ifelse(time<=time_max & time>=time_min, time, 0), fill = mean_color ), #data=subset(mean_data, time>time_min  & time<time_max), #aes(ymax=photometry), ymin=0,
            fill = mean_color, alpha=0.35) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Summary Statistics") +
  scale_fill_manual(values = c("black", "darkgray","#525252", "#E69F00", "#0868ac", "lightgrey" ) ) +
  theme_classic() + 
  geom_vline(xintercept=time_min,
             linetype="dashed",
             color = "darkgray",
             size = rel(0.3),
             alpha = 0.9) +
  geom_vline(xintercept=time_max,
             linetype="dashed",
             color = "darkgray",
             size = rel(0.3),
             alpha = 0.9) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.55) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_text(size=14, face="bold"), 
        axis.title.y = element_text(size=14, face="bold"),
        legend.position="bottom" ) + 
  guides(color= guide_legend(title="Prep Lick"), fill = "none") +
  # peak amplitude
  geom_segment(aes(x = -0.15, y = 0, xend = -0.15, yend = peaks$peak[2]),
               lineend = "round", linejoin = "bevel", colour = peak_color, 
               size = 0.5, 
              arrow = arrow(length = unit(0.5, "cm")) ) + 
  annotate(geom = "text", x = -0.35, y = peaks$peak[2] / 2 + 0.075, 
           label = "Peak Amplitude", fontface =2,
           color = peak_color, size = 3.35,  angle = 90) +
  # time to peak
  geom_segment(aes(x = 0.05, y = peaks$peak[2], xend = 1.8, yend = peaks$peak[2]),
               lineend = "round", linejoin = "bevel", colour = time_color,
               size = 0.5, 
               arrow = arrow(length = unit(0.5, "cm")) ) + 
  annotate(geom = "text", x = (1.825 - 0.3) / 2, y = peaks$peak[2] + 0.15, 
           label = "Time to Peak", fontface =2,
           color = time_color, size = 3.35) +
# time to peak
  annotate(geom = "text", x = 2.3, y = 0.4, label = "AUC", 
           color = mean_color, size = 3.35, fontface =2)
  

##########################################################################################################################################
#################################################################################
# variability within animal 1 across trials within individual session (and within lick+)
#################################################################################
# filter data
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                seshID == max(seshID), # last session so well trained
                ids == control[1], # arbitrarily choose the first control animal
                lickState == 1) %>%  # only look at lick+
  dplyr::mutate(time = as.numeric(time) / Hz - min_time) %>%
  dplyr::group_by(time) 

# average over trials
mean_dat1 <- mean_data  %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) # +/- sem
                   )

# filter and choose specific trials
mean_dat2 <- mean_data %>% 
  dplyr::filter(trial_num %in% 
                  as.integer(  seq( min(mean_data$trial_num),  
                                    max(mean_data$trial_num), 
                                    length = 10)) ) %>%
  dplyr::mutate(trial_num = trial_num - min(trial_num) + 1) %>%
  dplyr::rename(photometry = value)

# first plot average for animal 1 +/- s.e.m.
p4 <- 
  ggplot(NULL, aes(x = time, y = photometry)) + 
  geom_ribbon(data = mean_dat1, aes(ymin=y_low, ymax=y_high, fill = avg_color), alpha=0.5, colour = NA) + 
  geom_line(data = mean_dat2, size = 1, aes(colour = as.factor(trial_num))) + 
  scale_color_grey(start = 0.7, end = 0.1) +
  geom_line(data = mean_dat1, size = 1.5, colour = avg_color) + 
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Between-Trial Variability for Example Animal") +
  coord_cartesian(ylim = c(-2.5, 4)) +
  scale_fill_manual(values = c(avg_color) ) + 
  theme_classic() + 
  geom_vline(xintercept=151/ Hz - min_time,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  geom_vline(xintercept=300/ Hz - min_time, 
             linetype="dashed", 
             color = "gray", 
             size = rel(0.3),
             alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size=14, face="bold"),
        legend.position="right") + 
  guides(color= guide_legend(title="Trial"), fill = "none") +
  annotate(geom = "text", x = -0.75, y = 2, 
           label = latex2exp::TeX("Session Average", bold = TRUE),
           color = avg_color, size = 4)

#################################################################################
# variability within animal 1 across sessions (and within lick+)
#################################################################################
# last 5 sessions
# filter data average over trials to get Session Avg.s
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                seshID >= 8, # last common session (that all animals have) so well trained
                ids == control[1], # arbitrarily choose the first control animal
                lickState == 1) %>%  # only look at lick+
  dplyr::mutate(time = as.numeric(time) / Hz - min_time) %>%
  dplyr::group_by(time, seshID) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) # +/- sem
                  )

p5 <- mean_data %>%
  ggplot(aes(x = time, y = photometry, fill = as.factor(seshID) )) + 
  geom_line(size = 1, aes(colour = as.factor(seshID))) +  
  scale_color_grey(start = 0.7, end = 0.1) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Variability in Session Averages for Example Animal") +
  coord_cartesian(ylim = c(-.7, 2.25)) +
  theme_classic() + 
  geom_vline(xintercept=151/ Hz - min_time,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  geom_vline(xintercept=300/ Hz - min_time, 
             linetype="dashed", 
             color = "gray", 
             size = rel(0.3),
             alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position="right") + 
  guides(color= guide_legend(title="Session"), fill = "none")


#################################################################################
# variability for fixed session across animals (and within lick+)
#################################################################################
# pick session: find session with most number of animals thats late in training
sess_info <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                time == 10, # arbitrarily choose one time point on trial
                lickState == 1) %>%  # only look at lick+
  dplyr::group_by(seshID) %>%
  dplyr::summarise(n = n() )
sess_max <- sess_info$seshID[ which.max(sess_info$n) ] # chooses session 8 which is sufficiently late

# find animals with enough trials on this session (at least 10 (all have 60-90 except one animal that has 1 trial ))
id_sess_info <- dat %>% 
                  as_tibble() %>%
                  pivot_longer(cols = starts_with("photometry."),
                               names_to = "time",
                               names_prefix = "photometry.",
                               values_to = "value") %>% 
                  dplyr::filter(stimState == 0, # no opto
                                group == "control", # control animals
                                time == 10, # arbitrarily choose 
                                seshID == sess_max,
                                lickState == 1) %>%  # only look at lick+
                  dplyr::group_by(ids) %>%
                  dplyr::summarise(n = n() )
ids_inc <- id_sess_info$ids[ id_sess_info$n >= 10 ]

# filter data average over trials to get Session Avg.s for each animal
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                seshID == sess_max, # last session so well trained
                ids %in% ids_inc, # only animals with sufficient trials (only eliminates one animal)
                lickState == 1) %>%  # only look at lick+
  dplyr::mutate(time = as.numeric(time) / Hz - min_time,
                ids = as.integer(as.factor(ids))) %>%
  dplyr::group_by(time, ids) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE),
                   y_low = mean(value, na.rm = TRUE) - sd(value) / sqrt(n()), # +/- sem
                   y_high = mean(value, na.rm = TRUE) + sd(value) / sqrt(n()) )

# first plot average for animal 1 +/- s.e.m.
p6 <- mean_data %>%
  ggplot(aes(x = time, y = photometry, ymin=y_low, ymax=y_high, fill = as.factor(ids) )) + 
  geom_ribbon(aes(fill = avg_color), alpha=0.3, colour = NA) + 
  geom_line(size = 1, aes(colour = as.factor(ids))) + 
  scale_color_grey(start = 0.7, end = 0.1) +
  scale_fill_grey(start = 0.7, end = 0.1) +
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Variability in Session Averages Across Animals") +
  coord_cartesian(ylim = c(-0.75, 3.25)) +
  theme_classic() + 
  geom_vline(xintercept=151/ Hz - min_time,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  geom_vline(xintercept=300/ Hz - min_time, 
             linetype="dashed", 
             color = "gray", 
             size = rel(0.3),
             alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_text(size=14, face="bold"),
        legend.position="right") + 
  guides(color= guide_legend(title="Animal"), fill = "none")


#################################################################################
# variability within session across animals comparing Lick+/Lick-
#################################################################################
# pick session: find session with most number of animals thats late in training
sess_info <- dat %>% 
  as_tibble() %>%
  dplyr::filter(stimState == 0, # no opto
                group == "control" # control animals
                ) %>%  
  dplyr::group_by(seshID) %>%
  dplyr::summarise(n = n() )
sess_max <- sess_info$seshID[ which.max(sess_info$n) ] # chooses session 8 which is sufficiently late

# find animals with enough trials on this session (at least 10 (all have 60-90 except one animal that has 1 trial ))
id_sess_info <- dat %>% 
  as_tibble() %>%
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                seshID == sess_max) %>%  # only look at lick+
  dplyr::group_by(ids) %>%
  dplyr::summarise(n = n() )
ids_inc <- id_sess_info$ids[ id_sess_info$n >= 10 ]

# filter data average over trials to get Session Avg.s for each animal (Lick+/Lick-)
mean_data <- dat %>% 
  as_tibble() %>%
  pivot_longer(cols = starts_with("photometry."),
               names_to = "time",
               names_prefix = "photometry.",
               values_to = "value") %>% 
  dplyr::filter(stimState == 0, # no opto
                group == "control", # control animals
                seshID == sess_max, # last session so well trained
                ids %in% ids_inc # only animals with sufficient trials (only eliminates one animal)
                ) %>%  
  dplyr::mutate(time = as.numeric(time) / Hz - min_time,
                ids = as.integer(as.factor(ids))) %>%
  dplyr::group_by(time, ids, lickState) %>%
  dplyr::summarise(photometry = mean(value, na.rm = TRUE) )

# calculate differences in mean activity between lick+ and Lick-
mean_dataPlus <- mean_data %>% 
                    dplyr::filter(lickState == 1) %>%
                    dplyr::group_by(ids, time)

mean_dataMinus <- mean_data %>% 
                    dplyr::filter(lickState == 0) %>%
                    dplyr::group_by(ids, time)

mean_data <- mean_dataPlus # arbitrarily copy for dimension reasons
mean_data$photometry <- mean_dataPlus$photometry - mean_dataMinus$photometry # difference in means for each animal
rm(mean_dataPlus, mean_dataMinus)

# calculate mean of means (across animals)
mean_data_avg <- mean_data %>% 
                    dplyr::group_by(time) %>%
                    dplyr::summarise(y_low = mean(photometry, na.rm = TRUE) - sd(photometry, na.rm = TRUE) / sqrt(n()), # +/- sem
                                     y_high = mean(photometry, na.rm = TRUE) + sd(photometry, na.rm = TRUE) / sqrt(n()),
                                     photometry = mean(photometry, na.rm = TRUE))


# first plot average for animal 1 +/- s.e.m.
p7 <- 
  ggplot(NULL, aes(x = time, y = photometry)) + 
  geom_ribbon(data = mean_data_avg, aes(ymin=y_low, ymax=y_high, fill = avg_color), alpha=0.3) + 
  geom_line(data = mean_data, size = 1, aes(colour = as.factor(ids))) +  
  scale_color_grey(start = 0.7, end = 0.1) +
  geom_line(data = mean_data_avg, size = 1.5, colour = avg_color) + 
  xlab("Time from Cue Onset (sec)") +
  ylab("Photometry Signal") +
  labs(title = "Variability between Animals in Lick+/Lick- Difference") +
  coord_cartesian(ylim = c(-2.1, 2.25)) +
  scale_fill_manual(values = c(avg_color) ) + 
  theme_classic() + 
  geom_vline(xintercept=151/ Hz - min_time,
             linetype="dashed",
             color = "gray",
             size = rel(0.3),
             alpha = 0.7) +
  geom_vline(xintercept=300/ Hz - min_time, 
             linetype="dashed", 
             color = "gray", 
             size = rel(0.3),
             alpha = 0.7) +
  geom_hline(yintercept=0, 
             linetype="dashed", 
             color = "black", 
             size = rel(0.5),
             alpha = 0.7) +
  theme(plot.title = element_text(size=14, face="bold", hjust=0.5),
        axis.title.x = element_blank(), 
        axis.title.y = element_blank(),
        legend.position="right") + 
  guides(color= guide_legend(title="Animal"), fill = "none")  +
  annotate(geom = "text", x = -0.75, y = 1, 
           label = latex2exp::TeX("Average", bold = TRUE),
           color = avg_color, size = 4)
####################################################################
##########################
# join figures and save
##########################

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

library(patchwork)
set.seed(25)
plot_combined <- p6 + p7 + p4 + p5 + p2 + p1_combined + 
                    plot_layout(ncol = 2,
                                nrow = 3,
                                byrow = NULL,
                                widths = NULL,
                                heights = 10,
                                guides = NULL)

set.seed(25)
ggsave( "~/Desktop/NIMH Research/Photometry/fLME_methods_paper/Figures/intro_figure/intro_figure.pdf",
        plot = plot_combined,
        width = 14,
        height = 12)

