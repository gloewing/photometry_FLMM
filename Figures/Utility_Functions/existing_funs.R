# functions for photometry data pre-processing
# written by Gabriel Loewinger


# time aligning
time_align <- function(time_truth, data, name = "timestamps", save_time = FALSE){
  # time_truth - is a vector of timestamps for the (usually photometry) timepoints to align to (ground truth)
  # data - is the dataset with timestamps to change to align to the time_truth
  # "name"  -  is the column name in "data" variable
  # save_time  - whether to save original time variable
  
  data_new <- data.table::as.data.table( data[name] ) 
  tm <- data.table::data.table(time_temp = as.numeric(time_truth) ) # vector of times we want to align with (photometry times) -- ground truth time
  tm[, time_aligned := time_truth]
  
  data.table::setkeyv(data_new, name) # column name in data file that we want to align to be consistent with photometry (the ground truth time)
  data.table::setkeyv(tm, c('time_temp'))
  
  data_new <- as.data.frame( tm[data, roll='nearest'] )
  
  # delete original time variable
  if(!save_time){
    data_new <- subset(data_new, select = -time_temp )
    colnames(data_new)[ colnames(data_new) == "time_aligned" ] <- name
  }   
  
  return(data_new)
  
}


########################################################
# behavior-bout finder across who session
########################################################
# finds all bouts (no trial times given)
bout_finder <- function(x, Hz = 1000, ibi = 1){
  # Hz - sampling rate
  # x - behavior time series (vector of times when )
  # ibi - inter behavior interval in seconds (length of time between any press/lick etc before bout considered over)
  
  time_diff <- diff( c(0, x) )
  bout_interval <- ibi * Hz # number of samples in bout length between any successive behavioral event
  
  # find bout start times/indices
  bout_start_idx <- which(time_diff > ibi) # find start of bout indices 
  bout_start_time <- x[bout_start_idx]
  
  # find bout end times/indices
  num_events <- bout_end_time <- vector(length = length(bout_start_idx))
  for(i in 1:(length(bout_start_time) - 1) ){
    bout_i_start <- bout_start_time[i]
    tm_diff <- time_diff[x > bout_i_start & x < bout_start_time[i+1] ] # all inter behavior intervals starting after ith bout but before i+1th bout
    tm <- x[x >= bout_i_start & x < bout_start_time[i+1] ] # same as about but the times
    num_events[i] <- length(tm)
    bout_end_time[i] <- tm[ length(tm) ] # last time in this vector of times (by definition all will be within correct ibi)
  }
  
  # last bout only has one
  if(is.na(num_events[ length(num_events) ]) ){
    num_events[ length(num_events) ] <- 1
    bout_end_time[ length(num_events) ] <- bout_start_time[ length(num_events) ]
  }
  
  # return results
  res <- matrix(NA, ncol = 4, nrow = length(bout_start_time))
  colnames(res) <- c("bout_start_time", "bout_end_time", "num_events", "rate")
  res[,1] <- bout_start_time
  res[,2] <- bout_end_time
  res[,3] <- num_events
  res[,4] <- num_events / (bout_end_time - bout_start_time)
  
  return(res)
}


########################################################
# behavior-bout finder for trials, must give bout_start_time
########################################################

bout_finder_trial <- function(x, Hz = 1000, ibi = 1, bout_start_time = NULL, include = NULL){
  # Hz - sampling rate
  # x - behavior time series (vector of times when )
  # ibi - inter behavior interval in seconds (length of time between any press/lick etc before bout considered over)
  # bout_start_time - if it is null then just find all bouts, otherwise find bouts starting at those times (e.g., if the bouts are during a reward period)
  # include - Boolean vector of length of bout_start_time, can be used to exlcude trials and include NAs
  
  if(is.null(include))    include <- rep(TRUE, length(bout_start_time))
  
  time_diff <- diff( c(0, x) )
  bout_interval <- ibi * Hz # number of samples in bout length between any successive behavioral event
  
  # find bout end times/indices
  num_events <- bout_end_time <- vector(length = length(bout_start_time))
  for(i in 1:(length(bout_start_time)) ){
    if(include[i]){
      bout_i_start <- bout_start_time[i]
      bout_next <- ifelse(i < length(bout_start_time), bout_start_time[i+1], max(x) ) # if last bout, then just use trial end as last time
      tm_diff <- time_diff[x >= bout_i_start & x < bout_next ] # all inter behavior intervals starting after ith bout but before i+1th bout
      tm <- x[x >= bout_i_start & x < bout_next ] # same as about but the times
      num_events[i] <- length(tm)
      bout_end_time[i] <- tm[ length(tm) ] # last time in this vector of times (by definition all will be within correct ibi)
    }
  }
  
  # last bout only has one
  if(is.na(num_events[ length(num_events) ]) ){
    num_events[ length(num_events) ] <- 1
    bout_end_time[ length(num_events) ] <- bout_start_time[ length(num_events) ]
  }
  
  # return results
  res <- matrix(NA, ncol = 4, nrow = length(bout_start_time))
  colnames(res) <- c("bout_start_time", "bout_end_time", "num_events", "rate")
  res[include,1] <- bout_start_time[include]
  res[include,2] <- bout_end_time[include]
  res[include,3] <- num_events[include]
  res[include,4] <- (num_events / (bout_end_time - bout_start_time))[include]
  
  return(res)
  
}
