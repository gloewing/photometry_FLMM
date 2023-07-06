# simulate data from fLME - multiple random effects
photo_stimulate <- function(X = NULL, # design matrix
                            Z = NULL, # design matrix for random effects
                            N = NULL, # target sample size
                            include = NULL, # extra data to include
                            model, # fLME model object
                            fixed_smooth = FALSE, # whether to use smoothed betas for true fixed effects
                            IDs,  # col of IDs (N x 1 vector)
                            resid_var_indep = TRUE, # assume residual variance is indepedent across fn. domain as in paper
                            resid_var_subj = TRUE, # calculate error variance in subject-specific fashion as opposed to globally
                            resid_var_scale = 1, # scale the error variance
                            fixed_effects_draw = FALSE, # if use uncertainty around fixed effects
                            fixed_effects_scale = 1, # scale non-intercept of fixed effects
                            random_effects_scale = 1, # scale non-intercept random effects
                            rand_int = FALSE, # only use a random intecept (0 out and other random effects)
                            rand_int_sigma = FALSE, # whether to use random intercept Sigma for random slopes covariance matrix
                            ID_level_re = TRUE, # indicates that random effects are subject-level only (do not vary within subject)
                            fixed_shift = 0, # amount to shift fixed effect functional coefficient
                            fixed_shift_idx = 2 # index of fixed effect coefficient to shift
                            ){
  
  # draw fixed effects
  if(fixed_smooth)   fixed_eff <- model$betaTilde # use non-smoothed betas
  if(!fixed_smooth)  fixed_eff <- model$betaHat  # use smoothed betas
  
  fixed_eff[-1,] <- fixed_eff[-1,] * fixed_effects_scale
    
  # use model-saved design matrix
  if(is.null(X))   X <- model$designmat
  if(!is.null(include))   include_name <- colnames(include)
  
  # use model-saved design matrix
  if(is.null(Z))   Z <- model$Z
  RE_dims <- sapply(Z, ncol) # total number of random effect design matrix columns for each random effect terms
  qq <- sum(RE_dims) # total number of random effect columns
  q <- length(Z) # actual number of random effects
  
  IDs <- IDs_orig <- as.numeric(as.factor(IDs)) # ensure all IDs are 1:n sequence
  beta <- fixed_eff # initialize true fixed effect vector
  L <- ncol(beta) # dim of function
  p <- nrow(beta) # number of covaraites *including intercept*
  
  # sample size
  subjects <- subjects_orig <- sort( unlist(unique(IDs)) )
  n <- length(subjects) # sample size of observed data
  if(is.null(N))    N <- n
  
  # order design matrix based on IDs
  X <- X[order(IDs, decreasing = FALSE), ]
  Z <- lapply(Z, function(x) x[order(IDs, decreasing = FALSE),])
  
  # draw synthetic subject IDs from observed IDs and use the corresponding design matrix
  if(N != n){
    # subj_samp <- sample(subjects, N, replace = I(N > n) ) # sample subject IDs
    if(N > n){
      subj_samp <- sample(subjects, N - n, replace = TRUE ) # sample subject IDs
      subj_samp <- c(subjects, subj_samp)
    }else{
      subj_samp <- sample(subjects, N, replace = FALSE ) # sample subject IDs
    }
    idx_list <- lapply(subjects, function(x) which(IDs == x)) # find row indices corresponding to each ID
    samp_lst <- lapply(subj_samp, function(x) idx_list[[ x ]]) # element j has row indices of original dataset to be used for synthetic subject j 
    idx <- do.call(c, samp_lst ) # concatenate
    # n <- N
    X <- X[idx,] # use IDs to repeat the design matrix
    Z <- lapply(Z, function(x) x[idx,]) # NEED TO CHANGE
    include <- include[idx,]
    # IDs <- IDs[idx]
    subjects <- 1:n
    source_IDs <- IDs[idx] # each row has which observed source animal
    IDs <- do.call(c, sapply(seq_along(samp_lst), function(x) rep(x, length(samp_lst[[x]]))  ) ) # new IDs just repeated in order because of how samp_lst is structured
    
    rm(idx)
  }else{
    # if n == N then subj_samp is just the original subjects
    subj_samp <- subjects
    source_IDs <- IDs
  }
  
  Y <- matrix(0, nrow = nrow(X), ncol = L) # initialize outcome matrix
  colnames(Y) <- paste0("photometry.", 1:L)
  nn <- nrow(X)
  
  if(fixed_effects_draw){
    # draw true fixed effects from distribution that accounts for uncertainty in model estimates provided to simulator
    for(l in 1:p)   beta[l,] <- MASS::mvrnorm(n = 1, mu = fixed_eff[l,], Sigma = model$betaHat.var[,,l])
  }else{
    beta <- fixed_eff
  }
  
  # shift fixed effect coefficient up or down
  beta[fixed_shift_idx,] <- beta[fixed_shift_idx,] + fixed_shift
  
  # residual variance across functional domain from lfosr function
  # resid_var <- diag( model$var_random["var.Residual",] ) # --not smoothed use diagonal matrix of residual variance from lfosr function (so errors will be iid)
  resid_var <- diag( as.numeric( model$R ) ) # smoothed use diagonal matrix of residual variance from lfosr function (so errors will be iid)
  resid_cov <- NULL #diag(as.numeric(resid_var)) # this will be replaced below depending on resid_var_indep and resid_var_subj
  
  ## *** below is for non-independent ***
  if(!resid_var_indep & !resid_var_subj){
    # if you don't want to assume error are independent across functional domain
    # use residuals to calculate full covariance matrix
    resid_cov <- cov(model$residuals)
    # resid_cov <- refund::fbps(resid_cov)$Yhat # smooth - dont smooth because issues with 
    
    # trim eigenvalues to ensure matrix is PSD
    edcomp <- eigen(resid_cov) ## trim non-positive eigenvalues to ensure positive semidefinite
    eigen.positive <- which(edcomp$values > 0)
    
    if(length(eigen.positive) == qq){
      # nothing needed here because matrix is already PSD
    }else 
    if(length(eigen.positive) == 0){
      resid_cov <- tcrossprod(edcomp$vectors[,1]) * edcomp$values[1] 
    }else{
      # sum of outerproducts of eigenvectors scaled by eigenvalues for all positive eigenvalues
      resid_cov <- Reduce('+', lapply(eigen.positive, function(x)  tcrossprod(edcomp$vectors[,x]) * edcomp$values[x] ) ) 
    }
    
    rm(eigen.positive, edcomp)
  }else if(resid_var_indep & !resid_var_subj){
    resid_cov <- as.numeric( matrixStats::colVars(as.matrix( model$residual) ) )
    argvals <- 1:L
    
    # smooth across function
    resid_cov <- diag(as.numeric(smooth.spline(x = argvals, y = resid_cov)$y))
  }
  
  random_effects <- array(0, c(qq, L, nn)) # initialize all random effects as 0
  
  ########################
  # draw random effects
  ########################
  #random_var <- model$var_random["var.ID.(Intercept)",] # unsmoothed
  # random_var <- as.numeric(model$H) # smoothed  -- change to as.matrix() if not just random intercept
  trim_flag <- 0
  RE_dims <- c(0, cumsum(RE_dims))
  var_names <- dimnames(model$var_random)[[1]][-length(dimnames(model$var_random)[[1]])] # make sure you are using the right indicies of GHat[]
  var_nm_idx <- grep(paste0("^", "var."), var_names)
  
  for(re in 1:qq){
    
    # find which actual RE term corresponds to the column of Z (with expanded factor levels)
    if(q == 1){
      re_trm <- 1
    }else{
      re_trm <- which( sapply( 1:length(RE_dims-1), function(x) re >= (RE_dims[x] + 1) & re <= RE_dims[x+1] ) )
    }
    
    flag <- ifelse(trim_flag == re_trm, TRUE, FALSE) # see if already trimmed
    flag <- ifelse(rand_int_sigma & re_trm > 1, TRUE, FALSE) # whether to re-use random intercept Sigma
    if(!flag){
      trim_flag <- re_trm # set trim_flag to current random effect 
      
      # trim Sigma
      Sigma <- as.matrix(model$G[ var_nm_idx[re_trm] ,,] )
      edcomp <- eigen(Sigma)
      eigen.positive <- which(edcomp$values > 0)
      
      # check to see if we need to trim eigenvalues
      if(length(eigen.positive) == L){
        # do not need to do anything
      }else if(length(eigen.positive) == 0){
        Sigma <- tcrossprod(edcomp$vectors[,1]) * edcomp$values[1] #matrix(edcomp$vectors[,1], ncol = 1) %*% edcomp$values[1] %*% matrix(edcomp$vectors[,1], nrow = 1)
      }else{
        # sum of outerproducts of eigenvectors scaled by eigenvalues for all positive eigenvalues
        Sigma <- Reduce('+', lapply(eigen.positive, function(x)  tcrossprod(edcomp$vectors[,x]) * edcomp$values[x] ) ) #tcrossprod(edcomp$vectors[,1]) * edcomp$values[1] #Outer(x, y, oper = "*") matrix(edcomp$vectors[,1], ncol = 1) %*% edcomp$values[1] %*% matrix(edcomp$vectors[,1], nrow = 1) 
      }
      
      # Sigma[,i,j] <- GTilde[,j,i] <- diag(cov.trimmed) # save trimmed version of G covariance matrix: q x T x T (the diagonal of cov.trimmed is the cov(u_l (s_i), u_l(s_j) ) -- so same random effect but different time points
    }
      
    
    # iterate through columns of random effects
    if(re_trm != 1)     Sigma <- Sigma * random_effects_scale # scale random effects (for non random intercept)
    if(!rand_int | re_trm == 1)   random_effects[re,,] <- t( MASS::mvrnorm(nn, mu = rep(0, L), Sigma = Sigma ) )   # draw random effect vector for each RE column of Z (not necesssarily)

  }

  # subj_samp <- sample(subjects, N, replace = TRUE) # sample subject IDs
  # idx_list <- lapply(subjects, function(x) which(IDs == x))
  # samp_lst <- lapply(subj_samp, function(x) idx_list[[ x ]])
  # idx <- do.call(c, samp_lst )
  # IDs_new <- do.call(c, sapply(seq_along(samp_lst), function(x) rep(x, length(samp_lst[[x]])) ) ) # new IDs
  # n <- N
  # X <- X[idx,]
  # IDs <- IDs[idx] 
  Z_i <- as.matrix(do.call(cbind, Z))
  
  for(i in 1:N){
    idx <- which( IDs == i ) # rows corresponding to this subject in new dataset
    idx_orig <- which(IDs_orig == subj_samp[i]) # subject rows of observed data ID (that synthetic ID corresponds to)
    beta_i <- beta 
    # Z_i <- do.call(cbind, lapply(Z, function(x) x[idx,]) ) # use rows of original data (Z) and concatenate design matrix
    
    # # calculate residual variance
    if(!resid_var_indep & resid_var_subj){
      # use residuals to calculate full covariance matrix -- assumes non-independent errors across functional domain
      resid_cov <- cov(model$residuals[idx_orig,])
      # resid_cov <- refund::fbps(resid_cov)$Yhat # smooth
      
      # browser()
      # trim eigenvalues to ensure matrix is PSD
      # edcomp <- eigen(resid_cov) ## trim non-positive eigenvalues to ensure positive semidefinite
      # eigen.positive <- which(edcomp$values > 0)
      # if(length(eigen.positive) != 2){
      #   resid_cov <- matrix(edcomp$vectors[,1], ncol = 1) %*% edcomp$values[1] %*% matrix(edcomp$vectors[,1], nrow = 1)
      #   rm(eigen.positive, edcomp)
      #   }
      
      
    }else if(resid_var_indep & resid_var_subj){
      resid_cov <- as.numeric( matrixStats::colVars(as.matrix( model$residual[idx_orig,]) ) )
      argvals <- 1:L
      
      # smooth across function
      resid_cov <- diag(as.numeric(smooth.spline(x = argvals, y = resid_cov)$y))
    }
        
    
    if(is.null(resid_cov)){
      resid_cov <- diag(as.numeric(resid_var))
      message("use residual variance")
    }   
    # browser()
    eps <- MASS::mvrnorm(length(idx), mu = rep(0, L), Sigma = resid_cov * resid_var_scale) # error
    
    # calculate Z_{i,j}^T gamma_{i,j,l} for each i (note: Z_{i,j} is same for all l \in [L] :: i.e., same across all indices on functional domain)
    # re_i <- sapply( seq_along(idx_orig), function(ii) Z_i[ii, ] %*% random_effects[,,idx_orig[ii] ] ) # each observation is multiplied by a different draw of a random effect
    if(ID_level_re){
      # only one random effect draw per animal (arbitrarily take the first random effect draw)
      re_i <- sapply( idx, function(ii) Z_i[ii, ] %*% random_effects[,,idx[1] ] ) # each observation is multiplied by a different draw of a random effect
      
    }else{
      re_i <- sapply( idx, function(ii) Z_i[ii, ] %*% random_effects[,,ii ] ) # each observation is multiplied by a different draw of a random effect
    }
    
    # draw sample of outcome across functional domain for all observations for subject i
    Y[idx,] <- as.matrix(X[idx,]) %*% beta_i + t(re_i) + eps
  }

  # remove intercept column and return data with column names 
  if(is.null(include)){
    dat <- cbind(IDs, X[,-1], Y)
    colnames(dat) <- c("ID", colnames(X)[-1], colnames(Y) )
    
  }else{
    dat <- cbind(IDs, include, X[,-1], Y) 
    colnames(dat) <- c("ID", include_name, colnames(X)[-1], colnames(Y) )
    
  }
  

  return( list(data = as.data.frame(dat),
               beta = beta,
               random_effects = random_effects,
               source_IDs = source_IDs, # indices of observed ID for each row
               subj_samp = subj_samp, # IDs of sampled subjects
               resid_cov = resid_cov # covariance for error
               )
          )
  
}



####################################
# function to save on cluster
####################################

saveFn <- function(file, fileNm, iterNum, save.folder = NA){
  
  if( !is.na(save.folder) ){
    # set working directory if specified
    fileNm <- paste0(save.folder, "/", fileNm)
  }
  
  # check if file exists
  if(  file.exists(fileNm)  ){
    # if exists read in file and save this result to correspond to row
    res <- read.csv( fileNm )
    res[iterNum,] <- file[iterNum,]
    write.csv(res, fileNm, row.names = FALSE)
    
  }else{
    # if it does not exist (first iteration to complete) then save resMat
    write.csv(file, fileNm, row.names = FALSE)
  }
  
}

# cluster finder for permutation test -- only significant for certain lengths
cluster_finder <- function(r, min.length = 3){
  # take in binary vector r and return binary vector that 0s out any sequence of 
  # 1s that is not of at least length min.length
  m <- rep(0, length(r))
  
  runs=rle(r > 0)
  runs.lengths.cumsum = cumsum(runs$lengths)
  myruns = which(runs$values == TRUE & runs$lengths >= min.length)
  ends = runs.lengths.cumsum[myruns]
  
  newindex = ifelse(myruns>1, myruns-1, 0)
  starts = runs.lengths.cumsum[newindex] + 1
  if (0 %in% newindex) starts = c(1,starts)
  
  if(length(starts) > 0){
    for(x in 1:length(starts)){  m[seq(starts[x], ends[x])] <- 1 }
  }
  
  return(m)
}


# trim eigenvalues and return PSD matrix
# V is a n x n matrix
eigenval_trim <- function(V){
  edcomp <- base::eigen(V) ## trim non-positive eigenvalues to ensure positive semidefinite
  eigen.positive <- which(edcomp$values > 0)
  
  if(length(eigen.positive) == ncol(V)){
    # nothing needed here because matrix is already PSD
  }else if(length(eigen.positive) == 0){
    V <- tcrossprod(edcomp$vectors[,1]) * edcomp$values[1] 
  }else{
    # sum of outerproducts of eigenvectors scaled by eigenvalues for all positive eigenvalues
    V <- Reduce('+', lapply(eigen.positive, function(x)  tcrossprod(edcomp$vectors[,x]) * edcomp$values[x] ) ) 
  }
  
  return(V)
}



####################################
# function to save on cluster
####################################

saveFn_new <- function(file, fileNm, iterNum, iters, save.folder = NA){
  
  filePrefix <- paste0(fileNm, "_")
  
  if( !is.na(save.folder) ){
    # set working directory if specified
    fileNmIndiv <- paste0(save.folder, "/", filePrefix, iterNum)
  }
  
  # save individual file
  write.csv(file[iterNum,], fileNmIndiv, row.names = FALSE)
  
  fls <- list.files(save.folder) # all files
  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
  
  if(length(fls_nms) == iters){
    # if all have been saved
    
    for(i in 1:iters){
      ff <- paste0(save.folder, "/", filePrefix, i)
      res <- read.csv( ff )
      
      # first one use is file
      if(i == 1){
        mat <- res
      }else{
        mat[i,] <- res[i,]
      }
      
      setwd(save.folder)
      
      base::unlink( paste0(filePrefix, iterNum) ) # delete file
    }
    
    write.csv(mat, fileNm, row.names = FALSE)
    
  }
  
  # # check if file exists
  # if(  file.exists(fileNm)  ){
  #   # if exists read in file and save this result to correspond to row
  #   res <- read.csv( fileNm )
  #   res[iterNum,] <- file[iterNum,]
  #   write.csv(res, fileNm, row.names = FALSE)
  #   
  # }
}




####################################
# function to save on cluster
####################################

process_files <- function(fileNm, save.folder = NA){
  
  if(is.na(save.folder))  save.folder <- getwd()
  
  filePrefix <- paste0(fileNm, "_")
  
  fls <- list.files(save.folder) # all files
  fls_nms <- grep(paste0("^", filePrefix), fls) # files with names that match
  setwd(save.folder)
  
  for(i in 1:length(fls_nms)){
    
      ff <- fls[ fls_nms[i] ]
      res <- read.csv( ff )
      num <- gsub(filePrefix, "" ,ff)
      
      # first one use is file
      if(i == 1){
        mat <- res
        base::unlink( ff ) # delete file
      }else{
        mat[i,] <- res[num,]
        base::unlink( ff ) # delete file
      }
  }
  
  write.csv(mat, fileNm, row.names = FALSE)
  
}
    
