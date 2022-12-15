library(forcats)
library(plyr)
library(dplyr)
library(simIReff)
library(VineCopula)
library(rio)
library(stringr)
library(doParallel)
library(ggplot2)
library(ggpubr)
library(gtable)
library(grid)
library(ggplotify)
library(tools)
library(devEMF)
library(boot)
library(plotly)
library(glue)
library(moments)
library(corrplot)
theme_set(theme_bw())

# ===================== Beautify names & Re-Order factors for plotting ===================== 
beautify <- function(x) {
  dict <- list()
  # measures
  dict['ap'] <- 'AP'
  dict['ndcg20'] <- 'nDCG@20'
  dict['err20'] <- 'ERR@20'
  dict['p10'] <- 'P@10'
  dict['rr'] <- 'RR'
  # margins
  dict['norm'] <- 'Tr.Norm'
  dict['beta'] <- 'Beta'
  dict['nks'] <- 'Tr.Norm KS'
  dict['bks'] <- 'Beta KS'
  dict['bbinom'] <- 'Beta-Binom'
  dict['dks(1)'] <- 'Disc. KS'
  dict['dks(2)'] <- 'Disc. KS (2)'
  dict['dks(5)'] <- 'Disc. KS (5)'
  dict['dks(10)'] <- 'Disc. KS (10)'
  # copulas
  dict['Gaussian'] <- 'Gaussian'
  dict['t'] <- 'Student t'
  dict['Clayton'] <- 'Clayton'
  dict['Survival Clayton'] <- 'Clayton'
  dict['Rotated Clayton 90 degrees'] <- 'Clayton'
  dict['Rotated Clayton 270 degrees'] <- 'Clayton'
  dict['Gumbel'] <- 'Gumbel'
  dict['Survival Gumbel'] <- 'Gumbel'
  dict['Rotated Gumbel 90 degrees'] <- 'Gumbel'
  dict['Rotated Gumbel 270 degrees'] <- 'Gumbel'
  dict['Frank'] <- 'Frank'
  dict['Joe'] <- 'Joe'
  dict['Survival Joe'] <- 'Joe'
  dict['Rotated Joe 90 degrees'] <- 'Joe'
  dict['Rotated Joe 270 degrees'] <- 'Joe'
  dict['BB1'] <- 'BB1'
  dict['Survival BB1'] <- 'BB1'
  dict['BB6'] <- 'BB6'
  dict['Survival BB6'] <- 'BB6'
  dict['BB7'] <- 'BB7'
  dict['Survival BB7'] <- 'BB7'
  dict['Rotated BB7 90 degrees'] <- 'BB7'
  dict['Rotated BB7 270 degrees'] <- 'BB7'
  dict['BB8'] <- 'BB8'
  dict['Survival BB8'] <- 'BB8'
  dict['Rotated BB8 90 degrees'] <- 'BB8'
  dict['Rotated BB8 270 degrees'] <- 'BB8'
  dict['Tawn  type 1'] <- 'Tawn'
  dict['Tawn  type 2'] <- 'Tawn'
  dict['Rotated Tawn type 1 90 degrees'] <- 'Tawn'
  dict['Rotated Tawn type 2 90 degrees'] <- 'Tawn'
  dict['Rotated Tawn type 1 180 degrees'] <- 'Tawn'
  dict['Rotated Tawn type 2 180 degrees'] <- 'Tawn'
  dict['Rotated Tawn type 1 270 degrees'] <- 'Tawn'
  dict['Rotated Tawn type 2 270 degrees'] <- 'Tawn'
  dict['Independence'] <- 'Indep.'
  return(dict[[x]])
}

beautify2 <- function(x) {
  dict <- list()
  dict['Tawn  type 1'] <- 'Tawn 1'
  dict['Tawn  type 2'] <- 'Tawn 2'
  return(dict[[x]])
}

recode_and_reorder_help <- function(x, ...) {
  suppressWarnings(fct_recode(fct_relevel(x, ...), ...))
}

recode_and_reorder_help_rev <- function(x, ...) {
  suppressWarnings(fct_rev(fct_recode(fct_relevel(x, ...), ...)))
}

recode_and_reorder <- function(column){
  return(recode_and_reorder_help(column,
                                 `1` = "1",
                                 `2` = "2",
                                 `3` = "3",
                                 `4` = "4",
                                 `5` = "5",
                                 `6` = "6",
                                 `7` = "7",
                                 `8` = "8",
                                 `9` = "9",
                                 `10` = "10",
                                 `11` = "11",
                                 `12` = "12",
                                 `13` = "13",
                                 `14` = "14",
                                 `15` = "15",
                                 `16` = "16",
                                 `17` = "17",
                                 `18` = "18",
                                 `19` = "19",
                                 `20` = "20",
                                 `21` = "21",
                                 `22` = "22",
                                 `23` = "23",
                                 `24` = "24",
                                 `25` = "25",
                                 # measures
                                 `AP` = "ap",
                                 `P@10` = "p10",
                                 `RR` = "rr",
                                 `nDCG@20` = "ndcg20",
                                 `ERR@20` = "err20",
                                 # margins
                                 `Tr.Norm` = "norm",
                                 `Beta`  = "beta",
                                 `Tr.Norm KS`  = "nks",
                                 `Beta KS`  = "bks",
                                 `Beta-Binom` = "bbinom",
                                 `Disc. KS` = "dks(1)",
                                 `Disc. KS (2)` = "dks(2)",
                                 `Disc. KS (5)` = "dks(5)",
                                 `Disc. KS (10)` = "dks(10)",
                                 # copulas
                                 `Gaussian` = 'Gaussian',
                                 `Student t` = 't',
                                 `Frank` = 'Frank',
                                 
                                 `Clayton` = 'Clayton',
                                 `Clayton` = 'Survival Clayton',
                                 `Clayton` = 'Rotated Clayton 90 degrees',
                                 `Clayton` = 'Rotated Clayton 180 degrees',
                                 `Clayton` = 'Rotated Clayton 270 degrees',
                                 
                                 `Gumbel` = 'Gumbel',
                                 `Gumbel` = 'Survival Gumbel',
                                 `Gumbel` = 'Rotated Gumbel 90 degrees',
                                 `Gumbel` = 'Rotated Gumbel 180 degrees',
                                 `Gumbel` = 'Rotated Gumbel 270 degrees',
                                 
                                 `Joe` = 'Joe',
                                 `Joe` = 'Survival Joe',
                                 `Joe` = 'Rotated Joe 90 degrees',
                                 `Joe` = 'Rotated Joe 180 degrees',
                                 `Joe` = 'Rotated Joe 270 degrees',
                                 
                                 `BB1` = 'BB1',
                                 `BB1` = 'Survival BB1',
                                 `BB1` = 'Rotated BB1 90 degrees',
                                 `BB1` = 'Rotated BB1 180 degrees',
                                 `BB1` = 'Rotated BB1 270 degrees',
                                 
                                 `BB6` = 'BB6',
                                 `BB6` = 'Survival BB6',
                                 `BB6` = 'Rotated BB6 90 degrees',
                                 `BB6` = 'Rotated BB6 180 degrees',
                                 `BB6` = 'Rotated BB6 270 degrees',
                                 
                                 `BB7` = 'BB7',
                                 `BB7` = 'Survival BB7',
                                 `BB7` = 'Rotated BB7 90 degrees',
                                 `BB7` = 'Rotated BB7 180 degrees',
                                 `BB7` = 'Rotated BB7 270 degrees',
                                 
                                 `BB8` = 'BB8',
                                 `BB8` = 'Survival BB8',
                                 `BB8` = 'Rotated BB8 90 degrees',
                                 `BB8` = 'Rotated BB8 180 degrees',
                                 `BB8` = 'Rotated BB8 270 degrees',
                                 
                                 `Tawn` = 'Tawn  type 1',
                                 `Tawn` = 'Tawn  type 2',
                                 `Tawn` = 'Rotated Tawn type 1 90 degrees',
                                 `Tawn` = 'Rotated Tawn type 2 90 degrees',
                                 `Tawn` = 'Rotated Tawn type 1 180 degrees',
                                 `Tawn` = 'Rotated Tawn type 2 180 degrees',
                                 `Tawn` = 'Rotated Tawn type 1 270 degrees',
                                 `Tawn` = 'Rotated Tawn type 2 270 degrees',
                                 `Indep.` = 'Independence',
                                 # criteria
                                 `SHC` = 'SHC',
                                 `LogLik` = 'LogLik',
                                 `BIC` = 'BIC',
                                 `AIC (top-2)` = 'AIC (top-2)',
                                 `AIC` = 'AIC',
                                 `AIC (only Beta KS)` = 'AIC (only Beta KS)',
                                 `Worst` = 'Worst',
                                 `Best` = 'Best',
                                 `Random` = 'Random',
                                 #
                                 `Only Beta KS` = 'Only Beta KS',
                                 `Except Beta KS` = 'Except Beta KS'
  ))
}


# ===================== Data reading ===================== 

# Remove the bottom p% of runs
removeBottom <- function(dat, p = .25) {
  mu <- colMeans(dat)
  i <- mu >= quantile(mu, p)
  dat[,i]
}

# Remove duplicate runs by their per-topic score
removeDuplicates <- function(dat, tol = 1e-5) {
  toremove <- integer(0)
  for(i in seq(ncol(dat)-1)) {
    x <- dat[,i]
    for(j in seq(i+1,ncol(dat))) {
      y <- dat[,j]
      if(all(abs(x-y)<=tol))
        toremove <- c(toremove, j)
    }
  }
  if(length(toremove) > 0)
    dat[-toremove]
  else
    dat
}

# Read our (.csv formatted) evaluation data
read_evaluation_data <- function(measures, collections, remove_bottom=0.1) {
  dat <- vector(mode = "list")
  for (measure in measures) {
    for(collection in collections) {
      collmeas <- paste0(collection, "_", measure)
      path_scores <- file.path("data", paste0(collmeas, ".csv"))
      if(file.exists(path_scores)) {
        d <- import(path_scores)
        d <- removeDuplicates(d)
        d <- removeBottom(d, p = remove_bottom) # remove bottom X% of systems
        dat[[measure]][[collection]] <- d
      }
    }
  }
  return(dat)
}

# Convert raw evaluation data to our simple .csv formal
read_terabyte2006 <- function(terabyte2006_path='data/terabyte2006/', path_out='data/') {
  for (measure in c('P10', 'recip_rank', 'map')) {
    measure_change_name <- function(measure) {
      dict <- list()
      dict['map'] <- 'ap'
      dict['P10'] <- 'p10'
      dict['recip_rank'] <- 'rr'
      return(dict[[measure]])
    }
    
    m <- matrix(0.0, nrow = 150, ncol = 61)
    
    col <- 1
    for(f in list.files(terabyte2006_path)) {
      f_name <- paste0(terabyte2006_path, f, '/', f)
      df <- tryCatch(
        read.csv(file=f_name, sep='', skip=1452, header=FALSE),
        error=function(e) e
      )
      if(inherits(df, "error")) next
      df <- filter(df, V1 == measure)
      m[,col] <- df[,'V3']
      col <- col + 1
    }
    
    # Remove last row
    m<-m[1:149,]
    
    # change RR values, to be the closest in the set support('rr')... 
    # which is {0/1000, 1/1000, 2/1000, ..., 1000/1000}
    if(measure=='recip_rank') {
      m2 <- sapply(matchTol(m, support('rr'), tol = 0.001), function(i) {support('rr')[i]})
      m2 <- matrix(unlist(m2), ncol = 61, byrow = FALSE)
      m2 <- round(m2,4)
      
      print(paste0('RR modified values: ', sum(m-m2 != 0)))
      m <- m2
    }
    
    # Add column names "run1", "run2", ...
    col_names <- rep('', 61)
    for(i in 1:61) {
      col_names[i] <- paste0('run', i)
    }
    colnames(m) <- col_names
    
    # Export
    export(m, paste0(path_out, 'terabyte2006_', measure_change_name(measure), '.csv'), 'csv')
  }
}


# ================= Compute area/volume between curves/surfaces ================= 
# Compute area between 2 CDF 2d-curves of marginal distributions
# set only_mean = TRUE to get more values, which are required for the visualization of the deltas
cdf_cdf_diff <- function(cdf_1, cdf_2, measure, X_points = 1001, only_mean=TRUE) {
  if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
    X <- seq(0, 1, length.out = X_points)
  } else { # discrete case
    X <- support(measure)
  }
  
  if(only_mean) {
    return(mean(abs(cdf_1(X) - cdf_2(X))))
  } else {
    cdf_1_Y <- cdf_1(X)
    cdf_2_Y <- cdf_2(X)
    
    # Get highest Y values
    ind_max = which(cdf_1_Y < cdf_2_Y)
    highest_Y = cdf_1_Y
    highest_Y[ind_max] <- cdf_2_Y[ind_max] 
    
    # Get lowest Y values
    ind_min = which(cdf_1_Y > cdf_2_Y)
    lowest_Y = cdf_1_Y
    lowest_Y[ind_min] <- cdf_2_Y[ind_min]
    
    # Get mean absolute diffs
    mean_abs_diff <- mean(abs(cdf_1_Y - cdf_2_Y))
    
    ret <- NULL
    ret$mean <- mean_abs_diff
    ret$X <- X
    ret$lowest_Y <- lowest_Y
    ret$highest_Y <- highest_Y
    ret$cdf_1_Y <- cdf_1_Y
    ret$cdf_2_Y <- cdf_2_Y
    
    return(ret)
  }
}

# Compute volume between 2 joint CDFs. Each CDF is a 3d-surface (of a copula model).
joint_cdf_cdf_diff <- function(cdf_1, cdf_2, n = 101) {
  xy_range <- seq(0, 1, length.out = n)
  
  # n*n Z points
  Z1 <- sapply(xy_range, function(x) cdf_1(xy_range, rep(x, n)))
  Z2 <- sapply(xy_range, function(x) cdf_2(xy_range, rep(x, n)))
  
  return(mean(abs(Z1-Z2)))
}

# ================= Extensions of package VineCopula, simIReff ================= 

# Fixes an arithmetic bug in package VineCopula where CDF(x) can give NA instead of 0
BiCopCDF_2 <- function (x, y, cop) {
  ret <- BiCopCDF(x, y, cop)
  ret <- sapply(ret, 
                function(x) {
                  if (is.na(x)){
                    x <- 0
                  } 
                  return(x)
                })
  return(ret)
}

# Fixes an arithmetic bug in package simIReff where CDF(x) can give values 1
# where "1 > 1" resolves TRUE
peff_2 <- function (q, eff) {
  ret <- peff(q, eff)
  ret <- sapply(ret,
                function(x) {
                  if (x > 1){
                    x <- 1
                  } else if (x < 0) {
                    x <- 0
                  }
                  return(x)
                })
  return(ret)
}

# Fit all marginal models, sort them based on sort_by_criterion
fit_all_margins <- function(x, measure, silent = TRUE, sort_by_criterion = "AIC",
                            compute_splithalf_criterion = FALSE, N_TRIALS = 5) {
  methods <- c("AIC", "BIC", "logLik")
  
  # Fit all margins, continuous
  if(measure %in% c("ap", "ndcg20", "err20")) {
    effs <- as.list(rep(NA, 4))
    try(effs[[1]] <- effCont_norm(x), silent = silent)
    try(effs[[2]] <- effCont_beta(x), silent = silent)
    try(effs[[3]] <- effCont_nks(x), silent = silent)
    try(effs[[4]] <- effCont_bks(x), silent = silent)
  }# Fit all margins, Discrete
  else{ # "p10", "rr"
    effs <- as.list(rep(NA, 5))
    try(effs[[1]] <- effDisc_bbinom(x, support(measure)), silent = silent)
    try(effs[[2]] <- effDisc_dks(x, support(measure)), silent = silent)
    try(effs[[3]] <- effDisc_dks(x, support(measure), mult = 2), silent = silent)
    try(effs[[4]] <- effDisc_dks(x, support(measure), mult = 5), silent = silent)
    try(effs[[5]] <- effDisc_dks(x, support(measure), mult = 10), silent = silent)
  }
  effs <- effs[!is.na(effs)]
  if (length(effs) == 0)
    stop("Unable to fit a marginal distribution")
  
  # Compute some selection criterion values
  scores <- vector(mode = "list")
  for (method in methods){
    FUN <- get(method)
    scores[[method]] <- sapply(effs, FUN)
  }
  
  # Compute my selection criterion
  if(compute_splithalf_criterion) {
    methods <- c(methods, 'splithalf_criterion')
    model_names <- sapply(effs, function(eff) eff$model$type)
    scores[['splithalf_criterion']] <- 
      compute_splithalf_criterion_margins(x, measure, model_names, N_TRIALS)
  }
  
  # Get sorted indices
  scores_sorted_indices <- order(scores[[sort_by_criterion]], 
                                 decreasing = sort_by_criterion == "logLik")
  
  # Add all marginals to a list for export, sorted based on selected criterion
  margins <- vector("list", length(effs))
  for (j in 1:length(margins)) {
    i <- scores_sorted_indices[[j]]
    
    margins[[j]] <- effs[[i]]
    for(method in methods) {
      margins[[j]][[method]] <- scores[[method]][[i]] # save the corresponding criterion values to the margin object
    }
  }
  
  # ----------------------------------
  # Fix a very rare bug where some models that are fitted produces NaNs with their CDF
  # A way to determine such wrongly fitted models seems to be an infinite value on logLik
  # This only seems to happen only when the data that are used to fit are small, i.e. 12 query scores
  margins <- Filter(function(margin) is.finite(margin[['logLik']]), 
                    margins)
  if (length(margins) == 0)
    stop("Unable to fit a marginal distribution")
  # ----------------------------------
  
  return(margins)
}

# Fit all copula models, sort them based on sort_by_criterion
fit_all_copulas <- function(pseudoscores_s1_T1, pseudoscores_s2_T1, sort_by_criterion = "AIC",
                            compute_splithalf_criterion = FALSE, N_TRIALS = 5,
                            scores_s1_T1 = NULL, scores_s2_T1 = NULL, silent=TRUE) {
  
  copulas <- BiCopEstList(pseudoscores_s1_T1, pseudoscores_s2_T1)$models
  
  # Compute my selection criterion
  if(compute_splithalf_criterion) {
    model_names <- sapply(copulas, function(cop) cop$familyname)
    splithalf_criterion_values <- compute_splithalf_criterion_copulas(scores_s1_T1, scores_s2_T1, model_names, N_TRIALS)
    
    for(i in 1:length(copulas)) {
      copulas[[i]][['splithalf_criterion']] <- splithalf_criterion_values[[i]]
    }
  }
  
  sorted_indices <- order(sapply(copulas, function(copula) copula[[sort_by_criterion]]), 
                          decreasing = sort_by_criterion == "logLik")
  # Sort copulas
  copulas <- lapply(sorted_indices, function(x) copulas[[x]])
  return(copulas)
}

# Compute my selection criterion value
compute_splithalf_criterion_copulas <- function (x1, x2, model_names, N_TRIALS) {
  # Columns: the model_names
  # Rows: the trials
  # Initialize with N/A
  df <- data.frame(matrix(ncol = length(model_names), nrow = N_TRIALS))
  colnames(df) <- model_names
  
  # Do a few trial splits for better accuracy
  for(t in 1:N_TRIALS) {
    # Slit-half
    T1_data_indices = sample(seq_len(length(x1)), size = floor(0.5*length(x1)))
    # Paired scores of s1,s2 on T1
    scores_s1_T1 = x1[T1_data_indices]
    scores_s2_T1 = x2[T1_data_indices]
    # Paired scores of s1,s2 on T2
    scores_s1_T2 = x1[-T1_data_indices]
    scores_s2_T2 = x2[-T1_data_indices]
    
    # Get pseudoscores, using ECDF (modified to deal with ties) instead of CDF
    # because we want to evaluate the goodness of fit of the copula models
    # in isolation, meaning without using the model for the marginals
    u1 <- pobs(cbind(scores_s1_T1, scores_s2_T1), ties.method = "random") # compute pseudo-observations
    u2 <- pobs(cbind(scores_s1_T2, scores_s2_T2), ties.method = "random") # compute pseudo-observations
    pseudoscores_s1_T1 <- u1[,1]
    pseudoscores_s2_T1 <- u1[,2]
    pseudoscores_s1_T2 <- u2[,1]
    pseudoscores_s2_T2 <- u2[,2]
    
    # Fit all copulas (on T1 data)
    copulas <- fit_all_copulas(pseudoscores_s1_T1, pseudoscores_s2_T1)
    
    # Select only relevant copulas who appear in model_names
    copulas <- Filter(function(cop) cop$familyname %in% model_names, 
                      copulas)
    
    # Compute CDFs
    cdfs <- sapply(copulas, 
                   # For every copula
                   function(cop) {
                     # Get its CDF
                     cdf <- function (X, Y) { return(BiCopCDF_2(X, Y, cop)) }
                     return(cdf)
                   })
    
    # Define ECDFs
    ecdf_t1 <- function (X, Y) {return(sapply(1:length(X), function(i) mean(pseudoscores_s1_T1 <= X[i] & pseudoscores_s2_T1 <= Y[i])))}
    ecdf_t2 <- function (X, Y) {return(sapply(1:length(X), function(i) mean(pseudoscores_s1_T2 <= X[i] & pseudoscores_s2_T2 <= Y[i])))}
    
    # Compute Deltas
    delta_observed <- sapply(cdfs, function(cdf) joint_cdf_cdf_diff(cdf, ecdf_t2)) 
    delta_expected <- joint_cdf_cdf_diff(ecdf_t1, ecdf_t2)
    
    # Compute my criterion
    for(i in 1:length(copulas)) {
      model_name <- copulas[[i]]$familyname
      # normalized_delta <- (delta_observed[[i]] - delta_expected) / delta_expected
      
      if (model_name %in% model_names) {
        df[t, model_name] <- delta_observed[[i]]
      }
    }
  }
  
  my_criterion_values <- as.list(rep(NA, length(model_names)))
  for(i in 1:length(model_names)){
    model_name <- model_names[[i]]
    val <- mean(df[[model_name]], na.rm = TRUE)
    if(!is.nan(val)){
      my_criterion_values[[i]] <- val
    }
  }
  my_criterion_values <- unlist(my_criterion_values)
  
  return(my_criterion_values)
}

# Compute my selection criterion value
compute_splithalf_criterion_margins <- function (d, measure, model_names, N_TRIALS) {
  # Columns: the model_names
  # Rows: the trials
  # Initialize with N/A
  df <- data.frame(matrix(ncol = length(model_names), nrow = N_TRIALS))
  colnames(df) <- model_names
  
  # Do a few trial splits for better accuracy
  for(t in 1:N_TRIALS) {
    # Split-half
    T1_data_indices <- sample(seq_len(length(d)), size = floor(0.5*length(d)))
    
    T1_data <- d[T1_data_indices]
    T2_data <- d[-T1_data_indices]
    
    # Fit all margins
    margins <- fit_all_margins(T1_data, measure)
    
    # Select only relevant margins who appear in model_names
    margins <- Filter(function(margin) margin$model$type %in% model_names, 
                      margins)
    
    # Compute CDFs
    cdfs <- sapply(margins, 
                   # For every margin
                   function(margin) {
                     # Get its CDF
                     cdf <- function (X) { return(peff_2(X, margin))}
                     return(cdf)
                   })
    # Define ECDFs
    ecdf_t1 <- function (X) { return(sapply(X, function(x) sum(T1_data <= x) / length(T1_data))) }
    ecdf_t2 <- function (X) { return(sapply(X, function(x) sum(T2_data <= x) / length(T2_data))) }
    
    # Compute Deltas
    delta_observed <- sapply(cdfs, function(cdf) cdf_cdf_diff(cdf, ecdf_t2, measure)) 
    delta_expected <- cdf_cdf_diff(ecdf_t1, ecdf_t2, measure)
    
    # Compute my criterion
    for(i in 1:length(margins)) {
      model_name <- margins[[i]]$model$type
      # normalized_delta <- (delta_observed[[i]] - delta_expected) / delta_expected
      
      if (model_name %in% model_names) {
        df[t, model_name] <- delta_observed[[i]]
      }
    }
  }
  
  my_criterion_values <- as.list(rep(NA, length(model_names)))
  
  for(i in 1:length(model_names)){
    model_name <- model_names[[i]]
    val <- mean(df[[model_name]], na.rm = TRUE)
    # print(paste(sd(df[[model_name]], na.rm = TRUE), val))
    if(!is.nan(val)){
      my_criterion_values[[i]] <- val
    }
  }
  my_criterion_values <- unlist(my_criterion_values)
  
  return(my_criterion_values)
}


# ===================== Various Utilities ===================== 

# Converts a list of strings, to a single long string
list_to_str <- function(lst, delim) {
  return(paste(unlist(lst), collapse=delim))
}

# Converts a list of objects to a dataframe
list_to_df <- function(lst) {
  lst <- unlist(lst, recursive=FALSE, use.names=FALSE)
  df <- ldply (lst, data.frame)
  return(df)
}

# Converts a DataFrame's column to a different format. Example:
# Column df$model_names:  {row1: "norm|beta|nks|bks",           row2: "nks|bks|beta",         row3: ... }
#      converts to:       {row1: c("norm","beta","nks","bks"),  row2: c("nks","bks","beta"),  row3: ... }
# If `i` is set, only the i-th element of each list is returned
dfcolumn_str_to_strlist <- function(df_col, delim) {
  lapply(df_col, function(x) str_to_strlist(x, delim))
}

# Converts a DataFrame's column to a different format. Example:
# Column df$T1_data_indices:  {row1: "24|52|66|82",   row2: "3|5|16",   row3: ... }
#      converts to:           {row1: c(24,52,66,82),  row2: c(3,6,16),  row3: ... }
# If `i` is set, only the i-th element of each list is returned
dfcolumn_str_to_numlist <- function(df_col, delim) {
  lapply(df_col, function(x) str_to_numlist(x, delim))
}

# Converts a string (i.e. "word1|word2") to a list("word1", "word2")
str_to_strlist <- function(str, delim) {
  if (delim=='|') {
    delim <- '\\|'
  }
  return(str_split(str, pattern=delim)[[1]])
}

# Converts a string (i.e. "23|42") to a list(23, 42)
str_to_numlist <- function(str, delim) {
  if (delim=='|') {
    delim <- '\\|'
  }
  return(as.numeric(str_split(str, pattern=delim)[[1]]))
}

get_i <- function(lst, i) {
  if (i <= length(lst)) {
    return(lst[[i]])
  } else {
    return(NA)
  }
}

# Add text strip to the right of plot, so that rows have a large label 
add_text_strip_to_the_right <- function(tg, text) {
  posR <- subset(tg$layout, grepl("strip-r", name), select = t:r)
  width <- tg$widths[max(posR$r)]    # width of current right strips
  tg <- gtable_add_cols(tg, width, max(posR$r))
  stripR <- gTree(name = "Strip_right", children = gList(
    rectGrob(gp = gpar(lwd = 2, col = "black", fill = "grey90")),
    textGrob(text, rot = -90, gp = gpar(fontsize = 8.8, col = "grey0"))))
  tg <- gtable_add_grob(tg, stripR, t = min(posR$t), l = max(posR$r) + 1, b = max(posR$b), name = "strip-right")
  tg <- gtable_add_cols(tg, unit(1/5, "line"), max(posR$r))
  return(tg)
}

# Change a plot's row height, for label-rows being cut off
change_plot_row_height <- function(tg, row, cm) {
  tg$heights[tg$layout$b[grep(paste0('strip-r-', row), tg$layout$name)]] = unit(cm, 'cm')
  return(tg)
}

k_signf_digits <- function(x, k) {
  return(format(round(x, k), nsmall = k))
}

sample_my_results <- function(df, percentage) {
  if(percentage < 100) {
    df <- df %>% group_by(measure) %>% sample_frac(percentage / 100.0)
  }
  return(df)
}

percentage_to_str <- function(percentage) {
  if(percentage<100){
    return(paste0('_[perc=', percentage, ']'))
  }
  return('')
}

uniquify_file_name <- function(path_name) {
  if(!file.exists(path_name)) {
    return(path_name)
  }
  
  ext <- paste0('.', file_ext(path_name)) # get file name extension
  file_name <- path_name %>% str_remove(ext) # remove file name extension
  c <- 2
  while(file.exists(path_name)) {
    path_name <- paste0(file_name, '_', c, ext) # add the counter to file name
    c <- c + 1 
  }
  return(path_name)
}

get_all_similar_file_names <- function(path_name) {
  path_names <- vector(mode = "list")
  
  if(!file.exists(path_name)) {
    return(path_names)
  }
  
  path_names[[1]] <- path_name
  ext <- paste0('.', file_ext(path_name)) # find the extension
  file_name <- path_name %>% str_remove(ext) # remove extension
  c <- 2
  while(file.exists(path_name)) {
    path_name <- paste0(file_name, '_', c, ext) # concatenate the counter to the path name
    path_names[[c]] <- path_name
    c <- c + 1 
  }
  path_names <- path_names[-length(path_names)] # remove last path_name
  return(unlist(path_names))
}

# Combine all similar .csv files into a single _COMBINED.csv file
# path_name = 'output/copulas/results_extrapolate_trials=500.csv'
# group_by='measure'
combine_results_df <- function(path_name, group_by='measure') {
  ext <- paste0('.', file_ext(path_name)) # find the extension
  file_name <- path_name %>% str_remove(ext) %>% str_remove(ext) # remove extension
  n_trials_per_measure <- as.numeric(regmatches(path_name, gregexpr("[[:digit:]]+", path_name))[[1]])
  file_name <- file_name %>% str_remove(as.character(n_trials_per_measure)) # remove number of trials
  file_name <- file_name %>% str_remove('trials=')
  file_name
  
  
  path_names <- get_all_similar_file_names(path_name)
  path_names
  dfs <- lapply(path_names, function(x) import(x))
  df <- bind_rows(dfs)
  df <- df[order(df[[group_by]]),]
  n_rows <- nrow(df)
  
  combined_path_name <- paste0(file_name, 'n', n_rows, ext)
  combined_path_name
  if(file.exists(combined_path_name)) {
    stop('File: ', combined_path_name, ' already exists.')
  }
  
  rio::export(df, combined_path_name)
}
