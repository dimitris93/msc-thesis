source("R/utils.R")
stopImplicitCluster()
registerDoParallel(cores = parallel::detectCores())


# ============================= Compute results =============================
compute_results_df <- function(collections, measures, n_trials_per_measure, split_n, output_path,
                               sort_by_criterion = "AIC", compute_splithalf_criterion = FALSE, N_TRIALS = 5) {
  dir.create(output_path, recursive = TRUE)
  
  # Read evaluation data
  dat <- read_evaluation_data(measures, collections, 0.1)
  
  results <- vector(mode = "list")
  
  for(measure in measures) {
    # Count the number of runs per collection, for sampling
    dat_c <- sapply(dat[[measure]], length)
    
    # Run trials
    results[[measure]] <- foreach(trial = 1:n_trials_per_measure) %dopar% {
      source("R/utils.R")
      
      if(split_n[1] > split_n[2]) {
        stop('n2 needs to be larger than n1')
      }
      
      # Sample collection, proportional to number of runs in it
      collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
      # Sample run
      run <- sample(names(dat[[measure]][[collection]]), 1)
      d <- dat[[measure]][[collection]][[run]]
      
      # Split n1-n, n2-n (i.e., 25-99 and 50-99)
      T2_data_indices <- sample(seq_len(length(d)), size = split_n[2])
      T1_data_indices <- sample(T2_data_indices, size = split_n[1])
      T2_data <- d[T2_data_indices]
      T_data <- d[-T2_data_indices]
      T1_data <- d[T1_data_indices]
      
      # Fit all margins
      margins_t1 <- fit_all_margins(T1_data, 
                                    measure,
                                    sort_by_criterion = sort_by_criterion,
                                    compute_splithalf_criterion = compute_splithalf_criterion,
                                    N_TRIALS = N_TRIALS)
      margins_t2 <- fit_all_margins(T2_data, 
                                    measure,
                                    sort_by_criterion = sort_by_criterion,
                                    compute_splithalf_criterion = compute_splithalf_criterion,
                                    N_TRIALS = N_TRIALS)
      
      # Compute CDFs
      cdfs_t1 <- sapply(margins_t1,
                        # For every margin
                        function(margin) {
                          # Get its CDF
                          cdf <- function (X) { return(peff_2(X, margin))}
                          return(cdf)
                        })
      cdfs_t2 <- sapply(margins_t2,
                        # For every margin
                        function(margin) {
                          # Get its CDF
                          cdf <- function (X) { return(peff_2(X, margin))}
                          return(cdf)
                        })
      
      # Define ECDFs
      ecdf_t1 <- function (X) { return(sapply(X, function(x) sum(T1_data <= x) / length(T1_data))) }
      ecdf_t2 <- function (X) { return(sapply(X, function(x) sum(T2_data <= x) / length(T2_data))) }
      ecdf_t <- function (X) { return(sapply(X, function(x) sum(T_data <= x) / length(T_data))) }
      
      # Compute Deltas
      delta_observed_t1 <- sapply(cdfs_t1, function(cdf) cdf_cdf_diff(cdf, ecdf_t, measure))
      delta_observed_t2 <- sapply(cdfs_t2, function(cdf) cdf_cdf_diff(cdf, ecdf_t, measure))
      delta_expected_t1 <- cdf_cdf_diff(ecdf_t1, ecdf_t, measure)
      delta_expected_t2 <- cdf_cdf_diff(ecdf_t2, ecdf_t, measure)
      
      model_names_t1 <- sapply(margins_t1, function(x) x$model$type)
      model_names_t2 <- sapply(margins_t2, function(x) x$model$type)
      AIC_t1 <- sapply(margins_t1, function(x) x$AIC)
      AIC_t2 <- sapply(margins_t2, function(x) x$AIC)
      BIC_t1 <- sapply(margins_t1, function(x) x$BIC)
      BIC_t2 <- sapply(margins_t2, function(x) x$BIC)
      logLik_t1 <- sapply(margins_t1, function(x) x$logLik)
      logLik_t2 <- sapply(margins_t2, function(x) x$logLik)
      if(compute_splithalf_criterion) {
        splithalf_criterion_t1 <- sapply(margins_t1, function(x) x$splithalf_criterion)
        splithalf_criterion_t2 <- sapply(margins_t2, function(x) x$splithalf_criterion)
      }
      random_model_order_t1 <- sample(seq(1, length(margins_t1)))
      random_model_order_t2 <- sample(seq(1, length(margins_t2)))
      
      ret <- NULL
      ret$measure <- measure
      ret$collection <- collection
      ret$run <- run
      ret$T1_data_indices <- list_to_str(T1_data_indices, delim='|')
      ret$T2_data_indices <- list_to_str(T2_data_indices, delim='|')
      ret$model_names_t1 <- list_to_str(model_names_t1, delim='|')
      ret$model_names_t2 <- list_to_str(model_names_t2, delim='|')
      ret$AIC_t1 <- list_to_str(AIC_t1, delim='|')
      ret$AIC_t2 <- list_to_str(AIC_t2, delim='|')
      ret$BIC_t1 <- list_to_str(BIC_t1, delim='|')
      ret$BIC_t2 <- list_to_str(BIC_t2, delim='|')
      ret$logLik_t1 <- list_to_str(logLik_t1, delim='|')
      ret$logLik_t2 <- list_to_str(logLik_t2, delim='|')
      if(compute_splithalf_criterion) {
        ret$splithalf_criterion_t1 <- list_to_str(splithalf_criterion_t1, delim='|')
        ret$splithalf_criterion_t2 <- list_to_str(splithalf_criterion_t2, delim='|')
      }
      ret$random_model_order_t1 <- list_to_str(random_model_order_t1, delim='|')
      ret$random_model_order_t2 <- list_to_str(random_model_order_t2, delim='|')
      ret$delta_observed_t1 <- list_to_str(delta_observed_t1, delim='|')
      ret$delta_observed_t2 <- list_to_str(delta_observed_t2, delim='|')
      ret$delta_expected_t1 <- delta_expected_t1
      ret$delta_expected_t2 <- delta_expected_t2
      return(ret)
    }
  }
  df <- list_to_df(results)
  export(df, uniquify_file_name(paste0(output_path, '/results_extrapolate_trials=', n_trials_per_measure, '.csv')))
}


# ============================= Plot results ============================= 
n_total_splits = 150000
output_path = 'output/margins'
split_n = c(25, 50) # plot '25-25 split vs 50-50 split'


# plot_results(n_trials_per_measure, split_n = c(25, 50), output_path = 'output/margins')
df <- import(paste0(output_path, '/results_extrapolate_n', n_total_splits, '.csv'))

# Rename some columns
colnames(df)[colnames(df) == 'splithalf_criterion_t1'] <- 'SHC_t1'
colnames(df)[colnames(df) == 'logLik_t1'] <- 'LL_t1'

colnames(df)[colnames(df) == 'splithalf_criterion_t2'] <- 'SHC_t2'
colnames(df)[colnames(df) == 'logLik_t2'] <- 'LL_t2'

# Read columns of strings, as lists
df$model_names_t1 <- dfcolumn_str_to_strlist(df$model_names_t1, '|')
df$delta_observed_t1 <- dfcolumn_str_to_numlist(df$delta_observed_t1, '|')
df$AIC_t1 <- dfcolumn_str_to_numlist(df$AIC_t1, '|')
df$SHC_t1 <- dfcolumn_str_to_numlist(df$SHC_t1, '|')
df$BIC_t1 <- dfcolumn_str_to_numlist(df$BIC_t1, '|')
df$LL_t1 <- dfcolumn_str_to_numlist(df$LL_t1, '|')
df$random_model_order_t1 <- dfcolumn_str_to_numlist(df$random_model_order_t1, '|')
df$T1_data_indices <- dfcolumn_str_to_numlist(df$T1_data_indices, '|')

df$model_names_t2 <- dfcolumn_str_to_strlist(df$model_names_t2, '|')
df$delta_observed_t2 <- dfcolumn_str_to_numlist(df$delta_observed_t2, '|')
df$AIC_t2 <- dfcolumn_str_to_numlist(df$AIC_t2, '|')
df$SHC_t2 <- dfcolumn_str_to_numlist(df$SHC_t2, '|')
df$BIC_t2 <- dfcolumn_str_to_numlist(df$BIC_t2, '|')
df$LL_t2 <- dfcolumn_str_to_numlist(df$LL_t2, '|')
df$random_model_order_t2 <- dfcolumn_str_to_numlist(df$random_model_order_t2, '|')
df$T2_data_indices <- dfcolumn_str_to_numlist(df$T2_data_indices, '|')

# Calculate T1_data, T2_data, T_data
collections <- c("terabyte2006")
measures <- c("ap", "p10", "rr")
df$collection <- 'terabyte2006'
dat <- read_evaluation_data(measures, collections, 0.1)
df$T1_data <- apply(df, 1, function(x) return(list(dat[[x$measure]][[x$collection]][[x$run]][x$T1_data_indices])))
df$T1_data <- lapply(df$T1_data, function (x) unlist(x, use.names = F))
df$T2_data <- apply(df, 1, function(x) return(list(dat[[x$measure]][[x$collection]][[x$run]][x$T2_data_indices])))
df$T2_data <- lapply(df$T2_data, function (x) unlist(x, use.names = F))
df$T_data <- apply(df, 1, function(x) return(list(dat[[x$measure]][[x$collection]][[x$run]][-x$T2_data_indices])))
df$T_data <- lapply(df$T_data, function (x) unlist(x, use.names = F))
df$m <- df$measure
df$measure <- recode_and_reorder(df$measure) # do this after the calculation of T1_data
df$T1_data_perc_zeros <- sapply(df$T1_data, function(x) length(which(x==0))/length(x))
df$T1_data_count_zeros <- sapply(df$T1_data, function(x) length(which(x==0)))
df$T2_data_perc_zeros <- sapply(df$T2_data, function(x) length(which(x==0))/length(x))
df$T2_data_count_zeros <- sapply(df$T2_data, function(x) length(which(x==0)))
df$T_data_perc_zeros <- sapply(df$T_data, function(x) length(which(x==0))/length(x))
df$T_data_count_zeros <- sapply(df$T_data, function(x) length(which(x==0)))

######################################
################   T1 ################
######################################
# ================= Best AIC ==================
df$best_model_name_AIC_t1 <- recode_and_reorder(sapply(df$model_names_t1, function(x) x[[1]])) # 1st model is the best fit based on AIC
df$best_delta_observed_AIC_t1 <- sapply(df$delta_observed_t1, function(x) x[[1]]) # 1st model is the best fit based on AIC
df$GoF_AIC_t1 <- -(df$best_delta_observed_AIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best AIC, except beta ks ==================
df$best_model_name_AIC_no_bks_bestindex_t1 <- sapply(df$model_names_t1, function(x) if (x[[1]]=='bks') 2 else 1)
df$best_model_name_AIC_no_bks_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[x$best_model_name_AIC_no_bks_bestindex_t1[[1]]]]))
df$best_delta_observed_AIC_no_bks_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[x$best_model_name_AIC_no_bks_bestindex_t1[[1]]]])
df$GoF_AIC_no_bks_t1 <- -(df$best_delta_observed_AIC_no_bks_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= 2nd best AIC ==================
df$top2_model_name_AIC_t1 <- recode_and_reorder(sapply(df$model_names_t1, function(x) get_i(x, 2))) # 2nd model is the top-2 fit based on AIC
df$top2_delta_observed_AIC_t1 <- sapply(df$delta_observed_t1, function(x) get_i(x, 2)) # 2nd model is the top-2 fit based on AIC
df$GoF_top2_AIC_t1 <- -(df$top2_delta_observed_AIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= 3rd best AIC ==================
df$top3_model_name_AIC_t1 <- recode_and_reorder(sapply(df$model_names_t1, function(x) get_i(x, 3))) # 3rd model is the top-3 fit based on AIC
df$top3_delta_observed_AIC_t1 <- sapply(df$delta_observed_t1, function(x) get_i(x, 3)) # 3rd model is the top-3 fit based on AIC
df$GoF_top3_AIC_t1 <- -(df$top3_delta_observed_AIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best SHC ==================
df$best_model_name_SHC_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$SHC_t1)[[1]]]]))
df$best_delta_observed_SHC_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$SHC_t1)[[1]]]])
df$GoF_SHC_t1 <- -(df$best_delta_observed_SHC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best SHC, except beta ks ==================
df$best_model_name_SHC_no_bks_bestindex_t1 <- apply(df, 1, function(x) if (x$model_names_t1[[order(x$SHC_t1)[[1]]]] =='bks') order(x$SHC_t1)[[2]] else order(x$SHC_t1)[[1]])
df$best_model_name_SHC_no_bks_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[x$best_model_name_SHC_no_bks_bestindex_t1[[1]]]]))
df$best_delta_observed_SHC_no_bks_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[x$best_model_name_SHC_no_bks_bestindex_t1[[1]]]])
df$GoF_SHC_no_bks_t1 <- -(df$best_delta_observed_SHC_no_bks_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best BIC ==================
df$best_model_name_BIC_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$BIC_t1)[[1]]]]))
df$best_delta_observed_BIC_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$BIC_t1)[[1]]]])
df$GoF_BIC_t1 <- -(df$best_delta_observed_BIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best BIC, except beta ks ==================
df$best_model_name_BIC_no_bks_bestindex_t1 <- apply(df, 1, function(x) if (x$model_names_t1[[order(x$BIC_t1)[[1]]]] =='bks') order(x$BIC_t1)[[2]] else order(x$BIC_t1)[[1]])
df$best_model_name_BIC_no_bks_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[x$best_model_name_BIC_no_bks_bestindex_t1[[1]]]]))
df$best_delta_observed_BIC_no_bks_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[x$best_model_name_BIC_no_bks_bestindex_t1[[1]]]])
df$GoF_BIC_no_bks_t1 <- -(df$best_delta_observed_BIC_no_bks_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best LL ==================
df$best_model_name_LL_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$LL_t1, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$LL_t1, decreasing = TRUE)[[1]]]])
df$GoF_LL_t1 <- -(df$best_delta_observed_LL_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best LL, except beta ks ==================
df$best_model_name_LL_no_bks_bestindex_t1 <- apply(df, 1, function(x) if (x$model_names_t1[[order(x$LL_t1, decreasing = T)[[1]]]] =='bks') order(x$LL_t1, decreasing = T)[[2]] else order(x$LL_t1, decreasing = T)[[1]])
df$best_model_name_LL_no_bks_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[x$best_model_name_LL_no_bks_bestindex_t1[[1]]]]))
df$best_delta_observed_LL_no_bks_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[x$best_model_name_LL_no_bks_bestindex_t1[[1]]]])
df$GoF_LL_no_bks_t1 <- -(df$best_delta_observed_LL_no_bks_t1 - df$delta_expected_t1) / df$delta_expected_t1

# =============== Random model ===============
df$best_model_name_RAND_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[x$random_model_order_t1[[1]]]]))
df$best_delta_observed_RAND_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[x$random_model_order_t1[[1]]]])
df$GoF_RAND_t1 <- -(df$best_delta_observed_RAND_t1 - df$delta_expected_t1) / df$delta_expected_t1

# =============== Theoretical Best model (according to lowest Delta observed) ===============
df$best_model_name_ACTUAL_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$delta_observed_t1)[[1]]]]))
df$best_delta_observed_ACTUAL_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$delta_observed_t1)[[1]]]])
df$GoF_ACTUAL_t1 <- -(df$best_delta_observed_ACTUAL_t1 - df$delta_expected_t1) / df$delta_expected_t1

# =============== Theoretical worst model (according to highest Delta observed) ===============
df$best_model_name_WORST_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$delta_observed_t1, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_WORST_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$delta_observed_t1, decreasing = TRUE)[[1]]]])
df$GoF_WORST_t1 <- -(df$best_delta_observed_WORST_t1 - df$delta_expected_t1) / df$delta_expected_t1


#####################################
################  T2 ################
#####################################

# ================= Best AIC ==================
df$best_model_name_AIC_t2 <- recode_and_reorder(sapply(df$model_names_t2, function(x) x[[1]])) # 1st model is the best fit based on AIC
df$best_delta_observed_AIC_t2 <- sapply(df$delta_observed_t2, function(x) x[[1]]) # 1st model is the best fit based on AIC
df$GoF_AIC_t2 <- -(df$best_delta_observed_AIC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best AIC, except beta ks ==================
df$best_model_name_AIC_no_bks_bestindex_t2 <- sapply(df$model_names_t2, function(x) if (x[[1]]=='bks') 2 else 1)
df$best_model_name_AIC_no_bks_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[x$best_model_name_AIC_no_bks_bestindex_t2[[1]]]]))
df$best_delta_observed_AIC_no_bks_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[x$best_model_name_AIC_no_bks_bestindex_t2[[1]]]])
df$GoF_AIC_no_bks_t2 <- -(df$best_delta_observed_AIC_no_bks_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= 2nd best AIC ==================
df$top2_model_name_AIC_t2 <- recode_and_reorder(sapply(df$model_names_t2, function(x) get_i(x, 2))) # 2nd model is the top-2 fit based on AIC
df$top2_delta_observed_AIC_t2 <- sapply(df$delta_observed_t2, function(x) get_i(x, 2)) # 2nd model is the top-2 fit based on AIC
df$GoF_top2_AIC_t2 <- -(df$top2_delta_observed_AIC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= 3rd best AIC ==================
df$top3_model_name_AIC_t2 <- recode_and_reorder(sapply(df$model_names_t2, function(x) get_i(x, 3))) # 3rd model is the top-3 fit based on AIC
df$top3_delta_observed_AIC_t2 <- sapply(df$delta_observed_t2, function(x) get_i(x, 3)) # 3rd model is the top-3 fit based on AIC
df$GoF_top3_AIC_t2 <- -(df$top3_delta_observed_AIC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best SHC ==================
df$best_model_name_SHC_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$SHC_t2)[[1]]]]))
df$best_delta_observed_SHC_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$SHC_t2)[[1]]]])
df$GoF_SHC_t2 <- -(df$best_delta_observed_SHC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best SHC, except beta ks ==================
df$best_model_name_SHC_no_bks_bestindex_t2 <- apply(df, 1, function(x) if (x$model_names_t2[[order(x$SHC_t2)[[1]]]] =='bks') order(x$SHC_t2)[[2]] else order(x$SHC_t2)[[1]])
df$best_model_name_SHC_no_bks_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[x$best_model_name_SHC_no_bks_bestindex_t2[[1]]]]))
df$best_delta_observed_SHC_no_bks_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[x$best_model_name_SHC_no_bks_bestindex_t2[[1]]]])
df$GoF_SHC_no_bks_t2 <- -(df$best_delta_observed_SHC_no_bks_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best BIC ==================
df$best_model_name_BIC_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$BIC_t2)[[1]]]]))
df$best_delta_observed_BIC_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$BIC_t2)[[1]]]])
df$GoF_BIC_t2 <- -(df$best_delta_observed_BIC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best BIC, except beta ks ==================
df$best_model_name_BIC_no_bks_bestindex_t2 <- apply(df, 1, function(x) if (x$model_names_t2[[order(x$BIC_t2)[[1]]]] =='bks') order(x$BIC_t2)[[2]] else order(x$BIC_t2)[[1]])
df$best_model_name_BIC_no_bks_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[x$best_model_name_BIC_no_bks_bestindex_t2[[1]]]]))
df$best_delta_observed_BIC_no_bks_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[x$best_model_name_BIC_no_bks_bestindex_t2[[1]]]])
df$GoF_BIC_no_bks_t2 <- -(df$best_delta_observed_BIC_no_bks_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best LL ==================
df$best_model_name_LL_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$LL_t2, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$LL_t2, decreasing = TRUE)[[1]]]])
df$GoF_LL_t2 <- -(df$best_delta_observed_LL_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best LL, except beta ks ==================
df$best_model_name_LL_no_bks_bestindex_t2 <- apply(df, 1, function(x) if (x$model_names_t2[[order(x$LL_t2, decreasing = T)[[1]]]] =='bks') order(x$LL_t2, decreasing = T)[[2]] else order(x$LL_t2, decreasing = T)[[1]])
df$best_model_name_LL_no_bks_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[x$best_model_name_LL_no_bks_bestindex_t2[[1]]]]))
df$best_delta_observed_LL_no_bks_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[x$best_model_name_LL_no_bks_bestindex_t2[[1]]]])
df$GoF_LL_no_bks_t2 <- -(df$best_delta_observed_LL_no_bks_t2 - df$delta_expected_t2) / df$delta_expected_t2

# =============== Random model ===============
df$best_model_name_RAND_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[x$random_model_order_t2[[1]]]]))
df$best_delta_observed_RAND_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[x$random_model_order_t2[[1]]]])
df$GoF_RAND_t2 <- -(df$best_delta_observed_RAND_t2 - df$delta_expected_t2) / df$delta_expected_t2

# =============== Theoretical Best model (according to lowest Delta observed) ===============
df$best_model_name_ACTUAL_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$delta_observed_t2)[[1]]]]))
df$best_delta_observed_ACTUAL_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$delta_observed_t2)[[1]]]])
df$GoF_ACTUAL_t2 <- -(df$best_delta_observed_ACTUAL_t2 - df$delta_expected_t2) / df$delta_expected_t2

# =============== Theoretical worst model (according to highest Delta observed) ===============
df$best_model_name_WORST_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$delta_observed_t2, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_WORST_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$delta_observed_t2, decreasing = TRUE)[[1]]]])
df$GoF_WORST_t2 <- -(df$best_delta_observed_WORST_t2 - df$delta_expected_t2) / df$delta_expected_t2



# ---------- Plot 1 - Dobs 25, Dobs 50, Dexp 25 ----------
for(measure in measures) {
  df2 <- df[df[, "m"] == measure,]
  
  my_colors <- c('black', 'red')
  my_shapes <- list('', 8)
  my_lines <- c('solid', 'solid')
  my_labels <- c(bquote(Delta[obs]^{'25'}), bquote(Delta[obs]^{'50'}))
  
  mean_df <- df2 %>%
    group_by(best_model_name_AIC_t1) %>%
    summarise(mean_val = mean(best_delta_observed_AIC_t1))
  
  mean_df2 <- df2 %>%
    group_by(best_model_name_AIC_t1) %>%
    summarise(mean_val2 = mean(best_delta_observed_AIC_t2))
  
  legend_pos <- if(measure == 'p10') 'bottom' else 'none'
  h <- if(measure == 'p10') 5.3 else 4.732
  pdf(paste0(output_path, '/extrapolate_n150000_fig1_', measure,'.pdf'), width=3, height=h, onefile=F)
  p <- ggplot() +
    ggtitle(beautify(measure)) +
    geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
    
    geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[1]), linetype = my_lines[1]) +
    geom_vline(data=mean_df2, aes(xintercept=mean_val2, color=my_colors[2]), linetype = 'dotted') +
    
    geom_point(data=df2, aes(x=best_delta_observed_AIC_t2, y=best_model_name_AIC_t2, color=my_colors[2]),
               shape=my_shapes[2], stat="summary", fun="mean") +
    geom_linerange(data=df2, aes(x=best_delta_observed_AIC_t2, y=best_model_name_AIC_t2, color=my_colors[2]),
                   stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
    
    scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
    guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
    facet_grid(best_model_name_AIC_t1~., space = "free", scale = "free") +
    labs(y= "Best model (AIC) fitted on 50 topics", x='') +
    scale_x_continuous(limits = c(0, 0.13), breaks =  seq(0, 0.12, 0.04), oob = scales::oob_keep) +
    theme(plot.margin = margin(1, 1, 0, 1, "mm"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(size = 11),
          legend.position = legend_pos)
  # print(p)
  # dev.off()
  
  tg <- ggplotGrob(p)
  tg <- add_text_strip_to_the_right(tg, 'Best model (AIC) fitted on 25 topics')
  if(measure == 'p10') {
    tg <- change_plot_row_height(tg, row = 4, cm = 1.9)
  }
  if(measure == 'rr') {
    tg <- change_plot_row_height(tg, row = 5, cm = 2)
  }
  p <- as.ggplot(tg)
  
  print(p)
  dev.off()
}


# ---------- Table 1 ----------
mean_df <- df %>%
  group_by(measure) %>%
  summarise(dobs50_minus_dobs25 = mean(best_delta_observed_AIC_t2 - best_delta_observed_AIC_t1))
mean_df

mean_df <- df %>%
  group_by(measure) %>%
  summarise(dobs50_minus_dobs25_div_dobs25 = mean((best_delta_observed_AIC_t2 - best_delta_observed_AIC_t1) / best_delta_observed_AIC_t1))
mean_df

mean_df <- df %>%
  group_by(measure) %>%
  summarise(dobs50_minus_dobs25_div_dobs25 = mean((best_delta_observed_AIC_t2 - best_delta_observed_AIC_t1) / delta_expected_t1))
mean_df

mean_df <- filter(df, best_model_name_AIC_t1=='Beta KS') %>%
  group_by(measure) %>%
  summarise(dobs_51 = mean((best_delta_observed_AIC_t2 - best_delta_observed_AIC_t1) / best_delta_observed_AIC_t1))
mean_df


# ---------- Fig 2 -----------
pdf(paste0(output_path, '/extrapolate_n', nrow(df), '_fig2.pdf'), width=3, height=4.65)
my_colors <- c('blue', 'red', 'black')
my_shapes <- list(16, 8, '')
my_lines <- c('solid', 'solid', 'solid')
my_labels <- c(bquote(Delta[obs]^{'25'}), bquote(Delta[obs]^{'50'}), bquote(Delta[exp]^{'25'}))

mean_df <- df %>%
  group_by(measure) %>%
  summarise(mean_val = mean(delta_expected_t1))

p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  
  # T1
  geom_point(data=df, aes(x=best_delta_observed_AIC_t1, y='AIC', color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_AIC_t1, y='AIC', color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_BIC_t1, y='BIC', color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_BIC_t1, y='BIC', color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_LL_t1, y='LL', color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_LL_t1, y='LL', color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_SHC_t1, y='SHC', color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_SHC_t1, y='SHC', color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_ACTUAL_t1, y='Best-case', color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_ACTUAL_t1, y='Best-case', color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  # T2 
  geom_point(data=df, aes(x=best_delta_observed_AIC_t2, y='AIC', color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_AIC_t2, y='AIC', color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_BIC_t2, y='BIC', color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_BIC_t2, y='BIC', color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_LL_t2, y='LL', color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_LL_t2, y='LL', color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_SHC_t2, y='SHC', color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_SHC_t2, y='SHC', color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_ACTUAL_t2, y='Best-case', color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_ACTUAL_t2, y='Best-case', color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  # exp T1
  geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[3]), linetype = 'solid') +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= '', x = '') +
  theme(legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, -0.5, 1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()



# # Run this repeatedly and then combine results to one dataframe
# compute_results_df(split_n = c(25,50),
#                    collections = c("terabyte2006"),
#                    measures =  c("ap", "p10", "rr"),
#                    n_trials_per_measure = 1000,
#                    output_path = 'output/margins',
#                    sort_by_criterion = "AIC",
#                    compute_splithalf_criterion = TRUE,
#                    N_TRIALS = 10)

# # Combine to one dataframe
# combine_results_df('output/margins/results_extrapolate_trials=1000.csv')
