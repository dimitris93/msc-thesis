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
      # Sample 2 runs
      runs <- sample(names(dat[[measure]][[collection]]), 2)
      scores_s1 <- dat[[measure]][[collection]][[runs[[1]]]]
      scores_s2 <- dat[[measure]][[collection]][[runs[[2]]]]
      
      # Split n1-n, n2-n (i.e., 25-99 and 50-99)
      T2_data_indices = sample(seq_len(length(scores_s1)), size = split_n[2]) 
      T1_data_indices <- sample(T2_data_indices, size = split_n[1]) 
      # Paired scores of s1,s2 on T1
      scores_s1_T1 = scores_s1[T1_data_indices]
      scores_s2_T1 = scores_s2[T1_data_indices]
      # Paired scores of s1,s2 on T2
      scores_s1_T2 = scores_s1[T2_data_indices]
      scores_s2_T2 = scores_s2[T2_data_indices]
      # Paired scores of s1,s2 on T
      scores_s1_T = scores_s1[-T2_data_indices]
      scores_s2_T = scores_s2[-T2_data_indices]
      
      # Get pseudoscores, using ECDF (modified to deal with ties) instead of CDF
      # because we want to evaluate the goodness of fit of the copula models
      # in isolation, meaning without using the model for the marginals
      u1 <- pobs(cbind(scores_s1_T1, scores_s2_T1), ties.method = "random") # compute pseudo-observations
      u2 <- pobs(cbind(scores_s1_T2, scores_s2_T2), ties.method = "random") # compute pseudo-observations
      u <- pobs(cbind(scores_s1_T, scores_s2_T), ties.method = "random") # compute pseudo-observations
      
      pseudoscores_s1_T1 <- u1[,1]
      pseudoscores_s2_T1 <- u1[,2]
      pseudoscores_s1_T2 <- u2[,1]
      pseudoscores_s2_T2 <- u2[,2]
      pseudoscores_s1_T <- u[,1]
      pseudoscores_s2_T <- u[,2]
      
      # Fit all copulas
      copulas_t1 <- NULL
      try(copulas_t1 <- fit_all_copulas(pseudoscores_s1_T1, 
                                        pseudoscores_s2_T1, 
                                        sort_by_criterion = sort_by_criterion,
                                        compute_splithalf_criterion = compute_splithalf_criterion, 
                                        N_TRIALS = N_TRIALS,
                                        scores_s1_T1 = scores_s1_T1,
                                        scores_s2_T1 = scores_s2_T1), silent = TRUE)
      
      copulas_t2 <- NULL
      try(copulas_t2 <- fit_all_copulas(pseudoscores_s1_T2, 
                                        pseudoscores_s2_T2, 
                                        sort_by_criterion = sort_by_criterion,
                                        compute_splithalf_criterion = compute_splithalf_criterion, 
                                        N_TRIALS = N_TRIALS,
                                        scores_s1_T1 = scores_s1_T2,
                                        scores_s2_T1 = scores_s2_T2), silent = TRUE)
      
      
      if(is.null(copulas_t1)) {
        return(NULL)
      }
      if(is.null(copulas_t2)) {
        return(NULL)
      }
      
      # Compute CDFs
      cdfs_t1 <- sapply(copulas_t1, 
                        # For every copula
                        function(cop) {
                          # Get its CDF
                          cdf <- function (X, Y) { return(BiCopCDF_2(X, Y, cop)) }
                          return(cdf)
                        })
      cdfs_t2 <- sapply(copulas_t2, 
                        # For every copula
                        function(cop) {
                          # Get its CDF
                          cdf <- function (X, Y) { return(BiCopCDF_2(X, Y, cop)) }
                          return(cdf)
                        })
      
      # Define ECDFs
      ecdf_t1 <- function (X, Y) {return(sapply(1:length(X), function(i) mean(pseudoscores_s1_T1 <= X[i] & pseudoscores_s2_T1 <= Y[i])))}
      ecdf_t2 <- function (X, Y) {return(sapply(1:length(X), function(i) mean(pseudoscores_s1_T2 <= X[i] & pseudoscores_s2_T2 <= Y[i])))}
      ecdf_t <- function (X, Y) {return(sapply(1:length(X), function(i) mean(pseudoscores_s1_T <= X[i] & pseudoscores_s2_T <= Y[i])))}
      
      # Compute Deltas
      delta_observed_t1 <- sapply(cdfs_t1, function(cdf) joint_cdf_cdf_diff(cdf, ecdf_t)) 
      delta_observed_t2 <- sapply(cdfs_t2, function(cdf) joint_cdf_cdf_diff(cdf, ecdf_t)) 
      delta_expected_t1 <- joint_cdf_cdf_diff(ecdf_t1, ecdf_t)
      delta_expected_t2 <- joint_cdf_cdf_diff(ecdf_t2, ecdf_t)
      
      model_names_t1 <- sapply(copulas_t1, function(x) x$familyname)
      model_names_t2 <- sapply(copulas_t2, function(x) x$familyname)
      AIC_t1 <- sapply(copulas_t1, function(x) x$AIC)
      AIC_t2 <- sapply(copulas_t2, function(x) x$AIC)
      BIC_t1 <- sapply(copulas_t1, function(x) x$BIC)
      BIC_t2 <- sapply(copulas_t2, function(x) x$BIC)
      logLik_t1 <- sapply(copulas_t1, function(x) x$logLik)
      logLik_t2 <- sapply(copulas_t2, function(x) x$logLik)
      if(compute_splithalf_criterion) {
        splithalf_criterion_t1 <- sapply(copulas_t1, function(x) x$splithalf_criterion)
        splithalf_criterion_t2 <- sapply(copulas_t2, function(x) x$splithalf_criterion)
      }
      random_model_order_t1 <- sample(seq(1, length(copulas_t1)))
      random_model_order_t2 <- sample(seq(1, length(copulas_t2)))
      
      ret <- NULL
      ret$measure <- measure
      ret$collection <- collection
      ret$run1 <- runs[[1]]
      ret$run2 <- runs[[2]]
      ret$T1_data_indices <- list_to_str(T1_data_indices, delim='|')
      ret$T2_data_indices <- list_to_str(T2_data_indices, delim='|')
      ret$pseudoscores_s1_T1 <- list_to_str(pseudoscores_s1_T1, delim='|')
      ret$pseudoscores_s2_T1 <- list_to_str(pseudoscores_s2_T1, delim='|')
      ret$pseudoscores_s1_T2 <- list_to_str(pseudoscores_s1_T2, delim='|')
      ret$pseudoscores_s2_T2 <- list_to_str(pseudoscores_s2_T2, delim='|')
      ret$pseudoscores_s1_T <- list_to_str(pseudoscores_s1_T, delim='|')
      ret$pseudoscores_s2_T <- list_to_str(pseudoscores_s2_T, delim='|')
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
  write.csv(df, uniquify_file_name(paste0(output_path, '/results_extrapolate_trials=', n_trials_per_measure, '.csv')))
}


# # Run this repeatedly and then combine results to one dataframe
# compute_results_df(split_n = c(25,50),
#                    collections = c("terabyte2006"),
#                    measures =  c("ap", "p10", "rr"),
#                    n_trials_per_measure = 500,
#                    output_path = 'output/copulas',
#                    sort_by_criterion = "AIC",
#                    compute_splithalf_criterion = TRUE,
#                    N_TRIALS = 10)


# # Combine to one dataframe
# combine_results_df('output/copulas/results_extrapolate_trials=500.csv')



# ============================= Plot results =============================
n_total_splits = 150000
output_path = 'output/copulas'
split_n = c(25, 50) # plot '25-25 split vs 50-50 split'
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
df$m <- df$measure
df$measure <- recode_and_reorder(df$measure) # do this after the calculation of T1_data

######################################
################   T1 ################
######################################
# ================= Best AIC ==================
df$best_model_name_AIC_t1 <- recode_and_reorder(sapply(df$model_names_t1, function(x) x[[1]])) # 1st model is the best fit based on AIC
df$best_delta_observed_AIC_t1 <- sapply(df$delta_observed_t1, function(x) x[[1]]) # 1st model is the best fit based on AIC
df$GoF_AIC_t1 <- -(df$best_delta_observed_AIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best SHC ==================
df$best_model_name_SHC_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$SHC_t1)[[1]]]]))
df$best_delta_observed_SHC_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$SHC_t1)[[1]]]])
df$GoF_SHC_t1 <- -(df$best_delta_observed_SHC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best BIC ==================
df$best_model_name_BIC_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$BIC_t1)[[1]]]]))
df$best_delta_observed_BIC_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$BIC_t1)[[1]]]])
df$GoF_BIC_t1 <- -(df$best_delta_observed_BIC_t1 - df$delta_expected_t1) / df$delta_expected_t1

# ================= Best LL ==================
df$best_model_name_LL_t1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t1[[order(x$LL_t1, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL_t1 <- apply(df, 1, function(x) x$delta_observed_t1[[order(x$LL_t1, decreasing = TRUE)[[1]]]])
df$GoF_LL_t1 <- -(df$best_delta_observed_LL_t1 - df$delta_expected_t1) / df$delta_expected_t1

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

# ================= Best SHC ==================
df$best_model_name_SHC_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$SHC_t2)[[1]]]]))
df$best_delta_observed_SHC_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$SHC_t2)[[1]]]])
df$GoF_SHC_t2 <- -(df$best_delta_observed_SHC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best BIC ==================
df$best_model_name_BIC_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$BIC_t2)[[1]]]]))
df$best_delta_observed_BIC_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$BIC_t2)[[1]]]])
df$GoF_BIC_t2 <- -(df$best_delta_observed_BIC_t2 - df$delta_expected_t2) / df$delta_expected_t2

# ================= Best LL ==================
df$best_model_name_LL_t2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names_t2[[order(x$LL_t2, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL_t2 <- apply(df, 1, function(x) x$delta_observed_t2[[order(x$LL_t2, decreasing = TRUE)[[1]]]])
df$GoF_LL_t2 <- -(df$best_delta_observed_LL_t2 - df$delta_expected_t2) / df$delta_expected_t2

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



# ---------- Plot 1 - Dobs 25, Dobs 50 ----------
my_colors <- c('black', 'red')
my_shapes <- list(16, 8)
my_lines <- c('solid', 'solid')
my_labels <- c(bquote(Delta[obs]^{'25'}), bquote(Delta[obs]^{'50'}))

mean_df <- df %>%
  group_by(measure) %>%
  summarise(mean_val = mean(best_delta_observed_AIC_t1))

mean_df2 <- df %>%
  group_by(measure) %>%
  summarise(mean_val2 = mean(best_delta_observed_AIC_t2))

pdf(paste0(output_path, '/extrapolate_n150000_fig1.pdf'), width=3, height=5.3)
p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[1]), linetype = 'dotted') +
  geom_vline(data=mean_df2, aes(xintercept=mean_val2, color=my_colors[2]), linetype = 'dotted') +
  
  geom_point(data=df, aes(x=best_delta_observed_AIC_t1, y=best_model_name_AIC_t1, color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_AIC_t1, y=best_model_name_AIC_t1, color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_AIC_t2, y=best_model_name_AIC_t1, color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_AIC_t2, y=best_model_name_AIC_t1, color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= "Best model (AIC) fitted on 25 topics", x='') +
  # scale_x_continuous(limits = c(0, 0.13), breaks =  seq(0, 0.12, 0.04), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 1, -1, 1, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position = 'bottom')
print(p)
dev.off()


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
