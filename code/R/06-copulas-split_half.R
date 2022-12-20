source("R/utils.R")
stopImplicitCluster()
registerDoParallel(cores = parallel::detectCores())


# ============================= Compute results =============================
compute_results_df <- function(collections, measures, n_trials_per_measure, output_path, sort_by_criterion = "AIC",
                               compute_splithalf_criterion = FALSE, N_TRIALS = 5) {
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
      
      # Sample collection, proportional to number of runs in it
      collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
      # Sample 2 runs
      runs <- sample(names(dat[[measure]][[collection]]), 2)
      scores_s1 <- dat[[measure]][[collection]][[runs[[1]]]]
      scores_s2 <- dat[[measure]][[collection]][[runs[[2]]]]
      
      # Slit-half
      T1_data_indices = sample(seq_len(length(scores_s1)), size = floor(0.5*length(scores_s1)))
      # Paired scores of s1,s2 on T1
      scores_s1_T1 = scores_s1[T1_data_indices]
      scores_s2_T1 = scores_s2[T1_data_indices]
      # Paired scores of s1,s2 on T2
      scores_s1_T2 = scores_s1[-T1_data_indices]
      scores_s2_T2 = scores_s2[-T1_data_indices]
      
      # Get pseudoscores, using ECDF (modified to deal with ties) instead of CDF
      # because we want to evaluate the goodness of fit of the copula models
      # in isolation, meaning without using the model for the marginals
      u <- pobs(cbind(scores_s1_T1, scores_s2_T1, scores_s1_T2, scores_s2_T2), ties.method = "random") # compute pseudo-observations
      pseudoscores_s1_T1 <- u[,1]
      pseudoscores_s2_T1 <- u[,2]
      pseudoscores_s1_T2 <- u[,3]
      pseudoscores_s2_T2 <- u[,4]
      
      # Fit all copulas (on T1 data)
      # I am doing a try catch because of a rare bug situation with BiCopEstList
      # bug 1: In BiCop: The second parameter of the BB8 copula has to be in the interval [1e-4,1].
      # bug 2: In BiCop: The second parameter of the rotated BB8 copula has to be in the interval [-1,-1e-4].
      # bug 3: In BiCop: The parameter of the Frank copula has to be unequal to 0.
      # It looks like maximum likelihood optimization gives parameter values that are not good?
      # it is a pretty rare bug so I just ignore those cases and return NULL here. 
      copulas <- NULL
      try(copulas <- fit_all_copulas(pseudoscores_s1_T1, 
                                     pseudoscores_s2_T1, 
                                     sort_by_criterion = sort_by_criterion,
                                     compute_splithalf_criterion = compute_splithalf_criterion, 
                                     N_TRIALS = N_TRIALS,
                                     scores_s1_T1 = scores_s1_T1,
                                     scores_s2_T1 = scores_s2_T1), 
          silent = TRUE)
      if(is.null(copulas)) {
        return(NULL)
      }
      
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
      
      model_names <- sapply(copulas, function(x) x$familyname)
      AIC <- sapply(copulas, function(x) x$AIC)
      BIC <- sapply(copulas, function(x) x$BIC)
      logLik <- sapply(copulas, function(x) x$logLik)
      if(compute_splithalf_criterion) {
        splithalf_criterion <- sapply(copulas, function(x) x$splithalf_criterion)
      }
      random_model_order <- sample(seq(1, length(copulas)))
      
      ret <- NULL
      ret$measure <- measure
      ret$collection <- collection
      ret$run1 <- runs[[1]]
      ret$run2 <- runs[[2]]
      ret$T1_data_indices <- list_to_str(T1_data_indices, delim='|')
      ret$pseudoscores_s1_T1 <- list_to_str(pseudoscores_s1_T1, delim='|')
      ret$pseudoscores_s2_T1 <- list_to_str(pseudoscores_s2_T1, delim='|')
      ret$pseudoscores_s1_T2 <- list_to_str(pseudoscores_s1_T2, delim='|')
      ret$pseudoscores_s2_T2 <- list_to_str(pseudoscores_s2_T2, delim='|')
      ret$model_names <- list_to_str(model_names, delim='|')
      ret$AIC <- list_to_str(AIC, delim='|')
      ret$BIC <- list_to_str(BIC, delim='|')
      ret$logLik <- list_to_str(logLik, delim='|')
      if(compute_splithalf_criterion) {
        ret$splithalf_criterion <- list_to_str(splithalf_criterion, delim='|')
      }
      ret$random_model_order <- list_to_str(random_model_order, delim='|')
      ret$delta_observed <- list_to_str(delta_observed, delim='|')
      ret$delta_expected <- delta_expected
      return(ret)
    }
  }
  df <- list_to_df(results)
  export(df, uniquify_file_name(paste0(output_path, '/results_splithalf_trials=', n_trials_per_measure, '.csv')))
}


# # Run this repeatedly and then combine results to one dataframe
# compute_results_df(collections = c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013"),
#                    measures =  c("ap", "p10", "rr", "ndcg20", "err20"),
#                    n_trials_per_measure = 500,
#                    output_path = 'output/copulas',
#                    sort_by_criterion = "AIC",
#                    compute_splithalf_criterion = TRUE,
#                    N_TRIALS = 10)


# # Combine to one dataframe
# combine_results_df('output/copulas/results_splithalf_trials=1000.csv')



# ============================= Plot results =============================
n_total_splits = 250000
output_path = 'output/copulas'
df <- import(paste0(output_path, '/results_splithalf_n', n_total_splits, '.csv'))

# Rename some columns
colnames(df)[colnames(df) == 'splithalf_criterion'] <- 'SHC'
colnames(df)[colnames(df) == 'logLik'] <- 'LL'
# Read columns of strings, as lists
df$model_names <- dfcolumn_str_to_strlist(df$model_names, '|')
df$delta_observed <- dfcolumn_str_to_numlist(df$delta_observed, '|')
df$AIC <- dfcolumn_str_to_numlist(df$AIC, '|')
df$SHC <- dfcolumn_str_to_numlist(df$SHC, '|')
df$BIC <- dfcolumn_str_to_numlist(df$BIC, '|')
df$LL <- dfcolumn_str_to_numlist(df$LL, '|')
df$random_model_order <- dfcolumn_str_to_numlist(df$random_model_order, '|')
df$T1_data_indices <- dfcolumn_str_to_numlist(df$T1_data_indices, '|')
# Calculate T1_data 
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
measures <- c("ap", "p10", "rr", "ndcg20", "err20")
dat <- read_evaluation_data(measures, collections, 0.1)
df$s1_T1_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run1]][x$T1_data_indices])
df$s2_T1_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run2]][x$T1_data_indices])
df$s1_T2_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run1]][-x$T1_data_indices])
df$s2_T2_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run2]][-x$T1_data_indices])
df$pseudoscores_s1_T1 <- dfcolumn_str_to_numlist(df$pseudoscores_s1_T1, '|')
df$pseudoscores_s2_T1 <- dfcolumn_str_to_numlist(df$pseudoscores_s2_T1, '|')
df$pseudoscores_s1_T2 <- dfcolumn_str_to_numlist(df$pseudoscores_s1_T2, '|')
df$pseudoscores_s2_T2 <- dfcolumn_str_to_numlist(df$pseudoscores_s2_T2, '|')
df$m <- df$measure
df$measure <- recode_and_reorder(df$measure) # do this after the calculation of T1_data
df$s1_T1_data_perc_zeros <- sapply(df$s1_T1_data, function(x) length(which(x==0))/length(x))
df$s2_T1_data_perc_zeros <- sapply(df$s2_T1_data, function(x) length(which(x==0))/length(x))
df$s1_T2_data_perc_zeros <- sapply(df$s1_T2_data, function(x) length(which(x==0))/length(x))
df$s2_T2_data_perc_zeros <- sapply(df$s2_T2_data, function(x) length(which(x==0))/length(x))

is_BB1_or_BB6 <- c('BB1', 'Survival BB1', 'Rotated BB1 90 degrees', 'Rotated BB1 180 degrees', 'Rotated BB1 270 degrees',
                   'BB6', 'Survival BB6', 'Rotated BB6 90 degrees', 'Rotated BB6 180 degrees', 'Rotated BB6 270 degrees')

# ================= Best AIC ==================
df$best_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) x[[1]])) # 1st model is the best fit based on AIC
df$best_delta_observed_AIC <- sapply(df$delta_observed, function(x) x[[1]]) # 1st model is the best fit based on AIC
df$GoF_AIC <- -(df$best_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= 2nd best AIC ==================
df$top2_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) get_i(x, 2))) # 2nd model is the top-2 fit based on AIC
df$top2_delta_observed_AIC <- sapply(df$delta_observed, function(x) get_i(x, 2)) # 2nd model is the top-2 fit based on AIC
df$GoF_top2_AIC <- -(df$top2_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= 3rd best AIC ==================
df$top3_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) get_i(x, 3))) # 3rd model is the top-3 fit based on AIC
df$top3_delta_observed_AIC <- sapply(df$delta_observed, function(x) get_i(x, 3)) # 3rd model is the top-3 fit based on AIC
df$GoF_top3_AIC <- -(df$top3_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= Best AIC, except BB1, 6 ==================
df$best_model_name_AIC_no_bb1_6_bestindex <- sapply(df$model_names, function(x) which(!x %in% is_BB1_or_BB6)[[1]])
df$best_model_name_AIC_no_bb1_6  <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_AIC_no_bb1_6_bestindex[[1]]]]))
df$best_delta_observed_AIC_no_bb1_6  <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_AIC_no_bb1_6_bestindex[[1]]]])
df$GoF_AIC_no_bb1_6  <- -(df$best_delta_observed_AIC_no_bb1_6  - df$delta_expected) / df$delta_expected

# ================= Best SHC ==================
df$best_model_name_SHC <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC)[[1]]]]))
df$best_delta_observed_SHC <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC)[[1]]]])
df$GoF_SHC <- -(df$best_delta_observed_SHC - df$delta_expected) / df$delta_expected

# ================= Best BIC ==================
df$best_model_name_BIC <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$BIC)[[1]]]]))
df$best_delta_observed_BIC <- apply(df, 1, function(x) x$delta_observed[[order(x$BIC)[[1]]]])
df$GoF_BIC <- -(df$best_delta_observed_BIC - df$delta_expected) / df$delta_expected

# ================= Best LL ==================
df$best_model_name_LL <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$LL, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL <- apply(df, 1, function(x) x$delta_observed[[order(x$LL, decreasing = TRUE)[[1]]]])
df$GoF_LL <- -(df$best_delta_observed_LL - df$delta_expected) / df$delta_expected

# =============== Random model ===============
df$best_model_name_RAND <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$random_model_order[[1]]]]))
df$best_delta_observed_RAND <- apply(df, 1, function(x) x$delta_observed[[x$random_model_order[[1]]]])
df$GoF_RAND <- -(df$best_delta_observed_RAND - df$delta_expected) / df$delta_expected

# =============== Theoretical Best model (according to lowest Delta observed) ===============
df$best_model_name_ACTUAL <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$delta_observed)[[1]]]]))
df$best_delta_observed_ACTUAL <- apply(df, 1, function(x) x$delta_observed[[order(x$delta_observed)[[1]]]])
df$GoF_ACTUAL <- -(df$best_delta_observed_ACTUAL - df$delta_expected) / df$delta_expected

# =============== Theoretical worst model (according to highest Delta observed) ===============
df$best_model_name_WORST <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$delta_observed, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_WORST <- apply(df, 1, function(x) x$delta_observed[[order(x$delta_observed, decreasing = TRUE)[[1]]]])
df$GoF_WORST <- -(df$best_delta_observed_WORST - df$delta_expected) / df$delta_expected

# =========== Filter ===========
df <- filter(df, GoF_AIC!=-Inf)



# ---------- Print overall means  ----------
df2 <- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))

# aic
mean_df <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val = mean(GoF_AIC))
mean_df

mean_df2 <- df2 %>%
  summarise(mean_val2 = mean(GoF_AIC))
mean_df2

# shc
mean_df <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val = mean(GoF_SHC))
mean_df

mean_df2 <- df2 %>%
  summarise(mean_val2 = mean(GoF_SHC))
mean_df2


# ---------- Plot 1 - (observed, expected) [Best AIC] ----------
df2<- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))

pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig1.pdf'), width=3, height=7.8)
my_colors <- c('red', 'black')
my_shapes <- c(8, 16)
my_labels <- c(bquote(paste(Delta[obs])), bquote(paste(Delta[exp])))
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_point(data=df2, aes(x=delta_expected, y=best_model_name_AIC, color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2,aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_linerange(data=df2,aes(x=delta_expected, y=best_model_name_AIC,color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = '') +
  theme(legend.position = "bottom",
        plot.margin = margin(1, 1, -1, 1, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 2b - GoF [Best AIC no bb1 / 6] ----------
df2<- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))

mean_df <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val = mean(GoF_AIC_no_bb1_6))

pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig2b.pdf'), width=3, height=7.8)
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_vline(data=mean_df, aes(xintercept=mean_val), linetype = "dotted", color='black') +
  geom_point(data=df2, aes(x=GoF_AIC_no_bb1_6, y=best_model_name_AIC_no_bb1_6, color='black'), 
             stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_AIC_no_bb1_6, y=best_model_name_AIC_no_bb1_6, color='black'), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", labels=c(bquote('GoF')))  +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = '') +
  theme(legend.position = "bottom",
        plot.margin = margin(1, 1, -1, 1, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 3 - Frequency that each copula family is selected based on AIC ----------
df2 <- readRDS('data/margin_and_cop_frequencies.rds')
df2 <- filter(df2, name=='Copula')
df2$measure <- recode_and_reorder(df2$measure)
colnames(df2)[colnames(df2) == 'value'] <- 'best_model_name_AIC'
df2 <- df2[c("best_model_name_AIC","measure")]
df2[,'topics'] <- 50
df3 <- df[c("best_model_name_AIC","measure")]
df3[,'topics'] <- 25
df4 <- rbind(df2, df3)

mean_df <- filter(df4, topics==25) %>%
  group_by(measure, best_model_name_AIC) %>%
  summarise(n = n()) %>% 
  mutate(perc = (n / sum(n)*100)) 
mean_df[,'topics'] <- 25

mean_df2 <- filter(df4, topics==50) %>%
  group_by(measure, best_model_name_AIC) %>%
  summarise(n = n()) %>% 
  mutate(perc = (n / sum(n)*100)) 
mean_df2[,'topics'] <- 50

df5 <- rbind(mean_df, mean_df2)

pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig3.pdf'), width=9, height=2.25)
p <- ggplot() +
  geom_bar(data=df5, position="dodge2", stat='identity', aes(x=best_model_name_AIC, y=perc, fill=as.factor(topics), group=measure)) +
  facet_grid(.~measure, space = "free", scale = "free") +
  labs(y= 'Frequency', x = 'Best model (AIC)') +
  scale_fill_manual(name = "", values=c("gray70", "gray30"), labels=c("25 topics","50 topics")) +
  theme(plot.margin = margin(1, -1, 1, 1, "mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p)
dev.off()

filter(df5, best_model_name_AIC %in% c('BB1', 'BB6', 'Student t'))


# ---------- Table 1 - How often is BB1/BB6 selected as the best model by AIC?  ----------
df2 <- readRDS('data/margin_and_cop_frequencies.rds')
df2 <- filter(df2, name=='Copula')
df2$measure <- recode_and_reorder(df2$measure)
colnames(df2)[colnames(df2) == 'value'] <- 'best_model_name_AIC'
df2 <- df2[c("best_model_name_AIC","measure")]
df2[,'topics'] <- 50
df3 <- df[c("best_model_name_AIC","measure")]
df3[,'topics'] <- 25
df4 <- rbind(df2, df3)

mean_df <- filter(df4, topics==25) %>%
  group_by(best_model_name_AIC) %>%
  summarise(n = n()) %>% 
  mutate(perc = (n / sum(n)*100)) 
mean_df[,'topics'] <- 25

mean_df2 <- filter(df4, topics==50) %>%
  group_by(best_model_name_AIC) %>%
  summarise(n = n()) %>% 
  mutate(perc = (n / sum(n)*100)) 
mean_df2[,'topics'] <- 50

df5 <- rbind(mean_df, mean_df2)
print(filter(df5, best_model_name_AIC=='BB1'))
print(filter(df5, best_model_name_AIC=='BB6')) 


# ---------- Plot 4 - Compare selection criteria [Dobs] ----------
df2 <- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))

my_colors <- c('orange', 'gray', 'black', 'blue', 'darkgreen', 'red')
my_shapes <- c(3, 4, 16, 17, 15, 8)
my_labels <- c(bquote('Expectation' ~ (Delta[exp])), 'Best-case', 'AIC', 'BIC', 'LL', 'SHC')
pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig4.pdf'), width=7.5, height=1.6, pointsize = 10)
p1 <- ggplot() +
  
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  
  geom_point(data=df2, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=best_delta_observed_BIC, y=forcats::fct_rev(measure), color=my_colors[4]), 
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_BIC, y=forcats::fct_rev(measure), color=my_colors[4]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=best_delta_observed_LL, y=forcats::fct_rev(measure), color=my_colors[5]), 
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_LL, y=forcats::fct_rev(measure), color=my_colors[5]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=best_delta_observed_SHC, y=forcats::fct_rev(measure), color=my_colors[6]), 
             shape=my_shapes[6], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_SHC, y=forcats::fct_rev(measure), color=my_colors[6]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote(Delta[obs])) +
  # scale_x_continuous(limits = c(-0.29, 0.19), breaks =  c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, -1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11))
print(p1)
dev.off()


# ---------- Plot 4 - Compare selection criteria [GoF] ----------
df2 <- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))

my_colors <- c('gray', 'black', 'blue', 'darkgreen', 'red')
my_shapes <- c(4, 16, 17, 15, 8)
my_labels <- c('Best-case', 'AIC', 'BIC', 'LL', 'SHC')
pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig4b.pdf'), width=5.75, height=1.6, pointsize = 10)
p2 <- ggplot() +
  
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  
  geom_point(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=GoF_BIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_BIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=GoF_LL, y=forcats::fct_rev(measure), color=my_colors[4]), 
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_LL, y=forcats::fct_rev(measure), color=my_colors[4]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df2, aes(x=GoF_SHC, y=forcats::fct_rev(measure), color=my_colors[5]), 
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_SHC, y=forcats::fct_rev(measure), color=my_colors[5]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('GoF')) +
  # scale_x_continuous(limits = c(-0.29, 0.19), breaks =  c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position = 'none')
print(p2)
dev.off()


# ---------- Plot skewness ----------
clean_cop_name <- function(name) {
  name <- gsub("(^S|_?(90|180|270))", "", name)
  name <- gsub("Tawn2", "Tawn", name)
  name
}
recode_measure <- function(x) {
  recode_and_reorder_help(x,
                          AP = "ap",
                          `nDCG@20` = "ndcg20",
                          `ERR@20` = "err20",
                          `P@10` = "p10",
                          RR = "rr")
}
d <- lapply(c("ap", "p10", "rr", "ndcg20", "err20"), function(measure) {
  import(glue("data/skew_data/{measure}.csv")) %>%
    mutate(copula = clean_cop_name(family),
           copula = ifelse(copula == "Tawn", "Tawn", "Other copulas"),
           measure = recode_measure(measure))
}) %>% bind_rows()
jpeg(paste0(output_path, '/skewness.jpg'), width=14*0.66, height=4*0.66, pointsize = 11, res = 400, units = "in")
print(ggplot(d ,
             aes(skew_trec, skew_cop, color = copula)) +
        geom_point(alpha = .1) +
        facet_wrap(~measure, ncol = 5) +
        lims(x = c(-3, 3), y = c(-3, 3)) +
        labs(x = "Skewness of TREC data", y = "Skewness of simulated data", color = "Copula") +
        guides(color = guide_legend(override.aes = list(alpha = 1),
                                    reverse = TRUE)) +
        theme(legend.position = "bottom", 
              legend.title = element_blank(), 
              legend.text = element_text(size = 11))) 
dev.off()


# ---------- Plot 6 - skewness v2 ----------
df2 <- filter(df, delta_expected>0.001)
nrow(filter(df, !delta_expected>0.001))
df2$skewness_T1_T2 <- apply(df2, 1, function(x) x$s1_T1_data - x$s2_T1_data) 
df2$skewness_T1_T2 <- apply(df2, 1, function(x) skewness(x$skewness_T1_T2))   
df2 <- filter(df2, !is.na(skewness_T1_T2))
bins <- seq(-3, 3, 0.5)
df2$skewness_T1_T2_bins <- cut(df2$skewness_T1_T2, breaks=bins)
df2 <- filter(df2, !is.na(skewness_T1_T2_bins))
df3 <-filter(df2, delta_expected < 0.05)

pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig6a.pdf'), width=7, height=4)
my_colors <- c('red', 'black')
my_shapes <- c(8, 16)
my_labels <- c('observed', 'expected')
p <- ggplot() +
  # geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=best_delta_observed_AIC, y=skewness_T1_T2_bins, color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2,aes(x=best_delta_observed_AIC, y=skewness_T1_T2_bins, color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_point(data=df2, aes(x=delta_expected, y=skewness_T1_T2_bins, color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2,aes(x=delta_expected, y=skewness_T1_T2_bins, color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  facet_wrap(.~measure, ncol = 3, nrow = 2, scale = "fixed") +
  labs(x= bquote(Delta), y = 'Skewness') +
  theme(legend.position = "bottom",
        plot.margin = margin(1, 1, 0, 1, "mm"))
print(p)
dev.off()

pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig6b.pdf'), width=7, height=4)
my_colors <- c('black')
my_shapes <- c(16)
p <- ggplot() +
  geom_vline(data=df3, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df3, aes(x=GoF_AIC, y=skewness_T1_T2_bins, color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df3,aes(x=GoF_AIC, y=skewness_T1_T2_bins, color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  facet_wrap(.~measure, ncol = 3, nrow = 2, scale = "fixed") +
  labs(x= 'GoF', y = 'Skewness') +
  theme(legend.position = "none",
        plot.margin = margin(1, 1, 0, 1, "mm"))
print(p)
dev.off()


# ---------- Plot 7 - correlation matrices ----------
mean_df <- df %>%
  group_by(best_model_name_AIC, best_model_name_SHC) %>%
  summarise(n = n())
mean_df
m <- xtabs(n ~ best_model_name_AIC  + best_model_name_SHC, mean_df)
m2<-t(apply(m,1, function(x) x/sum(x)))
pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig7AIC_SHC.pdf'), width=3.5, height=3.5)
corrplot(m2,        # Correlation matrix
         method = "circle", # Correlation plot method
         type = "full",    # Correlation plot style (also "upper" and "lower")
         diag = TRUE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL,
         cl.lim = c(0,1),
         mar = c(0, 2, 2, 0)) 
title(ylab = "Model selected by AIC")
dev.off()

mean_df <- df %>%
  group_by(best_model_name_AIC, best_model_name_BIC) %>%
  summarise(n = n())
mean_df
m <- xtabs(n ~ best_model_name_AIC  + best_model_name_BIC, mean_df)
m2<-t(apply(m,1, function(x) x/sum(x)))
pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_fig7AIC_BIC.pdf'), width=3.5, height=3.5)
corrplot(m2,        # Correlation matrix
         method = "circle", # Correlation plot method
         type = "full",    # Correlation plot style (also "upper" and "lower")
         diag = TRUE,      # If TRUE (default), adds the diagonal
         tl.col = "black", # Labels color
         bg = "white",     # Background color
         title = "",       # Main title
         col = NULL,
         cl.lim = c(0,1),
         mar = c(0, 2, 2, 0)) 
title(ylab = "Model selected by AIC")
dev.off()
