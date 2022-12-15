source("R/utils.R")

n_total_splits = 250000
output_path = 'output/margins'
df <- import(paste0(output_path, '/results_splithalf_n', n_total_splits, '.csv'))
colnames(df)[colnames(df) == 'splithalf_criterion'] <- 'SHC'
colnames(df)[colnames(df) == 'logLik'] <- 'LL'
df$model_names <- dfcolumn_str_to_strlist(df$model_names, '|')
df$delta_observed <- dfcolumn_str_to_numlist(df$delta_observed, '|')
df$AIC <- dfcolumn_str_to_numlist(df$AIC, '|')
df$SHC <- dfcolumn_str_to_numlist(df$SHC, '|')
df$BIC <- dfcolumn_str_to_numlist(df$BIC, '|')
df$LL <- dfcolumn_str_to_numlist(df$LL, '|')
df$random_model_order <- dfcolumn_str_to_numlist(df$random_model_order, '|')
df$T1_data_indices <- dfcolumn_str_to_numlist(df$T1_data_indices, '|')
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
measures <- c("ap", "p10", "rr", "ndcg20", "err20")
dat <- read_evaluation_data(measures, collections, 0.1)
df$T1_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run]][x$T1_data_indices])
df$T2_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run]][-x$T1_data_indices])
df$m <- df$measure
df$measure <- recode_and_reorder(df$measure) # do this after the calculation of T1_data


n_values <- c(1, 2, 5, 10, 25)
# ntrials <- 250
########################################
# n_values <- c(1, 2)
ntrials <- 50

indices <- sample(c(1:nrow(df)), ntrials)
df2 <- df[indices,]


j=1
for(i in indices) {
  row <- df[i,]
  
  for(n in n_values) {
    margins <- fit_all_margins(row$T1_data[[1]], 
                               row$m[[1]], 
                               compute_splithalf_criterion = T, 
                               sort_by_criterion = 'AIC',
                               N_TRIALS = n)
    
    shc_values <- sapply(margins, function(x) x$splithalf_criterion)
    shc_values_str <- list_to_str(shc_values, delim='|')
    
    if(n==1) {
      df2[j,'SHC_n1']  <- shc_values_str
    }
    if(n==2) {
      df2[j,'SHC_n2']  <- shc_values_str
    }
    if(n==5) {
      df2[j,'SHC_n5']  <- shc_values_str
    }
    if(n==10) {
      df2[j,'SHC_n10']  <- shc_values_str
    }
    if(n==25) {
      df2[j,'SHC_n25']  <- shc_values_str
    }
    if(n==50) {
      df2[j,'SHC_n50']  <- shc_values_str
    }
    
  }
  
  df2[j,'SHC_n_testing_index']  <- i
  j <- j +1
}


saveRDS(df2, uniquify_file_name(paste0(output_path, '/new22results_splithalf_n', n_total_splits, '_estimate-shc-trials.rds')))

# combine_results_df('output/margins/results_splithalf_n250000_estimate-shc-trials.rds')


# df <- readRDS(paste0(output_path, '/results_splithalf_n250000_estimate-shc-trials.rds'))

df <- df2

df$SHC_n1 <- dfcolumn_str_to_numlist(df$SHC_n1, '|')
df$SHC_n2 <- dfcolumn_str_to_numlist(df$SHC_n2, '|')
df$SHC_n5 <- dfcolumn_str_to_numlist(df$SHC_n5, '|')
df$SHC_n10 <- dfcolumn_str_to_numlist(df$SHC_n10, '|')
df$SHC_n25 <- dfcolumn_str_to_numlist(df$SHC_n25, '|')


df$best_model_name_SHC_n1 <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC_n1)[[1]]]]))
df$best_delta_observed_SHC_n1 <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC_n1)[[1]]]])
df$GoF_SHC_n1 <- -(df$best_delta_observed_SHC_n1 - df$delta_expected) / df$delta_expected

df$best_model_name_SHC_n2 <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC_n2)[[1]]]]))
df$best_delta_observed_SHC_n2 <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC_n2)[[1]]]])
df$GoF_SHC_n2 <- -(df$best_delta_observed_SHC_n2 - df$delta_expected) / df$delta_expected

df$best_model_name_SHC_n5 <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC_n5)[[1]]]]))
df$best_delta_observed_SHC_n5 <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC_n5)[[1]]]])
df$GoF_SHC_n5 <- -(df$best_delta_observed_SHC_n5 - df$delta_expected) / df$delta_expected

df$best_model_name_SHC_n10 <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC_n10)[[1]]]]))
df$best_delta_observed_SHC_n10 <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC_n10)[[1]]]])
df$GoF_SHC_n10 <- -(df$best_delta_observed_SHC_n10 - df$delta_expected) / df$delta_expected

df$best_model_name_SHC_n25 <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC_n25)[[1]]]]))
df$best_delta_observed_SHC_n25 <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC_n25)[[1]]]])
df$GoF_SHC_n25 <- -(df$best_delta_observed_SHC_n25 - df$delta_expected) / df$delta_expected


my_colors <- c( 'red', 'gray', 'black', 'blue', 'darkgreen', 'orange', 'red')
my_shapes <- c(8, 3, 4, 16, 17, 15)
my_labels <- c(bquote('Expectation' ~ (Delta[exp])), 'n=1', 'n=2', 'n=5', 'n=10', 'n=25')


pdf(paste0(output_path, '/splithalf_n', n_total_splits, '_figdeterminenshc.pdf'), width=7.5, height=1.6, pointsize = 10)
p <- ggplot() +

  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +

  geom_point(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=best_delta_observed_SHC_n1, y=forcats::fct_rev(measure), color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC_n1, y=forcats::fct_rev(measure), color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=best_delta_observed_SHC_n2, y=forcats::fct_rev(measure), color=my_colors[3]),
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC_n2, y=forcats::fct_rev(measure), color=my_colors[3]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=best_delta_observed_SHC_n5, y=forcats::fct_rev(measure), color=my_colors[4]),
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC_n5, y=forcats::fct_rev(measure), color=my_colors[4]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=best_delta_observed_SHC_n10, y=forcats::fct_rev(measure), color=my_colors[5]),
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC_n10, y=forcats::fct_rev(measure), color=my_colors[5]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=best_delta_observed_SHC_n25, y=forcats::fct_rev(measure), color=my_colors[6]),
             shape=my_shapes[6], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC_n25, y=forcats::fct_rev(measure), color=my_colors[6]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote(Delta[obs])) +
  # scale_x_continuous(limits = c(-0.29, 0.19), breaks =  c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, -1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11))
print(p)
dev.off()




