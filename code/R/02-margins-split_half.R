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
      # Sample run
      run <- sample(names(dat[[measure]][[collection]]), 1) 
      d <- dat[[measure]][[collection]][[run]]
      
      # Split-half
      T1_data_indices <- sample(seq_len(length(d)),size = floor(0.5*length(d)))
      T1_data <- d[T1_data_indices]
      T2_data <- d[-T1_data_indices]
      
      # Fit all margins (on T1 data)
      margins <- fit_all_margins(T1_data, 
                                 measure,
                                 sort_by_criterion = sort_by_criterion,
                                 compute_splithalf_criterion = compute_splithalf_criterion,
                                 N_TRIALS = N_TRIALS)
      
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
      
      model_names <- sapply(margins, function(x) x$model$type)
      AIC <- sapply(margins, function(x) x$AIC)
      BIC <- sapply(margins, function(x) x$BIC)
      logLik <- sapply(margins, function(x) x$logLik)
      if(compute_splithalf_criterion) {
        splithalf_criterion <- sapply(margins, function(x) x$splithalf_criterion)
      }
      random_model_order <- sample(seq(1, length(margins)))
      
      ret <- NULL
      ret$measure <- measure
      ret$collection <- collection
      ret$run <- run
      ret$T1_data_indices <- list_to_str(T1_data_indices, delim='|')
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


# ============================= Plot results =============================
n_total_splits = 250000
output_path = 'output/margins'

# plot_results <- function(n_trials, output_path) {
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
df$T1_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run]][x$T1_data_indices])
df$T2_data <- apply(df, 1, function(x) dat[[x$measure]][[x$collection]][[x$run]][-x$T1_data_indices])
df$m <- df$measure
df$measure <- recode_and_reorder(df$measure) # do this after the calculation of T1_data
df$T1_data_perc_zeros <- sapply(df$T1_data, function(x) length(which(x==0))/length(x))
df$T1_data_count_zeros <- sapply(df$T1_data, function(x) length(which(x==0)))
df$T2_data_perc_zeros <- sapply(df$T2_data, function(x) length(which(x==0))/length(x))
df$T2_data_count_zeros <- sapply(df$T2_data, function(x) length(which(x==0)))

# ================= Best AIC ==================
df$best_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) x[[1]])) # 1st model is the best fit based on AIC
df$best_delta_observed_AIC <- sapply(df$delta_observed, function(x) x[[1]]) # 1st model is the best fit based on AIC
df$GoF_AIC <- -(df$best_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= Best AIC, except beta ks ==================
df$best_model_name_AIC_no_bks_bestindex <- sapply(df$model_names, function(x) if (x[[1]]=='bks') 2 else 1)
df$best_model_name_AIC_no_bks <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_AIC_no_bks_bestindex[[1]]]]))
df$best_delta_observed_AIC_no_bks <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_AIC_no_bks_bestindex[[1]]]])
df$GoF_AIC_no_bks <- -(df$best_delta_observed_AIC_no_bks - df$delta_expected) / df$delta_expected

# ================= 2nd best AIC ==================
df$top2_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) get_i(x, 2))) # 2nd model is the top-2 fit based on AIC
df$top2_delta_observed_AIC <- sapply(df$delta_observed, function(x) get_i(x, 2)) # 2nd model is the top-2 fit based on AIC
df$GoF_top2_AIC <- -(df$top2_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= 3rd best AIC ==================
df$top3_model_name_AIC <- recode_and_reorder(sapply(df$model_names, function(x) get_i(x, 3))) # 3rd model is the top-3 fit based on AIC
df$top3_delta_observed_AIC <- sapply(df$delta_observed, function(x) get_i(x, 3)) # 3rd model is the top-3 fit based on AIC
df$GoF_top3_AIC <- -(df$top3_delta_observed_AIC - df$delta_expected) / df$delta_expected

# ================= Best SHC ==================
df$best_model_name_SHC <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$SHC)[[1]]]]))
df$best_delta_observed_SHC <- apply(df, 1, function(x) x$delta_observed[[order(x$SHC)[[1]]]])
df$GoF_SHC <- -(df$best_delta_observed_SHC - df$delta_expected) / df$delta_expected

# ================= Best SHC, except beta ks ==================
df$best_model_name_SHC_no_bks_bestindex <- apply(df, 1, function(x) if (x$model_names[[order(x$SHC)[[1]]]] =='bks') order(x$SHC)[[2]] else order(x$SHC)[[1]])
df$best_model_name_SHC_no_bks <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_SHC_no_bks_bestindex[[1]]]]))
df$best_delta_observed_SHC_no_bks <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_SHC_no_bks_bestindex[[1]]]])
df$GoF_SHC_no_bks <- -(df$best_delta_observed_SHC_no_bks - df$delta_expected) / df$delta_expected

# ================= Best BIC ==================
df$best_model_name_BIC <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$BIC)[[1]]]]))
df$best_delta_observed_BIC <- apply(df, 1, function(x) x$delta_observed[[order(x$BIC)[[1]]]])
df$GoF_BIC <- -(df$best_delta_observed_BIC - df$delta_expected) / df$delta_expected

# ================= Best BIC, except beta ks ==================
df$best_model_name_BIC_no_bks_bestindex <- apply(df, 1, function(x) if (x$model_names[[order(x$BIC)[[1]]]] =='bks') order(x$BIC)[[2]] else order(x$BIC)[[1]])
df$best_model_name_BIC_no_bks <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_BIC_no_bks_bestindex[[1]]]]))
df$best_delta_observed_BIC_no_bks <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_BIC_no_bks_bestindex[[1]]]])
df$GoF_BIC_no_bks <- -(df$best_delta_observed_BIC_no_bks - df$delta_expected) / df$delta_expected

# ================= Best LL ==================
df$best_model_name_LL <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$LL, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_LL <- apply(df, 1, function(x) x$delta_observed[[order(x$LL, decreasing = TRUE)[[1]]]])
df$GoF_LL <- -(df$best_delta_observed_LL - df$delta_expected) / df$delta_expected

# ================= Best LL, except beta ks ==================
df$best_model_name_LL_no_bks_bestindex <- apply(df, 1, function(x) if (x$model_names[[order(x$LL, decreasing = T)[[1]]]] =='bks') order(x$LL, decreasing = T)[[2]] else order(x$LL, decreasing = T)[[1]])
df$best_model_name_LL_no_bks <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_LL_no_bks_bestindex[[1]]]]))
df$best_delta_observed_LL_no_bks <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_LL_no_bks_bestindex[[1]]]])
df$GoF_LL_no_bks <- -(df$best_delta_observed_LL_no_bks - df$delta_expected) / df$delta_expected

# =============== Random model ===============
df$best_model_name_RAND <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$random_model_order[[1]]]]))
df$best_delta_observed_RAND <- apply(df, 1, function(x) x$delta_observed[[x$random_model_order[[1]]]])
df$GoF_RAND <- -(df$best_delta_observed_RAND - df$delta_expected) / df$delta_expected

# =============== Theoretical Best model (according to lowest Delta observed) ===============
df$best_model_name_ACTUAL <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$delta_observed)[[1]]]]))
df$best_delta_observed_ACTUAL <- apply(df, 1, function(x) x$delta_observed[[order(x$delta_observed)[[1]]]])
df$GoF_ACTUAL <- -(df$best_delta_observed_ACTUAL - df$delta_expected) / df$delta_expected

# =============== Theoretical Best model, except Beta KS (according to lowest Delta observed) ===============
# df$best_model_name_ACTUAL_no_bks_bestindex <- apply(df, 1, function(x) if (x$model_names[[order(x$delta_observed, decreasing)[[1]]]] =='bks') order(x$delta_observed, decreasing)[[2]] else order(x$delta_observed, decreasing)[[1]])
# df$best_model_name_ACTUAL_no_bks <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[x$best_model_name_ACTUAL_no_bks_bestindex[[1]]]]))
# df$best_delta_observed_ACTUAL_no_bks <- apply(df, 1, function(x) x$delta_observed[[x$best_model_name_ACTUAL_no_bks_bestindex[[1]]]])
# df$GoF_ACTUAL_no_bks <- -(df$best_delta_observed_ACTUAL_no_bks - df$delta_expected) / df$delta_expected

# =============== Theoretical worst model (according to highest Delta observed) ===============
df$best_model_name_WORST <- recode_and_reorder(apply(df, 1, function(x) x$model_names[[order(x$delta_observed, decreasing = TRUE)[[1]]]]))
df$best_delta_observed_WORST <- apply(df, 1, function(x) x$delta_observed[[order(x$delta_observed, decreasing = TRUE)[[1]]]])
df$GoF_WORST <- -(df$best_delta_observed_WORST - df$delta_expected) / df$delta_expected



# ---------- Plot 0a - Visual comparison of marginal candidates ----------
# for (row_id in sample(1:50000, 20, replace=FALSE)) {  # 44739
# for (row_id in sample(150001:200000, 50, replace=FALSE)) {
for (row_id in c(44739, 152379)) {  # 44739, 152379
  row <- df[row_id,]
  T_data <- c(row$T1_data[[1]], row$T2_data[[1]])
  m_effs <- fit_all_margins(T_data, row$m)
  m_pdfs <- sapply(m_effs, function(m_eff) function (X) { return(deff(X, m_eff)) })
  
  if(row$m %in% c("ap", "ndcg20", "err20")) { # continuous case
    X <- seq(0, 1, length.out = 1001)
    measure_type <- 'Continuous'
  } else { # discrete case
    X <- support(row$m)
    measure_type <- 'Discrete'
  }
  
  model_names <- sapply(m_effs, function(m_eff) beautify(m_eff$model$type))
  my_colors <- c('blue', 'darkgreen', 'red', 'purple', 'yellow')
  ltys <- c(1,1,1,1,1)
  lwds <- c(1,1,1,1,1)
  my_ylim <- if (measure_type=='Continuous') c(0,4) else c(0,0.25)
  my_type <- if (measure_type=='Continuous') 'l' else 'o' 
  my_ylab <- if (measure_type=='Continuous') 'Density' else 'Mass' 
  
  if(measure_type=='Continuous' && length(model_names) !=4 ||
     measure_type=='Discrete' && length(model_names) !=3) {
    next
  }
  
  print(measure_type)
  print(row$collection)
  print(row$run)
  print(sort(T_data))
  print('--------')
  
  pdf(file=paste0(output_path, '/families_row', row_id, '_', row$m, '.pdf'), width=3, height=3, pointsize = 10)
  par(mar=c(4.2,4.2,1.9,0.4))
  if(measure_type=='Continuous') {
    hist(x=T_data, probability = TRUE, main = measure_type, xlab = row$measure, ylab='Density', border = "darkgrey",
         xlim = c(0,1), ylim=my_ylim)
  }
  else {
    plot(x=X, y=tabulate(T_data*10+1,11) / length(T_data), main = measure_type, xlab = row$measure, ylab=my_ylab, border = "darkgrey",
         xlim = c(0,1), ylim=my_ylim, pch = 19, col = "darkgrey")
    lines(x=X, y=tabulate(T_data*10+1,11) / length(T_data), col='darkgrey', type=my_type, lwd=1, lty=2)
  }
  
  for (i in 1:length(model_names)){
    lines(x=X, y=m_pdfs[[i]](X), col=my_colors[i], type=my_type, lwd=1)
  }
  box()
  legend(x='topright', legend=model_names, col=my_colors, lty=ltys, lwd=lwds, bty = "n")
  dev.off()
}


# ---------- Plot 0b - Visual comparison of copula candidates ----------
row1 <- df[44739,]
row2 <- df[152379,]
scores_s1 <- read.csv(paste0("data/", row1$collection,"_",row1$m,".csv"))[,c(row1$run)]  
scores_s2 <- read.csv(paste0("data/", row1$collection,"_",row1$m,".csv"))[,c(row2$run)]  
u <- pobs(cbind(scores_s1, scores_s2), ties.method = "random") # compute pseudo-observations
pseudoscores_s1 <- u[,1]
pseudoscores_s2 <- u[,2]
n <- 100
xy_range <- seq(0, 1, length.out = n)
cops <- fit_all_copulas(pseudoscores_s1, pseudoscores_s2, sort_by_criterion = "AIC") # fit copula

cop_1 = cops[[1]]
c_cdf_1 <- function (x, y) { return(BiCopCDF(x, y, cop_1)) }
Z_1 <- sapply(xy_range, function(x) c_cdf_1(xy_range, rep(x, n)))

cop_2 = cops[[14]]
c_cdf_2 <- function (x, y) { return(BiCopCDF(x, y, cop_2)) }
Z_2 <- sapply(xy_range, function(x) c_cdf_2(xy_range, rep(x, n)))

BiCopPDF(0.234, 0.523, cop_2)
BiCopPDF(qnorm(0.234), qnorm(0.523), cop_2)

cop_3 = cops[[21]]
c_cdf_3 <- function (x, y) { return(BiCopCDF(x, y, cop_3)) }
Z_3 <- sapply(xy_range, function(x) c_cdf_3(xy_range, rep(x, n)))

c_cdf_data_1 <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z_1)))
c_cdf_data_2 <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z_2)))
c_cdf_data_3 <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z_3)))
ecdf <- function (x, y) {return(sapply(1:length(x), function(i) mean(pseudoscores_s1 <= x[i] & pseudoscores_s2 <= y[i])))}
Z <- sapply(xy_range, function(x) ecdf(xy_range, rep(x, n)))
ecdf_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z))) 

print(cop_1$familyname)
print(cop_2$familyname)
print(cop_3$familyname)

# Contour Gaussian vs Tawn 1
pdf(file=paste0('output/copulas/families_copula_', row1$m, '.pdf'), width=3, height=3, pointsize = 10)
par(mar=c(4.2,4.2,0.4,0.4))
plot(qnorm(pseudoscores_s1), qnorm(pseudoscores_s2), col = "darkgrey", pch = 19, xlim = c(-3,3), ylim = c(-3,3),
     xlab = "System 52", ylab = "System 29", main = "")
axis(1, at=-3:3)
axis(2, at = -3:3)
box()
contour(cop_1, drawlabels = FALSE, col = 'red', add = TRUE)
contour(cop_2, drawlabels = FALSE, col = 'blue', add = TRUE)
legend(x='topleft', bty = "n", col = c('red','blue'), lwd = 1, legend = c(beautify2(cop_1$familyname), cop_2$familyname))
dev.off()

# Contour Independence
pdf(file=paste0('output/copulas/families_copula_indep_', row1$m, '.pdf'), width=3, height=3, pointsize = 10)
par(mar=c(4.2,4.2,0.4,0.4))
plot(qnorm(pseudoscores_s1), qnorm(pseudoscores_s2), col = "darkgrey", pch = 19, xlim = c(-3,3), ylim = c(-3,3),
     xlab = "System 52", ylab = "System 29", main = "")
axis(1, at=-3:3)
axis(2, at = -3:3)
box()
contour(cop_3, drawlabels = FALSE, col = 'darkgreen', add = TRUE)
legend(x='topleft', bty = "n", col = c('darkgreen'), lwd = 1, legend = c(cop_3$familyname))
dev.off()


# # testing "coupled with standard norm"
# c_pdf_1 <- function (x, y) { return(BiCopPDF(x, y, cop_1)) }
# xy_range_actual <- seq(-3, 3, length.out = n)
# xy_range_actual
# xy_range_pseudo <- pnorm(xy_range_actual)
# xy_range_pseudo
# c_PDF_vals_on_pseudo <- sapply(xy_range_pseudo, function(x) c_pdf_1(xy_range_pseudo, rep(x, n)))
# 
# contour(x = xy_range_actual,
#         y = xy_range_actual,
#         z = c_PDF_vals_on_pseudo*tcrossprod(dnorm(xy_range_actual)))
# 
# plot_ly(x=xy_range_actual, y=xy_range_actual, z=c_PDF_vals_on_pseudo*tcrossprod(dnorm(xy_range_actual)), type = "contour")


# # 3D
# fnt1 <- list(family = "Arial", size = 25)
# fnt2 <- list(family = "Arial", size = 15)
# p_1_sample <- plot_ly(width=1000, height=1000, showscale=F, showlegend=T) %>%
#   # config(mathjax = 'cdn')  %>%
#   add_trace(x=c_cdf_data_1$x, y=c_cdf_data_1$y, z=c_cdf_data_1$z, intensity=c_cdf_data_1$z, name=cop_1$familyname,
#             type="mesh3d", opacity=0.8, colorscale = list(c(0,1),c("red","red"))) %>%
#   add_trace(x=c_cdf_data_2$x, y=c_cdf_data_2$y, z=c_cdf_data_2$z, intensity=c_cdf_data_2$z, name=cop_2$familyname,
#             type="mesh3d", opacity=1, colorscale = list(c(0,1),c("blue","blue"))) %>%
#   add_trace(x=ecdf_data$x, y=ecdf_data$y, z=ecdf_data$z, intensity=ecdf_data$z, name='Empirical',
#             type="mesh3d", opacity=0.8, colorscale = list(c(0,1),c("darkgrey","darkgrey"))) %>%
#   # add_markers(x=ecdf_data$x, y=ecdf_data$y, z=ecdf_data$z,   name='',
#   #             opacity=1, marker=list(size=1, color='darkgrey')) %>%
#   layout(scene = list(xaxis=list(title=list(text='System 52', font=fnt1, standoff = 150), showticklabels = T, tickfont = fnt2),
#                       yaxis=list(title=list(text='System 29', font=fnt1), showticklabels = T, tickfont = fnt2),
#                       zaxis=list(title=list(text='Cumulative Probability', font=fnt1), showticklabels = T, tickfont = fnt2),
#                       camera=list(eye = list(x = 4, y = 1.25, z = 1.25))))
# # layout(scene = list(xaxis=list(title=list(text='System 52'), showticklabels = T),
# #                     yaxis=list(title=list(text='System 29'), showticklabels = T),
# #                     zaxis=list(title=list(text='Cumulative Probability'), showticklabels = T),
# #                     camera=list(eye = list(x = 4, y = 1.25, z = 1.25)))) 
# print(p_1_sample)


# # Attempt to use rgl, etc...
# dir.create('output/3d-copula-plots', recursive = TRUE)
# write.csv(c_cdf_data_1,'output/3d-copula-plots/vis-comp-joe.csv')
# write.csv(c_cdf_data_2,'output/3d-copula-plots/vis-comp-gaus.csv')
# write.csv(ecdf_data,'output/3d-copula-plots/vis-comp-ecdf.csv')
# 
# library(orientlib)
# library(rgl)
# open3d()
# view3d(theta = -45, phi = 25)
# axes3d( edges=c("x--", "y--", "z-+") )
# # axis3d('x', pos = c(0, 0, 0))
# persp3d(x=xy_range, y=xy_range, z=Z_1, alpha=0.5,
#       main="Perspective Plot of a Cone",
#       zlab = "Height",
#       col = "red", shade = 0.001, add = TRUE, axes=FALSE)
# persp3d(x=xy_range, y=xy_range, z=Z_1, alpha=0.5,
#         front = "lines", back = "lines", 
#         lit = FALSE, add = TRUE, axes=FALSE)
# persp3d(x=xy_range, y=xy_range, z=Z_2,
#         main="Perspective Plot of a Cone",
#         zlab = "Height",
#         col = "blue", shade = 0.001, add = TRUE, axes=FALSE)
# persp3d(x=xy_range, y=xy_range, z=Z_2, 
#         front = "lines", back = "lines", 
#         lit = FALSE, add = TRUE, axes=FALSE)
# # title3d('', '', 'X', 'Y')
# mtext3d("X", "x--", line=6, dist=2)
# mtext3d("Y", "y--", line=6)
# mtext3d("Z", "z-+", line=4)
#
# rglToBase(par3d(no.readonly=TRUE)$userMatrix)


# ---------- Plot 1 - (observed, expected) [Best AIC] ----------
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig1.pdf'), width=3.25, height=4.65)
my_colors <- c('red', 'black')
my_shapes <- c(8, 16)
my_labels <- c(bquote(paste(Delta[obs])), bquote(paste(Delta[exp])))
p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df, aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_point(data=df, aes(x=delta_expected, y=best_model_name_AIC, color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df,aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_linerange(data=df,aes(x=delta_expected, y=best_model_name_AIC,color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = '') +
  theme(legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, -0.5, 1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 2 - GoF [Best AIC] ----------
mean_df <- df %>%
  group_by(measure) %>%
  summarise(mean_val = mean(GoF_AIC))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig2.pdf'), width=3.25, height=4.65)
p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df, aes(x=GoF_AIC, y=best_model_name_AIC, color='black'), 
             stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_AIC, y=best_model_name_AIC, color='black'), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_vline(data=mean_df, aes(xintercept=mean_val), linetype = "dotted", color='black') +
  scale_color_identity(name = "", guide = "legend", labels=c(bquote('GoF')))  +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = '') +
  theme(legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, -0.5, 1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 3 - Frequency that each marginal family is selected based on AIC ----------
df2 <- readRDS('data/margin_and_cop_frequencies.rds')
df2 <- filter(df2, name=='Margin')
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

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig3.pdf'), width=6.5, height=2.25)
p <- ggplot() +
  geom_bar(data=df5, position="dodge2", stat='identity', aes(x=best_model_name_AIC, y=perc, fill=as.factor(topics), group=measure)) +
  facet_grid(.~measure, space = "free", scale = "free") +
  labs(y= 'Frequency', x = 'Best model (AIC)') +
  scale_fill_manual(name = "", values=c("gray70", "gray30"), labels=c("25 topics","50 topics")) +
  theme(plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
print(p)
dev.off()


# ---------- Table 1 - How often is Beta KS selected as the best model by AIC?  ----------
df2 <- readRDS('data/margin_and_cop_frequencies.rds')
df2 <- filter(df2, name=='Margin')
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

nrow(filter(df2, measure=='AP'))

print(filter(df5, best_model_name_AIC=='Beta KS'))


#  ---------------- Plot 4 - only bks - overall top1 vs top2 vs top 3 ----------------
df2 <- filter(df, best_model_name_AIC=='Beta KS')
df3 <- filter(df2, !is.na(top2_model_name_AIC))
df4 <- filter(df3, !is.na(top3_model_name_AIC))
my_colors <- c('red', 'blue', 'black', 'gray')
my_shapes <- list(8, 17, 16, 4)
my_lines <- c('solid', 'solid', 'solid', 'solid')
my_labels <- c('Beta KS', 'top-2 AIC', 'top-3 AIC', 'Best-case')

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig4.pdf'), width=4.6, height=1.5)
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  # Beta KS
  geom_point(data=df2, aes(x=best_delta_observed_AIC, y=measure, color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_AIC, y=measure, color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) + 
  # Top-2 AIC
  geom_point(data=df3, aes(x=top2_delta_observed_AIC, y=measure, color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=top2_delta_observed_AIC, y=measure, color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  # Top-3 AIC
  geom_point(data=df4, aes(x=top3_delta_observed_AIC, y=measure, color=my_colors[3]),
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df4, aes(x=top3_delta_observed_AIC, y=measure, color=my_colors[3]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  # Best case
  geom_point(data=df2, aes(x=best_delta_observed_ACTUAL, y=measure, color=my_colors[4]),
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_ACTUAL, y=measure, color=my_colors[4]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
  labs(y = bquote(paste('')), x = bquote(Delta[obs])) +
  # scale_x_continuous(limits = c(-1.6, 0.3), breaks =  seq(-1.5, 0, 0.5), oob = scales::oob_keep) +
  scale_y_discrete(limits=rev) +
  theme(legend.position = "right",
        plot.margin = margin(1.15, -16, 0, -1, "mm"),
        legend.margin = margin(-10, 0, 0, 0),
        legend.box.margin = margin(-5, 50, 0, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


#  ---------------- Plot 5 - only bks - detailed top1 vs top2 vs top 3 ----------------
df2 <- filter(df, best_model_name_AIC=='Beta KS')
df3 <- filter(df2, !is.na(top2_model_name_AIC))
df4 <- filter(df3, !is.na(top3_model_name_AIC))
my_colors <- c('red', 'blue', 'black')
my_shapes <- list('', 17, 16)
my_lines <- c('solid', 'solid', 'solid')
my_labels <- c('Beta KS', 'top-2 AIC', 'top-3 AIC')

mean_df <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val = mean(best_delta_observed_AIC))

mean_df2 <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val2 = mean(best_delta_observed_ACTUAL))

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig5.pdf'), width=5, height=2.8)
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  # Beta KS
  geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[1]), linetype = my_lines[1]) +
  # Top-2 AIC
  geom_point(data=df3, aes(x=top2_delta_observed_AIC, y=top2_model_name_AIC, color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=top2_delta_observed_AIC, y=top2_model_name_AIC, color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  # Top-3 AIC
  geom_point(data=df4, aes(x=top3_delta_observed_AIC, y=top3_model_name_AIC, color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df4, aes(x=top3_delta_observed_AIC, y=top3_model_name_AIC, color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
  facet_grid(measure~., scale = "free", space= "free") +
  labs(y = bquote(paste('Alternative model')), x = bquote(Delta[obs])) +
  # scale_x_continuous(limits = c(-1.6, 0.3), breaks =  seq(-1.5, 0, 0.5), oob = scales::oob_keep) +
  theme(legend.position = "right",
        plot.margin = margin(1, -16, 0, 1, "mm"),
        legend.margin = margin(-10, 0, 0, 0),
        legend.box.margin = margin(-5, 50, 0, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()


#  ---------------- Plot 6a-b-c - Exclude bks [dobs, dexp] ----------------
#  ======  a  ====== 
df2 <- filter(df, !is.na(best_model_name_AIC_no_bks))
df2 <- filter(df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))
df3 <- filter(df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))

mean_df <- df %>%
  group_by(measure) %>%
  summarise(mean_val = mean(best_delta_observed_AIC))
mean_df <- filter(mean_df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))
mean_df2 <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val2 = mean(best_delta_observed_AIC_no_bks))

my_colors <- c('black', 'blue')
my_shapes <- list(16, 8)
my_lines <- c('solid', 'solid')
my_labels <- c('Including Beta KS', 'Excluding Beta KS')

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig6a.pdf'), width=3, height=3.6)
p1 <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[1]), linetype = 'solid') +
  geom_vline(data=mean_df2, aes(xintercept=mean_val2), linetype = 'solid', color=my_colors[2]) +
  
  # to keep order
  geom_point(data=df3, aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]),
             shape=16, stat="summary", fun="mean") +
  
  geom_point(data=df2, aes(x=best_delta_observed_AIC_no_bks, y=best_model_name_AIC_no_bks, color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_AIC_no_bks, y=best_model_name_AIC_no_bks, color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df3, aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]),
             shape=16, stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=best_delta_observed_AIC, y=best_model_name_AIC, color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  # scale_x_continuous(limits = c(-1.6, 0.25), breaks =  seq(-1.5, 0, 0.25), oob = scales::oob_keep) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = bquote(Delta[obs])) +
  theme(legend.position = "bottom",
        legend.direction = 'vertical',
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
        legend.margin = margin(-15, 0, 0, 0),
        # legend.box.margin = margin(-25, 0, 15, 0),
        legend.text = element_text(size = 11))
print(p1)
dev.off()

#  ======  b  ======  
df2 <- filter(df, !is.na(best_model_name_AIC_no_bks))
df2 <- filter(df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))
df3 <- filter(df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))

mean_df <- df %>%
  group_by(measure) %>%
  summarise(mean_val = mean(GoF_AIC))
mean_df <- filter(mean_df, measure %in%c('AP', 'nDCG@20', 'ERR@20'))
mean_df2 <- df2 %>%
  group_by(measure) %>%
  summarise(mean_val2 = mean(GoF_AIC_no_bks))

print(mean_df)
print(mean_df2)

my_colors <- c('black', 'blue')
my_shapes <- list(16, 8)
my_lines <- c('solid', 'solid')
my_labels <- c('Including Beta KS', 'Excluding Beta KS')

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig6b.pdf'), width=3, height=3.6)
p2 <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_vline(data=mean_df, aes(xintercept=mean_val, color=my_colors[1]), linetype = 'solid') +
  geom_vline(data=mean_df2, aes(xintercept=mean_val2), linetype = 'solid', color=my_colors[2]) +
  
  # to keep order
  geom_point(data=df3, aes(x=GoF_AIC, y=best_model_name_AIC, color=my_colors[1]),
             shape=16, stat="summary", fun="mean") +
  
  geom_point(data=df2, aes(x=GoF_AIC_no_bks, y=best_model_name_AIC_no_bks, color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_AIC_no_bks, y=best_model_name_AIC_no_bks, color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  
  geom_point(data=df3, aes(x=GoF_AIC, y=best_model_name_AIC, color=my_colors[1]),
             shape=16, stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=GoF_AIC, y=best_model_name_AIC, color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  # scale_x_continuous(limits = c(-1.6, 0.25), breaks =  seq(-1.5, 0, 0.25), oob = scales::oob_keep) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(shape = my_shapes, linetype = my_lines))) +
  facet_grid(measure~., space = "free", scale = "free") +
  labs(y= 'Best model (AIC)', x = 'GoF') +
  theme(legend.position = "bottom",
        legend.direction = 'vertical',
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, "mm"),
        legend.margin = margin(-15, 0, 0, 0),
        # legend.box.margin = margin(-25, 0, 15, 0),
        legend.text = element_text(size = 11))
print(p2)
dev.off()

#  ======  c  ====== 
p4 <- ggarrange(p1, p2, common.legend = TRUE, legend = "bottom", widths = c(3.25,3.25), heights = c(3.6,3.6), ncol = 2, align = 'h') +
  theme(plot.margin = margin(0,0,1.5,0, "mm"))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig6c.pdf'), width=3.25*2, height=3.6, pointsize = 10)
print(p4)
dev.off()


# ---------- Plot - explore Beta KS: good vs bad fits ----------
df2 <- filter(df, best_model_name_AIC=='Beta KS')
df2 <- filter(df2, !is.na(top2_model_name_AIC))
df2 <- filter(df2, !is.na(top3_model_name_AIC))
df2 <- df2[order(df2$GoF_AIC),] # sort data by lowest GoF

output_path_2 <- paste0(output_path, '/explore_bks_performance/')
dir.create(output_path_2, recursive = TRUE)

n = 30
for(i in c(1:n, nrow(df2):(nrow(df2)-n))) {
  # Select row
  T1_data <- df2[i,]$T1_data[[1]]
  T2_data <- df2[i,]$T2_data[[1]]
  m <- df2[i,]$m
  GoF_AIC <- df2[i,]$GoF_AIC

  all_margins <- fit_all_margins(T1_data, m)
  top_1_margin <- all_margins[[1]] # 1st on the list is the best fit
  top_2_margin <- all_margins[[2]]
  top_3_margin <- all_margins[[3]]
  cdf_top_1_margin <- function (X) { return(peff_2(X, top_1_margin)) }
  cdf_top_2_margin <- function (X) { return(peff_2(X, top_2_margin)) }
  cdf_top_3_margin <- function (X) { return(peff_2(X, top_3_margin)) }
  ecdf_t1 <- function (X) { return(sapply(X, function(x) sum(T1_data <= x) / length(T1_data))) }
  ecdf_t2 <- function (X) { return(sapply(X, function(x) sum(T2_data <= x) / length(T2_data))) }
  delta_observed_top_1_margin <- cdf_cdf_diff(cdf_top_1_margin, ecdf_t2, m, only_mean=F)
  delta_observed_top_2_margin <- cdf_cdf_diff(cdf_top_2_margin, ecdf_t2, m, only_mean=F)
  delta_observed_top_3_margin <- cdf_cdf_diff(cdf_top_3_margin, ecdf_t2, m, only_mean=F)
  delta_expected <- cdf_cdf_diff(ecdf_t1, ecdf_t2, m, only_mean=F)

  pdf(file=paste0(output_path_2, if (i <= n) 'worst' else 'best_', i, '_', m, '.pdf'),
      width=4.1, height=4.1, pointsize = 10)
  par(mar=c(4.2,4.2,1.9,0.4))
  plot(x=delta_observed_top_1_margin$X, y=delta_observed_top_1_margin$cdf_1_Y,
       xlab='X', ylab='Cumulative Probability',
       type='l', xlim=c(0,1), ylim=c(0,1))
  title(bquote(bold(paste(Delta['obs'], ' = ', .(paste0(round(delta_observed_top_1_margin$mean , 3))),
                          # ' ', .(paste0(round(delta_observed_top_2_margin$mean , 3)))
  ))))
  # title(paste0(if (i <= n) 'Example - Bad Beta KS' else 'Example - Good Beta KS'))
  axis(side=1, at=seq(0, 1, by=0.2))
  axis(side=2, at=seq(0, 1, by=0.2))

  # Absolute deltas between the 2 curves (CDF_T1, ECDF_T2)
  segments(x0=delta_observed_top_1_margin$X, y0=delta_observed_top_1_margin$lowest_Y,
           x1=delta_observed_top_1_margin$X, y1=delta_observed_top_1_margin$highest_Y,
           col='lightgray', lty=1, lwd=1)

  # CDF_T1 top-1
  lines(x=delta_observed_top_1_margin$X, y=delta_observed_top_1_margin$cdf_1_Y, col='red', lty=1, lwd=2)
  # CDF_T1 top-2
  lines(x=delta_observed_top_2_margin$X, y=delta_observed_top_2_margin$cdf_1_Y, col='blue', lty=1, lwd=2)
  # # ECDF_T1
  # lines(x=delta_expected$X, y=delta_expected$cdf_1_Y, col='purple', lty=1, lwd=2,)
  # ECDF_T1
  lines(x=delta_expected$X, y=delta_expected$cdf_1_Y, col='black', lty=1, lwd=2,)
  # ECDF_T2
  lines(x=delta_expected$X, y=delta_expected$cdf_2_Y, col='darkgreen', lty=1, lwd=2,)
  legend(x="bottomright",
         legend=c('',
                  paste0('top-1 AIC (', beautify(top_1_margin$model$type), ')',
                         '\n  AIC=', round(top_1_margin$AIC, 1),
                         # '  BIC =', round(top_1_margin$BIC, 1),
                         '  LL=', round(top_1_margin$logLik, 1), '  edf=', round(top_1_margin$df, 1)),
                  '',
                  paste0('\ntop-2 AIC (', beautify(top_2_margin$model$type), ')',
                         '\n  AIC=', round(top_2_margin$AIC, 1),
                         # '  BIC =', round(top_2_margin$BIC, 1),
                         '  LL=', round(top_2_margin$logLik, 1), '  edf=', round(top_2_margin$df, 1)),
                  as.expression(bquote(F[1])),
                  as.expression(bquote(F[2]))),
         col=c('red',NA , 'blue', NA, 'black', 'darkgreen'), lty=c(1,1,1,1,1,1), lwd=c(2,2,2,2,2,2))
  dev.off()
}


# ---------- Plot 7 - Compare Beta KS vs other cases  ----------
df2 <- filter(df, best_model_name_AIC=='Beta KS')
df3 <- filter(df, best_model_name_AIC!='Beta KS')
df3 <- filter(df3, measure%in%c('AP', 'nDCG@20', 'ERR@20'))

my_colors <- c('red', 'black')
my_shapes <- c(8, 16)
my_labels <- c('Corner cases', 'Other cases')

p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=T1_data_perc_zeros, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=T1_data_perc_zeros, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_point(data=df3, aes(x=T1_data_perc_zeros, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=T1_data_perc_zeros, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('Percentage of zeros in the data')) +
  # scale_x_continuous(limits = c(-0.5, 0.2), breaks =  seq(-0.5, 0.2, 0.25), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 4.5, 1.5, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position = 'right')
print(p)

p2 <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_point(data=df3, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('Best-case ' ~ Delta[obs])) +
  # scale_x_continuous(limits = c(-0.45, 0.175), breaks =  seq(-0.5, 0.2, 0.1), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 4.5, 1.5, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position='none')

p3 <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_point(data=df3, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df3, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('Best-case GoF')) +
  # scale_x_continuous(limits = c(-0.45, 0.175), breaks =  seq(-0.5, 0.2, 0.1), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 4.5, 1.5, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position='none')

p4 <- ggarrange(p, p2, p3, common.legend = TRUE, legend = "right", widths = c(4.5,4.5,4.5), heights = c(1.2,1.2,1.2), nrow = 3, align = "v") +
  theme(plot.margin = margin(0,0,-1,0, "mm"))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig7.pdf'), width=4.5, height=3.6, pointsize = 10)
print(p4)
dev.off()


# ---------- Plot 8 - Motivate comparison of selection criteria  ----------
my_colors <- c('black', 'blue', 'gray', 'red')
my_shapes <- c(16, 17, 4, 8)
my_labels <- c('AIC', 
               'AIC (Excluding Beta KS)', 
               'Best-case',
               bquote('Expectation (' * Delta[exp] * ')'))

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig8.pdf'), width=7.5, height=1.6, pointsize = 10)
p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  
  geom_point(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[4]), 
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[4]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote(Delta[obs])) +
  # scale_x_continuous(limits = c(-0.45, 0.175), breaks =  seq(-0.5, 0.2, 0.1), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, -1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 8b - Motivate comparison of selection criteria  ----------
my_colors <- c('black', 'blue', 'gray')
my_shapes <- c(16, 17, 4)
my_labels <- c('AIC', 
               'AIC (Excluding Beta KS)',
               'Best-case')

pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig8b.pdf'), width=5.34, height=1.6, pointsize = 10)
p <- ggplot() +
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  
  geom_point(data=df, aes(x=GoF_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = 'GoF') +
  # scale_x_continuous(limits = c(-0.45, 0.175), breaks =  seq(-0.5, 0.2, 0.1), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position = 'none')
print(p)
dev.off()


# ---------- Plot 9 - Compare selection criteria  ----------
my_colors <- c('gray', 'black', 'blue', 'darkgreen', 'red')
my_shapes <- c(4, 16, 17, 15, 8)
my_labels <- c('Best-case', 'AIC', 'BIC', 'LL', 'SHC')
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig9a_gof.pdf'), width=6, height=1.6, pointsize = 10)
p <- ggplot() +

  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +

  geom_point(data=df, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_AIC, y=forcats::fct_rev(measure), color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=GoF_BIC, y=forcats::fct_rev(measure), color=my_colors[3]),
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_BIC, y=forcats::fct_rev(measure), color=my_colors[3]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=GoF_LL, y=forcats::fct_rev(measure), color=my_colors[4]),
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_LL, y=forcats::fct_rev(measure), color=my_colors[4]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df, aes(x=GoF_SHC, y=forcats::fct_rev(measure), color=my_colors[5]),
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=GoF_SHC, y=forcats::fct_rev(measure), color=my_colors[5]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('GoF')) +
  scale_x_continuous(limits = c(-0.29, 0.19), breaks =  c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, -1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 9b - Compare selection criteria (without BKS) ---------
df2 <- filter(df, measure%in%c('AP', 'nDCG@20', 'ERR@20'))
my_colors <- c('gray', 'black', 'blue', 'darkgreen', 'red')
my_shapes <- c(4, 16, 17, 15, 8)
my_labels <- c('Best-case', 'AIC', 'BIC', 'LL', 'SHC')
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig9b_gof.pdf'), width=4.78, height=1.15, pointsize = 10)
p <- ggplot() +

  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +

  geom_point(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=GoF_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=GoF_BIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[3]),
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_BIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[3]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=GoF_LL_no_bks, y=forcats::fct_rev(measure), color=my_colors[4]),
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_LL_no_bks, y=forcats::fct_rev(measure), color=my_colors[4]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=GoF_SHC_no_bks, y=forcats::fct_rev(measure), color=my_colors[5]),
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_SHC_no_bks, y=forcats::fct_rev(measure), color=my_colors[5]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote('GoF')) +
  scale_x_continuous(limits = c(-0.29, 0.19), breaks =  c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, 1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11),
        legend.position = 'none')
print(p)
dev.off()


# ---------- Plot 9a - Compare selection criteria (without BKS) ---------
my_colors <- c('gray', 'black', 'blue', 'purple', 'darkgreen', 'red')
my_shapes <- c(4, 16, 17, 18, 15, 8)
my_labels <- c('Best-case', 'AIC', 'BIC', 'LL', 'SHC', bquote('Expectation' ~ (Delta[exp])))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig9a.pdf'), width=7.5, height=1.6, pointsize = 10)
p <- ggplot() +
  
  geom_vline(data=df, aes(xintercept=0), linetype = "dashed") +
  
  geom_point(data=df, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[2]), 
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_AIC, y=forcats::fct_rev(measure), color=my_colors[2]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_BIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_BIC, y=forcats::fct_rev(measure), color=my_colors[3]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_LL, y=forcats::fct_rev(measure), color=my_colors[4]), 
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_LL, y=forcats::fct_rev(measure), color=my_colors[4]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=best_delta_observed_SHC, y=forcats::fct_rev(measure), color=my_colors[5]), 
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=best_delta_observed_SHC, y=forcats::fct_rev(measure), color=my_colors[5]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  geom_point(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[6]), 
             shape=my_shapes[6], stat="summary", fun="mean") +
  geom_linerange(data=df, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[6]), 
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote(Delta[obs])) +
  scale_x_continuous(limits = c(-0, 0.084), breaks =  c(0, 0.02, 0.04, 0.06, 0.08), oob = scales::oob_keep) +
  theme(plot.margin = margin(1, -1, 0, -3.5, "mm"),
        axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
        legend.text = element_text(size = 11))
print(p)
dev.off()


# ---------- Plot 9b - Compare selection criteria (without BKS) ---------
df2 <- filter(df, measure%in%c('AP', 'nDCG@20', 'ERR@20'))
my_colors <- c('gray', 'black', 'blue', 'purple', 'darkgreen', 'red')
my_shapes <- c(4, 16, 17, 18, 15, 8)
my_labels <- c('Best-case', 'AIC', 'BIC', 'LL', 'SHC', bquote('Expectation' ~ (Delta[exp])))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_fig9b.pdf'), width=5.75, height=1.15, pointsize = 10)
p <- ggplot() +

  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +

  geom_point(data=df2, aes(x=best_delta_observed_ACTUAL , y=forcats::fct_rev(measure), color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_ACTUAL, y=forcats::fct_rev(measure), color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=best_delta_observed_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_AIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=best_delta_observed_BIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[3]),
             shape=my_shapes[3], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_BIC_no_bks, y=forcats::fct_rev(measure), color=my_colors[3]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=best_delta_observed_LL_no_bks, y=forcats::fct_rev(measure), color=my_colors[4]),
             shape=my_shapes[4], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_LL_no_bks, y=forcats::fct_rev(measure), color=my_colors[4]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=best_delta_observed_SHC_no_bks, y=forcats::fct_rev(measure), color=my_colors[5]),
             shape=my_shapes[5], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_SHC_no_bks, y=forcats::fct_rev(measure), color=my_colors[5]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  geom_point(data=df2, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[6]),
             shape=my_shapes[6], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=delta_expected, y=forcats::fct_rev(measure), color=my_colors[6]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +

  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= '', x = bquote(Delta[obs])) +
  scale_x_continuous(limits = c(-0, 0.084), breaks =  c(0, 0.02, 0.04, 0.06, 0.08), oob = scales::oob_keep) +
    theme(plot.margin = margin(1, 1, 0, -3.5, "mm"),
          axis.text.x = element_text(angle = 0, vjust = 0.5, hjust = 0.5),
          legend.text = element_text(size = 11),
          legend.position = 'none')
  print(p)
  dev.off()


# ---------- Plot a1, a2 - Correlate zeros with Dobs, GoF ---------
# df <- sample_n(df, 2000)
  
df2 <- filter(df, T1_data_count_zeros <13)
df2$T1_data_count_zeros <-
  factor(df2$T1_data_count_zeros,
         levels = as.character(sort(unique(df2$T1_data_count_zeros))))

my_colors <- c('red', 'black')
my_shapes <- c(8, 16)
my_labels <- c(bquote(paste(Delta[obs])), bquote(paste(Delta[exp])))
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_figa1.pdf'), width=3, height=2.5, pointsize = 10)
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=best_delta_observed_AIC , y=T1_data_count_zeros, color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=best_delta_observed_AIC, y=T1_data_count_zeros, color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  geom_point(data=df2, aes(x=delta_expected, y=T1_data_count_zeros, color=my_colors[2]),
             shape=my_shapes[2], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=delta_expected, y=T1_data_count_zeros, color=my_colors[2]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= 'Number of zeros in the data', x = '') +
  # scale_x_continuous(limits = c(-0, 0.084), breaks =  c(0, 0.02, 0.04, 0.06, 0.08), oob = scales::oob_keep) +
  theme(legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, -0.5, 1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()

my_colors <- c('black')
my_shapes <- c(16)
my_labels <- c('GoF')
pdf(paste0(output_path, '/splithalf_n', nrow(df), '_figa2.pdf'), width=3, height=2.5, pointsize = 10)
p <- ggplot() +
  geom_vline(data=df2, aes(xintercept=0), linetype = "dashed") +
  geom_point(data=df2, aes(x=GoF_AIC, y=T1_data_count_zeros, color=my_colors[1]),
             shape=my_shapes[1], stat="summary", fun="mean") +
  geom_linerange(data=df2, aes(x=GoF_AIC, y=T1_data_count_zeros, color=my_colors[1]),
                 stat="summary", fun.data="mean_cl_boot", fun.args=list(conf.int=.95)) +
  scale_color_identity(name = "", guide = "legend", breaks = my_colors, labels = my_labels) +
  guides(colour = guide_legend(override.aes = list(shape = my_shapes))) +
  labs(y= 'Number of zeros in the data', x = '') +
  # scale_x_continuous(limits = c(-0, 0.084), breaks =  c(0, 0.02, 0.04, 0.06, 0.08), oob = scales::oob_keep) +
  theme(legend.position = "bottom",
        plot.margin = margin(1.5, 1.5, -0.5, 1.5, "mm"),
        legend.margin = margin(0, 0, 0, 0),
        legend.box.margin = margin(-10, 0, 5, 0),
        legend.text = element_text(size = 11))
print(p)
dev.off()

# # Run this repeatedly and then combine results to one dataframe
# compute_results_df(collections = c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013"),
#                    measures =  c("ap", "p10", "rr", "ndcg20", "err20"),
#                    n_trials_per_measure = 1250,
#                    output_path = 'output/margins',
#                    sort_by_criterion = "AIC",
#                    compute_splithalf_criterion = TRUE,
#                    N_TRIALS = 10)

# # Combine to one dataframe
# combine_results_df('output/margins/results_splithalf_trials=1250.csv')
