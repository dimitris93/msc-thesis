source("R/utils.R")
set.seed(1234)


# ============================= Plot example Delta calculations (10 plots) =============================
output_path <- 'output/margins'
measures <-  c("ap", "p10", "rr", "ndcg20", "err20")
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")

dir.create(output_path, recursive = TRUE)

dat <- read_evaluation_data(measures, collections, 0.1)
for(measure in measures) {
  # Read evaluation data
  dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
  
  # Sample collection, proportional to number of runs in it
  collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
  run <- sample(names(dat[[measure]][[collection]]), 1) # sample run
  d <- dat[[measure]][[collection]][[run]]
  
  # Slit-half
  T1_data_indices <- sample(seq_len(length(d)),size = floor(0.5*length(d)))
  T1_data_indices <- sort(T1_data_indices)
  T2_data_indices <- setdiff(seq_len(length(d)),T1_data_indices)
  T1_data <- d[T1_data_indices]
  T2_data <- d[T2_data_indices]
  
  # Fit all margins (on T1 data)
  all_margins <- fit_all_margins(T1_data, measure)
  eff_t1 <- all_margins[[1]] # 1st on the list is the best fit
  
  # Compute CDF
  cdf_t1 <- function (X) { return(peff_2(X, eff_t1)) }
  
  # Define ECDFs
  ecdf_t1 <- function (X) { return(sapply(X, function(x) sum(T1_data <= x) / length(T1_data))) }
  ecdf_t2 <- function (X) { return(sapply(X, function(x) sum(T2_data <= x) / length(T2_data))) }
  
  # Compute Deltas
  delta_observed <- cdf_cdf_diff(cdf_t1, ecdf_t2, measure, only_mean = FALSE)
  delta_expected <- cdf_cdf_diff(ecdf_t1, ecdf_t2, measure, only_mean = FALSE)
  
  legend_position <- "bottomright"
  red_line_width <- if (measure %in% c("ap", "ndcg20", "err20")) 1 else 2
  
  for(delta in c('obs', 'exp')) {
    # Compute Deltas
    if(delta == 'obs') {
      label_1 <- bquote(paste(F[1]^{"*"}, ' (', .(beautify(eff_t1$model$type)), ')'))
      plot_data <- delta_observed
    }
    else if(delta == 'exp') {
      label_1 <- expression(F[1])
      plot_data <- delta_expected
    }
    
    pdf(file=paste0(output_path, '/example_D', delta, '_', measure, '.pdf'), width=3, height=3, pointsize = 10)
    par(mar=c(4.2,4.2,1.9,0.4))
    plot(x=plot_data$X, y=plot_data$cdf_1_Y,
         xlab='X', ylab='Cumulative Probability',
         type='l', xlim=c(0,1), ylim=c(0,1))
    title(bquote(bold(paste(Delta[.(delta)], ' = ', .(paste0(round(plot_data$mean, 3)))))))
    axis(side=1, at=seq(0, 1, by=0.2))
    axis(side=2, at=seq(0, 1, by=0.2))
    # Absolute deltas between the 2 curves (CDF_T1, ECDF_T2)
    segments(x0=plot_data$X, y0=plot_data$lowest_Y, x1=plot_data$X, y1=plot_data$highest_Y,
             col="red", lty=1, lwd=red_line_width)
    # CDF_T1
    lines(x=plot_data$X, y=plot_data$cdf_1_Y, col='blue', lty=1, lwd=2)
    # ECDF_T2
    lines(x=plot_data$X, y=plot_data$cdf_2_Y, col='darkgreen', lty=1, lwd=2,)
    legend(x=legend_position, legend=c(label_1, expression(F[2]), 'Abs. differences'),
           col=c('blue', 'darkgreen', 'red'), lty=c(1,1,1), lwd=c(2,2,2))
    dev.off()
  }
  
  dobs <- round(delta_observed$mean, 3)
  dexp <- round(delta_expected$mean, 3)
  GoF <- round(- (dobs-dexp) / dexp, 3)
  
  
  print(paste0('Example: {', beautify(measure), ', ', run, ', ', collection, '}. GoF = ', GoF ,'.'))
}
