# ==================================================================================== 
# Plots for diagrams 1, 2 & 3
# ==================================================================================== 
source("R/utils.R")
set.seed(8235238)
output_path <- 'output/margins'
measures <-  c("ap")
measure <-  c("ap")
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
dir.create(output_path, recursive = TRUE)

dat <- read_evaluation_data(measures, collections, 0.1)
dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
runs <- sample(names(dat[[measure]][[collection]]), 2)
runs_num <- c(gsub("[^[:digit:].]", "",  runs[1]), gsub("[^[:digit:].]", "",  runs[2])) # sampled runs in numeric

scores_s1 <- dat[[measure]][[collection]][[runs[[1]]]]
scores_s2 <- dat[[measure]][[collection]][[runs[[2]]]]
# split-half
T1_data_indices = sample(seq_len(length(scores_s1)), size = floor(0.5*length(scores_s1))) 
scores_s1_T1 = scores_s1[T1_data_indices]
scores_s2_T1 = scores_s2[T1_data_indices]
scores_s1_T2 = scores_s1[-T1_data_indices]
scores_s2_T2 = scores_s2[-T1_data_indices]

m_eff_A <- fit_all_margins(scores_s1, measure)[[1]]
m_eff_B <- fit_all_margins(scores_s2, measure)[[1]]
m_cdf_A <- function (X) { return(peff_2(X, m_eff_A)) }
m_cdf_B <- function (X) { return(peff_2(X, m_eff_B)) }
m_ecdf_A <- function (X) { return(sapply(X, function(x) sum(scores_s1 <= x) / length(scores_s1))) }
m_ecdf_B <- function (X) { return(sapply(X, function(x) sum(scores_s2 <= x) / length(scores_s2))) }
m_eff_A_T1 <- fit_all_margins(scores_s1_T1, measure)[[1]]
m_eff_B_T1 <- fit_all_margins(scores_s2_T1, measure)[[1]]
m_eff_A_T2 <- fit_all_margins(scores_s1_T2, measure)[[1]]
m_eff_B_T2 <- fit_all_margins(scores_s2_T2, measure)[[1]]
m_cdf_A_T1 <- function (X) { return(peff_2(X, m_eff_A_T1)) }
m_cdf_B_T1 <- function (X) { return(peff_2(X, m_eff_B_T1)) }
m_cdf_A_T2 <- function (X) { return(peff_2(X, m_eff_A_T2)) }
m_cdf_B_T2 <- function (X) { return(peff_2(X, m_eff_B_T2)) }
m_ecdf_A_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T1 <= x) / length(scores_s1_T1))) }
m_ecdf_B_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T1 <= x) / length(scores_s2_T1))) }
m_ecdf_A_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T2 <= x) / length(scores_s1_T2))) }
m_ecdf_B_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T2 <= x) / length(scores_s2_T2))) }

u <- pobs(cbind(scores_s1, scores_s2), ties.method = "random") # compute pseudo-observations
pseudoscores_s1 <- u[,1]
pseudoscores_s2 <- u[,2]

u1 <- pobs(cbind(scores_s1_T1, scores_s2_T1), ties.method = "random") # compute pseudo-observations
u2 <- pobs(cbind(scores_s1_T2, scores_s2_T2), ties.method = "random") # compute pseudo-observations
pseudoscores_s1_T1 <- u1[,1]
pseudoscores_s2_T1 <- u1[,2]
pseudoscores_s1_T2 <- u2[,1]
pseudoscores_s2_T2 <- u2[,2]

cop <- BiCopSelect(pseudoscores_s1, pseudoscores_s2, selectioncrit = "AIC") # fit copula
c_cdf <- function (x, y) { return(BiCopCDF_2(x, y, cop)) }

if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
  X <- seq(0, 1, length.out = 1001)
} else { # discrete case
  X <- support(measure)
}

# for CvM, KS
delta_observed <- cdf_cdf_diff(m_cdf_A_T1, m_ecdf_A_T2, measure, only_mean = FALSE)
delta_expected <- cdf_cdf_diff(m_ecdf_A_T1, m_ecdf_A_T2, measure, only_mean = FALSE)
KS_index <- order(delta_observed$highest_Y - delta_observed$lowest_Y,decreasing=T)[1]


emf(file=paste0(output_path, '/diag1_m_cdf_A__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_cdf_A(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag1_m_cdf_B__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_cdf_B(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag1_m_ecdf_A__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_ecdf_A(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag1_m_ecdf_B__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_ecdf_B(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0('output/copulas/diag1_copula_PDF_contour_', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.4,0.4,0.4,0.4))
plot(NA,NA,xlim = c(-3,3), ylim = c(-3,3), xaxt='n', yaxt='n', xlab='', ylab='',bty='L')
contour(cop, drawlabels = FALSE, col = 'darkgreen', add = TRUE, lwd = 1.25)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag2&3_m_cdf_A_T1__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_cdf_A_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag2&3_m_ecdf_A_T2__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_ecdf_A_T2(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='darkgreen', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag3_m_ecdf_A_T1__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_ecdf_A_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='purple', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag2_delta_KS_(m_cdf_A_T1_,_m_ecdf_A_T2)__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_cdf_A_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
lines(x=X, y=m_cdf_A_T1(X), col='blue', lty=1, lwd=2,)
lines(x=X, y=m_ecdf_A_T2(X), col='darkgreen', lty=1, lwd=2,)
arrows(x0=delta_observed$X[KS_index], y0=delta_observed$lowest_Y[KS_index], 
       x1=delta_observed$X[KS_index], y1=delta_observed$highest_Y[KS_index],
       col="red", lty=1, lwd=2, code=3, length=0.1)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag2&3_delta_cvm_(m_cdf_A_T1_,_m_ecdf_A_T2)__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_cdf_A_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
segments(x0=delta_observed$X, y0=delta_observed$lowest_Y, 
         x1=delta_observed$X, y1=delta_observed$highest_Y,
         col="red", lty=1, lwd=1)
lines(x=X, y=m_cdf_A_T1(X), col='blue', lty=1, lwd=2,)
lines(x=X, y=m_ecdf_A_T2(X), col='darkgreen', lty=1, lwd=2,)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag3_delta_cvm_(m_ecdf_A_T1_,_m_ecdf_A_T2)__', measure, '.emf'), width=1.8, height=1.8, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=m_ecdf_A_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='purple', frame.plot = FALSE)
segments(x0=delta_expected$X, y0=delta_expected$lowest_Y, x1=delta_expected$X, y1=delta_expected$highest_Y,
         col="red", lty=1, lwd=1)
lines(x=X, y=m_ecdf_A_T1(X), col='purple', lty=1, lwd=2,)
lines(x=X, y=m_ecdf_A_T2(X), col='darkgreen', lty=1, lwd=2,)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()



# ==================================================================================== 
# Plots for presenting KS, CvM in early Chapter 3
# ==================================================================================== 
set.seed(612368983)
output_path <- 'output/margins'
measures <-  c("ap")
measure <-  c("ap")
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
dir.create(output_path, recursive = TRUE)

dat <- read_evaluation_data(measures, collections, 0.1)
dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
runs <- sample(names(dat[[measure]][[collection]]), 2)
runs_num <- c(gsub("[^[:digit:].]", "",  runs[1]), gsub("[^[:digit:].]", "",  runs[2])) # sampled runs in numeric

scores_s1 <- dat[[measure]][[collection]][[runs[[1]]]]
scores_s2 <- dat[[measure]][[collection]][[runs[[2]]]]
# split-half
T1_data_indices = sample(seq_len(length(scores_s1)), size = floor(0.5*length(scores_s1))) 
scores_s1_T1 = scores_s1[T1_data_indices]
scores_s2_T1 = scores_s2[T1_data_indices]
scores_s1_T2 = scores_s1[-T1_data_indices]
scores_s2_T2 = scores_s2[-T1_data_indices]

m_eff_A <- fit_all_margins(scores_s1, measure)[[1]]
m_eff_B <- fit_all_margins(scores_s2, measure)[[1]]
m_cdf_A <- function (X) { return(peff_2(X, m_eff_A)) }
m_cdf_B <- function (X) { return(peff_2(X, m_eff_B)) }
m_ecdf_A <- function (X) { return(sapply(X, function(x) sum(scores_s1 <= x) / length(scores_s1))) }
m_ecdf_B <- function (X) { return(sapply(X, function(x) sum(scores_s2 <= x) / length(scores_s2))) }
m_eff_A_T1 <- fit_all_margins(scores_s1_T1, measure)[[1]]
m_eff_B_T1 <- fit_all_margins(scores_s2_T1, measure)[[1]]
m_eff_A_T2 <- fit_all_margins(scores_s1_T2, measure)[[1]]
m_eff_B_T2 <- fit_all_margins(scores_s2_T2, measure)[[1]]
m_cdf_A_T1 <- function (X) { return(peff_2(X, m_eff_A_T1)) }
m_cdf_B_T1 <- function (X) { return(peff_2(X, m_eff_B_T1)) }
m_cdf_A_T2 <- function (X) { return(peff_2(X, m_eff_A_T2)) }
m_cdf_B_T2 <- function (X) { return(peff_2(X, m_eff_B_T2)) }
m_ecdf_A_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T1 <= x) / length(scores_s1_T1))) }
m_ecdf_B_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T1 <= x) / length(scores_s2_T1))) }
m_ecdf_A_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T2 <= x) / length(scores_s1_T2))) }
m_ecdf_B_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T2 <= x) / length(scores_s2_T2))) }


if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
  X <- seq(0, 1, length.out = 1001)
} else { # discrete case
  X <- support(measure)
}

# for CvM, KS
delta_observed <- cdf_cdf_diff(m_cdf_A_T1, m_cdf_A_T2, measure, only_mean = FALSE)
delta_expected <- cdf_cdf_diff(m_ecdf_A_T1, m_cdf_A_T2, measure, only_mean = FALSE)
KS_index <- order(delta_observed$highest_Y - delta_observed$lowest_Y,decreasing=T)[1]


pdf(file=paste0(output_path, '/present_delta_KS_(m_cdf_A_T1_,m_cdf_A_T2)__', measure, '.pdf'), width=2.75, height=2.75, pointsize = 10)
par(mar=c(4.1,4.1,1.3,0.3))
plot(x=X, y=m_cdf_A_T1(X),
     xlab='X', ylab='Cumulative Probability',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
lines(x=X, y=m_cdf_A_T1(X), col='blue', lty=1, lwd=2,)
lines(x=X, y=m_cdf_A_T2(X), col='darkgreen', lty=1, lwd=2,)
arrows(x0=delta_observed$X[KS_index], y0=delta_observed$lowest_Y[KS_index], 
       x1=delta_observed$X[KS_index], y1=delta_observed$highest_Y[KS_index],
       col="red", lty=1, lwd=2, code=3, length=0.1)
axis(side=1, at=seq(0, 1, by=0.2))
axis(side=2, at=seq(0, 1, by=0.2))
box(bty='o')
title('KS')
legend(x='bottomright', legend=c(expression(F^{"*"}), expression(F)),
       col=c('blue', 'darkgreen'), lty=c(1,1), lwd=c(2,2))
dev.off()

pdf(file=paste0(output_path, '/present_delta_cvm_(m_cdf_A_T1_,m_cdf_A_T2)__', measure, '.pdf'), width=2.75, height=2.75, pointsize = 10)
par(mar=c(4.1,4.1,1.3,0.3))
plot(x=X, y=m_cdf_A_T1(X),
     xlab='X', ylab='Cumulative Probability',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
segments(x0=delta_observed$X, y0=delta_observed$lowest_Y, 
         x1=delta_observed$X, y1=delta_observed$highest_Y,
         col="red", lty=1, lwd=1)
lines(x=X, y=m_cdf_A_T1(X), col='blue', lty=1, lwd=2,)
lines(x=X, y=m_cdf_A_T2(X), col='darkgreen', lty=1, lwd=2,)
axis(side=1, at=seq(0, 1, by=0.2))
axis(side=2, at=seq(0, 1, by=0.2))
box(bty='o')
title('CvM')
legend(x='bottomright', legend=c(expression(F^{"*"}), expression(F)),
       col=c('blue', 'darkgreen'), lty=c(1,1), lwd=c(2,2))
dev.off()



# ==================================================================================== 
# Plots for presenting explaining the need for a Delta_expected in early Chapter 3
# ==================================================================================== 
set.seed(125114)
output_path <- 'output/margins'
measures <-  c("ap")
measure <-  c("ap")
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
dir.create(output_path, recursive = TRUE)

dat <- read_evaluation_data(measures, collections, 0.1)
dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
runs <- sample(names(dat[[measure]][[collection]]), 2)
runs_num <- c(gsub("[^[:digit:].]", "",  runs[1]), gsub("[^[:digit:].]", "",  runs[2])) # sampled runs in numeric

scores_s1 <- dat[[measure]][[collection]][[runs[[1]]]]
scores_s2 <- dat[[measure]][[collection]][[runs[[2]]]]
# split-half
T1_data_indices = sample(seq_len(length(scores_s1)), size = floor(0.5*length(scores_s1))) 
scores_s1_T1 = scores_s1[T1_data_indices]
scores_s2_T1 = scores_s2[T1_data_indices]
scores_s1_T2 = scores_s1[-T1_data_indices]
scores_s2_T2 = scores_s2[-T1_data_indices]

fit_all_margins(scores_s1_T1, measure)[[1]]$model$type
fit_all_margins(scores_s1_T1, measure)[[2]]$model$type
fit_all_margins(scores_s1_T1, measure)[[3]]$model$type
fit_all_margins(scores_s1_T1, measure)[[4]]$model$type

m_eff_A <- fit_all_margins(scores_s1, measure)[[1]]
m_eff_B <- fit_all_margins(scores_s2, measure)[[1]]
m_cdf_A <- function (X) { return(peff_2(X, m_eff_A)) }
m_cdf_B <- function (X) { return(peff_2(X, m_eff_B)) }
m_ecdf_A <- function (X) { return(sapply(X, function(x) sum(scores_s1 <= x) / length(scores_s1))) }
m_ecdf_B <- function (X) { return(sapply(X, function(x) sum(scores_s2 <= x) / length(scores_s2))) }
m_eff_A_T1 <- fit_all_margins(scores_s1_T1, measure)[[1]]
m_eff_A_T1_top2 <- fit_all_margins(scores_s1_T1, measure)[[2]]
m_eff_B_T1 <- fit_all_margins(scores_s2_T1, measure)[[1]]
m_eff_A_T2 <- fit_all_margins(scores_s1_T2, measure)[[1]]
m_eff_B_T2 <- fit_all_margins(scores_s2_T2, measure)[[1]]
m_cdf_A_T1 <- function (X) { return(peff_2(X, m_eff_A_T1)) }
m_cdf_A_T1_top2 <- function (X) { return(peff_2(X, m_eff_A_T1_top2)) }
m_cdf_B_T1 <- function (X) { return(peff_2(X, m_eff_B_T1)) }
m_cdf_A_T2 <- function (X) { return(peff_2(X, m_eff_A_T2)) }
m_cdf_B_T2 <- function (X) { return(peff_2(X, m_eff_B_T2)) }
m_ecdf_A_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T1 <= x) / length(scores_s1_T1))) }
m_ecdf_B_T1 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T1 <= x) / length(scores_s2_T1))) }
m_ecdf_A_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s1_T2 <= x) / length(scores_s1_T2))) }
m_ecdf_B_T2 <- function (X) { return(sapply(X, function(x) sum(scores_s2_T2 <= x) / length(scores_s2_T2))) }

if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
  X <- seq(0, 1, length.out = 1001)
} else { # discrete case
  X <- support(measure)
}

# for CvM, KS
delta_observed <- cdf_cdf_diff(m_cdf_A_T1, m_ecdf_A_T2, measure, only_mean = FALSE)
delta_observed_top2 <- cdf_cdf_diff(m_cdf_A_T1_top2, m_ecdf_A_T2, measure, only_mean = FALSE)
delta_expected <- cdf_cdf_diff(m_ecdf_A_T1, m_ecdf_A_T2, measure, only_mean = FALSE)
plot_data <- delta_observed


pdf(file=paste0(output_path, '/need_for_Dexp_(', measure, ')1.pdf'), width=3, height=3, pointsize = 10)
par(mar=c(4.2,4.2,1.9,0.4))
plot(x=plot_data$X, y=plot_data$cdf_1_Y,
     xlab='X', ylab='Cumulative Probability',
     type='l', xlim=c(0,1), ylim=c(0,1))
title(bquote(bold(paste(Delta[obs], ' = ', .(paste0(round(plot_data$mean, 3)))))))
axis(side=1, at=seq(0, 1, by=0.2))
axis(side=2, at=seq(0, 1, by=0.2))
segments(x0=plot_data$X, y0=plot_data$lowest_Y, x1=plot_data$X, y1=plot_data$highest_Y,
         col="red", lty=1, lwd=2)
lines(x=plot_data$X, y=plot_data$cdf_2_Y, col='darkgreen', lty=1, lwd=2,)
lines(x=plot_data$X, y=plot_data$cdf_1_Y, col='blue', lty=1, lwd=2)
legend(x='bottomright', legend=c(bquote(paste(F[1]^{"*"}, ' (', .(beautify(m_eff_A_T1$model$type)), ')')),
                                 expression(F[2])),
       col=c('blue','darkgreen'), lty=c(1,1), lwd=c(2,2))
dev.off()

pdf(file=paste0(output_path, '/need_for_Dexp_(', measure, ')2.pdf'), width=3, height=3, pointsize = 10)
par(mar=c(4.2,4.2,1.9,0.4))
plot(x=plot_data$X, y=delta_observed_top2$cdf_1_Y,
     xlab='X', ylab='Cumulative Probability',
     type='l', xlim=c(0,1), ylim=c(0,1))
title(bquote(bold(paste(Delta[obs], ' = ', .(paste0(round(delta_observed_top2$mean, 3)))))))
axis(side=1, at=seq(0, 1, by=0.2))
axis(side=2, at=seq(0, 1, by=0.2))
segments(x0=plot_data$X, y0=delta_observed_top2$lowest_Y, x1=delta_observed_top2$X, y1=delta_observed_top2$highest_Y,
         col="red", lty=1, lwd=2)
lines(x=plot_data$X, y=plot_data$cdf_2_Y, col='darkgreen', lty=1, lwd=2,)
lines(x=plot_data$X, y=delta_observed_top2$cdf_1_Y, col='blue', lty=1, lwd=2)
legend(x='bottomright', legend=c(bquote(paste(F[1]^{"*"}, ' (', .(beautify(m_eff_A_T1_top2$model$type)), ')')),
                                 expression(F[2])),
       col=c('blue','darkgreen'), lty=c(1,1), lwd=c(2,2))
dev.off()

pdf(file=paste0(output_path, '/need_for_Dexp_(', measure, ')3.pdf'), width=3, height=3, pointsize = 10)
par(mar=c(4.2,4.2,1.9,0.4))
plot_data <- delta_expected
plot(x=plot_data$X, y=plot_data$cdf_1_Y,
     xlab='X', ylab='Cumulative Probability',
     type='l', xlim=c(0,1), ylim=c(0,1))
title(bquote(bold(paste(Delta[exp], ' = ', .(paste0(round(plot_data$mean, 3)))))))
axis(side=1, at=seq(0, 1, by=0.2))
axis(side=2, at=seq(0, 1, by=0.2))
segments(x0=plot_data$X, y0=plot_data$lowest_Y, x1=plot_data$X, y1=plot_data$highest_Y,
         col="red", lty=1, lwd=2)
lines(x=plot_data$X, y=plot_data$cdf_2_Y, col='darkgreen', lty=1, lwd=2,)
lines(x=plot_data$X, y=plot_data$cdf_1_Y, col='blue', lty=1, lwd=2)
legend(x='bottomright', legend=c(bquote(paste(F[1])),
                                 expression(F[2])),
       col=c('blue','darkgreen'), lty=c(1,1), lwd=c(2,2))
dev.off()



# ==============================================================================
# Plots for Diagram 4
# ==============================================================================
source("R/utils.R")
library(plotly)
set.seed(8235238)
output_path <- 'output/margins'
dir.create(output_path, recursive = TRUE)
measures <-  c("ap")
measure <-  c("ap")
collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
  X <- seq(0, 1, length.out = 1001)
} else { # discrete case
  X <- support(measure)
}

dat <- read_evaluation_data(measures, collections, 0.1)
dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
run <- sample(names(dat[[measure]][[collection]]), 1)
runs_num <- gsub("[^[:digit:].]", "",  run) # sampled runs in numeric
scores <- dat[[measure]][[collection]][[run]]
scores <- sample(scores, 25)

for (i in c(1,2)) {
  # split-half
  T1_data_indices = sample(seq_len(length(scores)), size = floor(0.5*length(scores))) 
  scores_T1 = scores[T1_data_indices]
  scores_T2 = scores[-T1_data_indices]
  
  length(scores_T1)
  length(scores_T2)
  
  eff_T1 <- fit_all_margins(scores_T1, measure)
  cdf_T1_top1 <- function (X) { return(peff_2(X, eff_T1[[1]])) }
  cdf_T1_top2 <- function (X) { return(peff_2(X, eff_T1[[2]])) }
  cdf_T1_top3 <- function (X) { return(peff_2(X, eff_T1[[3]])) }
  cdf_T1_top4 <- function (X) { return(peff_2(X, eff_T1[[4]])) }
  ecdf_T2 <- function (X) { return(sapply(X, function(x) sum(scores_T2 <= x) / length(scores_T2))) }
  
  
  print(paste(beautify(eff_T1[[1]]$model$type), cdf_cdf_diff(cdf_T1_top1, ecdf_T2, measure, only_mean = T)))
  print(paste(beautify(eff_T1[[2]]$model$type), cdf_cdf_diff(cdf_T1_top2, ecdf_T2, measure, only_mean = T)))
  print(paste(beautify(eff_T1[[3]]$model$type), cdf_cdf_diff(cdf_T1_top3, ecdf_T2, measure, only_mean = T)))
  print(paste(beautify(eff_T1[[4]]$model$type), cdf_cdf_diff(cdf_T1_top4, ecdf_T2, measure, only_mean = T)))
  
  emf(file=paste0(output_path, '/diag4_cdf_1__', measure, '_', beautify(eff_T1[[1]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(x=X, y=cdf_T1_top1(X),
       xlab='', ylab='',
       type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
  text(x=0.7, y=0.075, col='black',
       labels=bquote(paste(Delta[obs], ' = ', .(paste0(round(cdf_cdf_diff(cdf_T1_top1, ecdf_T2, measure, only_mean = T), 2))))))
  text(x=0.7, y=0.25, col='black',
       labels=beautify(eff_T1[[1]]$model$type))
  axis(1, lwd.tick=0)
  axis(1, lwd.tick=0)
  axis(2, lwd.tick=0)
  box(bty='L')
  dev.off()
  
  emf(file=paste0(output_path, '/diag4_cdf_2__', measure, '_', beautify(eff_T1[[2]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(x=X, y=cdf_T1_top2(X),
       xlab='', ylab='',
       type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
  text(x=0.7, y=0.075, col='black',
       labels=bquote(paste(Delta[obs], ' = ', .(paste0(round(cdf_cdf_diff(cdf_T1_top2, ecdf_T2, measure, only_mean = T), 2))))))
  text(x=0.7, y=0.25, col='black',
       labels=beautify(eff_T1[[2]]$model$type))
  axis(1, lwd.tick=0)
  axis(1, lwd.tick=0)
  axis(2, lwd.tick=0)
  box(bty='L')
  dev.off()
  
  emf(file=paste0(output_path, '/diag4_cdf_3__', measure, '_', beautify(eff_T1[[3]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(x=X, y=cdf_T1_top3(X),
       xlab='', ylab='',
       type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
  text(x=0.7, y=0.075, col='black',
       labels=bquote(paste(Delta[obs], ' = ', .(paste0(round(cdf_cdf_diff(cdf_T1_top3, ecdf_T2, measure, only_mean = T), 2))))))
  text(x=0.7, y=0.25, col='black',
       labels=beautify(eff_T1[[3]]$model$type))
  axis(1, lwd.tick=0)
  axis(1, lwd.tick=0)
  axis(2, lwd.tick=0)
  box(bty='L')
  dev.off()
  
  emf(file=paste0(output_path, '/diag4_cdf_4__', measure, '_', beautify(eff_T1[[4]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(x=X, y=cdf_T1_top4(X),
       xlab='', ylab='',
       type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
  text(x=0.7, y=0.075, col='black',
       labels=bquote(paste(Delta[obs], ' = ', .(paste0(round(cdf_cdf_diff(cdf_T1_top4, ecdf_T2, measure, only_mean = T), 2))))))
  text(x=0.7, y=0.25, col='black',
       labels=beautify(eff_T1[[4]]$model$type))
  axis(1, lwd.tick=0)
  axis(1, lwd.tick=0)
  axis(2, lwd.tick=0)
  box(bty='L')
  dev.off()
  
  emf(file=paste0(output_path, '/diag4_ecdf_4__', measure, '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
  par(mar=c(0.3,0.3,0.3,0.3))
  plot(x=X, y=ecdf_T2(X),
       xlab='', ylab='',
       type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='darkgreen', frame.plot = FALSE)
  text(x=0.7, y=0.075, col='black',
       labels='Empirical')
  axis(1, lwd.tick=0)
  axis(1, lwd.tick=0)
  axis(2, lwd.tick=0)
  box(bty='L')
  dev.off()
}


eff_T1 <- fit_all_margins(scores, measure, N_TRIALS = 10, compute_splithalf_criterion = T)
cdf_T1_top1 <- function (X) { return(peff_2(X, eff_T1[[1]])) }
cdf_T1_top2 <- function (X) { return(peff_2(X, eff_T1[[2]])) }
cdf_T1_top3 <- function (X) { return(peff_2(X, eff_T1[[3]])) }
cdf_T1_top4 <- function (X) { return(peff_2(X, eff_T1[[4]])) }
ecdf_T2 <- function (X) { return(sapply(X, function(x) sum(scores_T2 <= x) / length(scores_T2))) }


print(paste(beautify(eff_T1[[1]]$model$type), cdf_cdf_diff(cdf_T1_top1, ecdf_T2, measure, only_mean = T)))
print(paste(beautify(eff_T1[[2]]$model$type), cdf_cdf_diff(cdf_T1_top2, ecdf_T2, measure, only_mean = T)))
print(paste(beautify(eff_T1[[3]]$model$type), cdf_cdf_diff(cdf_T1_top3, ecdf_T2, measure, only_mean = T)))
print(paste(beautify(eff_T1[[4]]$model$type), cdf_cdf_diff(cdf_T1_top4, ecdf_T2, measure, only_mean = T)))
i<-3

emf(file=paste0(output_path, '/diag4_cdf_1__', measure, '_', beautify(eff_T1[[1]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1_top1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
text(x=0.7, y=0.075, col='black',
     labels=bquote(paste(SHC, ' = ', .(paste0(round(eff_T1[[1]]$splithalf_criterion, 2))))))
text(x=0.7, y=0.25, col='black',
     labels=beautify(eff_T1[[1]]$model$type))
axis(1, lwd.tick=0)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag4_cdf_2__', measure, '_', beautify(eff_T1[[2]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1_top2(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
text(x=0.7, y=0.075, col='black',
     labels=bquote(paste(SHC, ' = ', .(paste0(round(eff_T1[[2]]$splithalf_criterion, 2))))))
text(x=0.7, y=0.25, col='black',
     labels=beautify(eff_T1[[2]]$model$type))
axis(1, lwd.tick=0)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag4_cdf_3__', measure, '_', beautify(eff_T1[[3]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1_top3(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
text(x=0.7, y=0.075, col='black',
     labels=bquote(paste(SHC, ' = ', .(paste0(round(eff_T1[[3]]$splithalf_criterion, 2))))))
text(x=0.7, y=0.25, col='black',
     labels=beautify(eff_T1[[3]]$model$type))
axis(1, lwd.tick=0)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag4_cdf_4__', measure, '_', beautify(eff_T1[[4]]$model$type), '_', i, '.emf'), width=1.2, height=1.2, pointsize = 10)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1_top4(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
text(x=0.7, y=0.075, col='black',
     labels=bquote(paste(SHC, ' = ', .(paste0(round(eff_T1[[4]]$splithalf_criterion, 2))))))
text(x=0.7, y=0.25, col='black',
     labels=beautify(eff_T1[[4]]$model$type))
axis(1, lwd.tick=0)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()



# ==============================================================================
# Plots for Diagram 5
# ==============================================================================
source("R/utils.R")
set.seed(6235) # 6235
output_path <- 'output/margins'
dir.create(output_path, recursive = TRUE)
measures <-  c("ap", 'p10', 'rr')
measure <-  c("ap")
collections <- c("terabyte2006")
collection <- "terabyte2006"
split_n = c(25,50)
if(measure %in% c("ap", "ndcg20", "err20")) { # continuous case
  X <- seq(0, 1, length.out = 1001)
} else { # discrete case
  X <- support(measure)
}

dat <- read_evaluation_data(measures, collections, 0.1)
dat_c <- sapply(dat[[measure]], length) # count of runs per collection, for sampling
run <- sample(names(dat[[measure]][[collection]]), 1)
runs_num <- gsub("[^[:digit:].]", "",  run) # sampled runs in numeric
d <- dat[[measure]][[collection]][[run]]

# Split n1-n, n2-n (i.e., 25-99 and 50-99)
T2_data_indices <- sample(seq_len(length(d)), size = split_n[2])
T1_data_indices <- sample(T2_data_indices, size = split_n[1])
T2_data <- d[T2_data_indices]
T_data <- d[-T2_data_indices]
T1_data <- d[T1_data_indices]

eff_T1 <- fit_all_margins(T1_data, measure)[[1]]
eff_T2 <- fit_all_margins(T2_data, measure)[[1]]
eff_T <- fit_all_margins(T_data, measure)[[1]]

cdf_T1 <- function (X) { return(peff_2(X, eff_T1)) }
cdf_T2 <- function (X) { return(peff_2(X, eff_T2)) }
cdf_T <- function (X) { return(peff_2(X, eff_T)) }
ecdf_T1 <- function (X) { return(sapply(X, function(x) sum(T1_data <= x) / length(T1_data))) }
ecdf_T2 <- function (X) { return(sapply(X, function(x) sum(T2_data <= x) / length(T2_data))) }
ecdf_T <- function (X) { return(sapply(X, function(x) sum(T_data <= x) / length(T_data))) }

delta_obs_t1 <- cdf_cdf_diff(cdf_T1, ecdf_T, measure, only_mean = FALSE)
delta_obs_t2 <- cdf_cdf_diff(cdf_T2, ecdf_T, measure, only_mean = FALSE)
delta_exp_t1 <- cdf_cdf_diff(ecdf_T1, ecdf_T, measure, only_mean = FALSE)

size=1.4

emf(file=paste0(output_path, '/diag5_cdf_T1__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='blue', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_cdf_T2__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T2(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='red', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_ecdf_T1__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=ecdf_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='purple', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_ecdf_T__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=ecdf_T(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", lty=1, lwd=2, col='darkgreen', frame.plot = FALSE)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_delta_cvm_(cdf_T1,ecdf_T)__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
segments(x0=delta_obs_t1$X, y0=delta_obs_t1$lowest_Y, 
         x1=delta_obs_t1$X, y1=delta_obs_t1$highest_Y,
         col="gray", lty=1, lwd=1)
lines(x=X, y=cdf_T1(X), col='blue', lty=1, lwd=2,)
lines(x=X, y=ecdf_T(X), col='darkgreen', lty=1, lwd=2,)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_delta_cvm_(cdf_T2,ecdf_T)__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=cdf_T2(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
segments(x0=delta_obs_t2$X, y0=delta_obs_t2$lowest_Y, 
         x1=delta_obs_t2$X, y1=delta_obs_t2$highest_Y,
         col="gray", lty=1, lwd=1)
lines(x=X, y=cdf_T2(X), col='red', lty=1, lwd=2,)
lines(x=X, y=ecdf_T(X), col='darkgreen', lty=1, lwd=2,)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()

emf(file=paste0(output_path, '/diag5_delta_cvm_(ecdf_T1,ecdf_T)__', measure, '.emf'), width=size, height=size, pointsize = 9)
par(mar=c(0.3,0.3,0.3,0.3))
plot(x=X, y=ecdf_T1(X),
     xlab='', ylab='',
     type='l', xlim=c(0,1), ylim=c(0,1), xaxt = "n", yaxt = "n", 
     lty=1, lwd=2, col='blue', frame.plot = FALSE)
segments(x0=delta_exp_t1$X, y0=delta_exp_t1$lowest_Y, 
         x1=delta_exp_t1$X, y1=delta_exp_t1$highest_Y,
         col="gray", lty=1, lwd=1)
lines(x=X, y=ecdf_T1(X), col='purple', lty=1, lwd=2,)
lines(x=X, y=ecdf_T(X), col='darkgreen', lty=1, lwd=2,)
axis(1, lwd.tick=0)
axis(2, lwd.tick=0)
box(bty='L')
dev.off()
