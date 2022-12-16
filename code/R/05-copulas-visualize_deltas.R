source("R/utils.R")


# ============================= Plot for diag 6 =============================
set.seed(5124)
path_output <- 'output/copulas'
dir.create(path_output, recursive = TRUE)

# measures <-  c("ap", "p10", "rr", "ndcg20", "err20")
# collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
measures <-  c("ap")
collections <- c("adhoc7")
measure <- sample(measures, 1)

# Read evaluation data
dat <- read_evaluation_data(measures, collections, 0.1)
# Count the number of runs per collection, for sampling
dat_c <- sapply(dat[[measure]], length) 
# Sample collection, proportional to number of runs in it
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
# Sample 2 runs
runs <- sample(names(dat[[measure]][[collection]]), 2)
runs_num <- c(gsub("[^[:digit:].]", "",  runs[1]), gsub("[^[:digit:].]", "",  runs[2])) # sampled runs in numeric
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

# Fit copula model
cop_t1 <- BiCopSelect(pseudoscores_s1_T1, pseudoscores_s2_T1, selectioncrit = "AIC") # fit copula
cop_t2 <- BiCopSelect(pseudoscores_s1_T2, pseudoscores_s2_T2, selectioncrit = "AIC") # fit copula

# Compute/define CDFs/ECDFs
cdf_t1 <- function (x, y) { return(BiCopCDF_2(x, y, cop_t1)) }
cdf_t2 <- function (x, y) { return(BiCopCDF_2(x, y, cop_t2)) }
ecdf_t1 <- function (x, y) {return(sapply(1:length(x), function(i) mean(pseudoscores_s1_T1 <= x[i] & pseudoscores_s2_T1 <= y[i])))}
ecdf_t2 <- function (x, y) {return(sapply(1:length(x), function(i) mean(pseudoscores_s1_T2 <= x[i] & pseudoscores_s2_T2 <= y[i])))}

# All copulas are continuous, so we take 100 equally spaced values within the [0,1] range as x-range and y-range
n <- 100
xy_range <- seq(0, 1, length.out = n)

# Z1 is a matrix such that i.e, Z1[x0, y0] == cdf_t1(x0, y0)
Z1 <- sapply(xy_range, function(x) cdf_t1(xy_range, rep(x, n)))
Z11 <- sapply(xy_range, function(x) cdf_t2(xy_range, rep(x, n)))
Z2 <- sapply(xy_range, function(x) ecdf_t1(xy_range, rep(x, n)))
Z22 <- sapply(xy_range, function(x) ecdf_t2(xy_range, rep(x, n)))
cdf_t1_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z1)))
cdf_t2_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z11)))
ecdf_t1_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z2)))
ecdf_t2_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z22)))
diff1 <- joint_cdf_cdf_diff(cdf_t1, ecdf_t2)
diff2 <- joint_cdf_cdf_diff(ecdf_t1, ecdf_t2)


p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=cdf_t1_data$x, y=cdf_t1_data$y, z=cdf_t1_data$z, intensity =cdf_t1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("blue","blue")), name=TeX(paste('CDF_{T1}\\text{(', beautify(cop_t1$familyname), ')}'))) %>%
  layout(font=list(size = 14),
         margin=list(l=0, r=0, t=0, b=0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Z'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0))
print(p)


p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=ecdf_t1_data$x, y=ecdf_t1_data$y, z=ecdf_t1_data$z, intensity =ecdf_t1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("purple","purple")), name=TeX(paste('ECDF_{T1}'))) %>%
  layout(font=list(size = 14),
         margin=list(l=0, r=0, t=0, b=0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Z'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0))
print(p)


p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=ecdf_t2_data$x, y=ecdf_t2_data$y, z=ecdf_t2_data$z, intensity =ecdf_t2_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("darkgreen","darkgreen")), name=TeX(paste('ECDF_{T2}'))) %>%
  layout(font=list(size = 14),
         margin=list(l=0, r=0, t=0, b=0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Z'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0))
print(p)


p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=cdf_t1_data$x, y=cdf_t1_data$y, z=cdf_t1_data$z, intensity =cdf_t1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("blue","blue")), name=TeX(paste('CDF_{T1}\\text{(', beautify(cop_t1$familyname), ')}'))) %>%
  add_trace(x=ecdf_t2_data$x, y=ecdf_t2_data$y, z=ecdf_t2_data$z, intensity =ecdf_t2_data$z, type="mesh3d",  opacity=1,  colorscale = list(c(0,1),c("darkgreen","darkgreen")), name=TeX('ECDF_{T2}')) %>%
  layout(title=list(text = TeX(paste0('\\Large{\\frac{\\text{ Example: {', beautify(measure), ', run', runs_num[1], '&', runs_num[2], ', ', collection, '}}}',
                                      '{\\text{  }\\Delta_{1}(CDF_{T1}, ECDF_{T2}) = ', round(diff1, 4), '}}')),
                    y=0.94), # move title down
         font=list(size = 14),
         margin=list(l=0, r=0, t=0, b=0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Z'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0))
print(p)


p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=ecdf_t1_data$x, y=ecdf_t1_data$y, z=ecdf_t1_data$z, intensity =ecdf_t1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("purple","purple")), name=TeX('ECDF_{T1}')) %>%
  add_trace(x=ecdf_t2_data$x, y=ecdf_t2_data$y, z=ecdf_t2_data$z, intensity =ecdf_t2_data$z, type="mesh3d",  opacity=1,  colorscale = list(c(0,1),c("darkgreen","darkgreen")), name=TeX('ECDF_{T2}')) %>%
  layout(title=list(text = TeX(paste0('\\Large{\\frac{\\text{ Example: {', beautify(measure), ', run', runs_num[1], '&', runs_num[2], ', ', collection, '}}}',
                                      '{\\text{  }\\Delta_{2}(ECDF_{T1}, ECDF_{T2}) = ', round(diff2, 4), '}}')),
                    y=0.94), # move title down
         font=list(size = 14),
         margin=list(t = 0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Z'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0))
print(p)



# ============================= Plot for showing an example random split in chapter 4 =============================
source("R/utils.R")
set.seed(345345)
path_output <- 'output/copulas'
dir.create(path_output, recursive = TRUE)

# measures <-  c("ap", "p10", "rr", "ndcg20", "err20")
# collections <- c("adhoc5", "adhoc6", "adhoc7", "adhoc8", "web2010", "web2011", "web2012", "web2013")
measures <-  c("ap")
collections <- c("adhoc7")
measure <- sample(measures, 1)

# Read evaluation data
dat <- read_evaluation_data(measures, collections, 0.1)
# Count the number of runs per collection, for sampling
dat_c <- sapply(dat[[measure]], length) 
# Sample collection, proportional to number of runs in it
collection <- sample(names(dat[[measure]]), 1, prob = dat_c / sum(dat_c))
# Sample 2 runs
runs <- sample(names(dat[[measure]][[collection]]), 2)
runs_num <- c(gsub("[^[:digit:].]", "",  runs[1]), gsub("[^[:digit:].]", "",  runs[2])) # sampled runs in numeric
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

# Fit copula model
cops <- fit_all_copulas(pseudoscores_s1_T1, pseudoscores_s2_T1)
cop_1 <- cops[[1]]
cop_2 <- cops[[2]]

# Compute/define CDFs/ECDFs
cdf_1 <- function (x, y) { return(BiCopCDF_2(x, y, cop_1)) }
cdf_2 <- function (x, y) { return(BiCopCDF_2(x, y, cop_2)) }
ecdf_t1 <- function (x, y) {return(sapply(1:length(x), function(i) mean(pseudoscores_s1_T1 <= x[i] & pseudoscores_s2_T1 <= y[i])))}
ecdf_t2 <- function (x, y) {return(sapply(1:length(x), function(i) mean(pseudoscores_s1_T2 <= x[i] & pseudoscores_s2_T2 <= y[i])))}

# All copulas are continuous, so we take 100 equally spaced values within the [0,1] range as x-range and y-range
n <- 100
xy_range <- seq(0, 1, length.out = n)

# Z1 is a matrix such that i.e, Z1[x0, y0] == cdf_t1(x0, y0)
Z1 <- sapply(xy_range, function(x) cdf_1(xy_range, rep(x, n)))
Z11 <- sapply(xy_range, function(x) cdf_2(xy_range, rep(x, n)))
Z2 <- sapply(xy_range, function(x) ecdf_t1(xy_range, rep(x, n)))
Z22 <- sapply(xy_range, function(x) ecdf_t2(xy_range, rep(x, n)))
cdf_1_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z1)))
cdf_2_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z11)))
ecdf_t1_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z2)))
ecdf_t2_data <- data.frame(x=rep(xy_range, n), y=rep(xy_range, each=n), z=unlist(as.list(Z22)))
diff1 <- joint_cdf_cdf_diff(cdf_1, ecdf_t2)
diff2 <- joint_cdf_cdf_diff(cdf_2, ecdf_t2)
diff3 <- joint_cdf_cdf_diff(ecdf_t1, ecdf_t2)


# Print 
beautify(measure)
runs_num[1]
runs_num[2]
collection


# Cop 1
p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=cdf_1_data$x, y=cdf_1_data$y, z=cdf_1_data$z, intensity =cdf_1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("blue","blue")), name=TeX(paste0('F_{1}^{*}\\text{ (', beautify(cop_1$familyname), ')}'))) %>%
  add_trace(x=ecdf_t2_data$x, y=ecdf_t2_data$y, z=ecdf_t2_data$z, intensity =ecdf_t2_data$z, type="mesh3d",  opacity=1,  colorscale = list(c(0,1),c("darkgreen","darkgreen")), name=TeX('F_{2}')) %>%
  layout(title=list(font=list(size = 29), text = TeX(paste0('\\bf{\\Delta_{obs} = ', round(diff1, 3), '}')),
                    y=0.94), # move title down
         font=list(size = 14),
         margin=list(l=0, r=0, t=0, b=0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Joint Cumulative Probability'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0, font=list(size = 28)))
print(p)


# 2 empiricals
p <- plot_ly(showscale = FALSE, showlegend=TRUE, width=650, height=550) %>%
  config(mathjax = 'cdn')  %>%
  add_trace(x=ecdf_t1_data$x, y=ecdf_t1_data$y, z=ecdf_t1_data$z, intensity =ecdf_t1_data$z, type="mesh3d", opacity=1,  colorscale = list(c(0,1),c("purple","purple")), name=TeX('F_{1}')) %>%
  add_trace(x=ecdf_t2_data$x, y=ecdf_t2_data$y, z=ecdf_t2_data$z, intensity =ecdf_t2_data$z, type="mesh3d",  opacity=1,  colorscale = list(c(0,1),c("darkgreen","darkgreen")), name=TeX('F_{2}')) %>%
  layout(title=list(font=list(size = 29), text = TeX(paste0('\\bf{\\Delta_{exp} = ', round(diff3, 3), '}')),
                    y=0.94), # move title down
         font=list(size = 14),
         margin=list(t = 0), # add top margin for title if needed
         scene = list(xaxis=list(title='X'),
                      yaxis=list(title='Y'),
                      zaxis=list(title='Joint Cumulative Probability'),
                      camera=list(eye = list(x = 1.25, y = -1.25, z = 0.25))),
         legend = list(orientation = "h", xanchor = "center", x = 0.5, y=0, font=list(size = 28)))
print(p)
