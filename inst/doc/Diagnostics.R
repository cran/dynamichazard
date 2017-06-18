## ----setup, include=FALSE------------------------------------------------
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "... output abbreviated ..."
   if(is.list(lines)){
      x <- lapply(lines, function(z) x[z])
      for(i in 1:(length(x) - 1))
        x[[i]] <- c(x[[i]], more)
      x <- unlist(x)
   }
   else if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(more, x[lines], more)
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(echo = TRUE, warning = F, message = F, 
                      dpi = 144, cache = T, dev = "png")
knitr::opts_knit$set(warning = F, message = F)
options(digits = 3)

knit_print.data.frame <- function(
  x, ..., row.names = F, digits = getOption("digits")){
  res <- paste(c(
    "", "", knitr::kable(x, ..., row.names = row.names, digits = digits)), 
    collapse = "\n")
  knitr::asis_output(res)
}
knit_print.matrix <- function(x, ..., digits = getOption("digits"), row.names = F, escape = T){
  col.names <- rep("", ncol(x))
  col.names <- if(is.null(colnames(x)))
    rep("", ncol(x)) else colnames(x)
  # The long table does not do well with the ( or [ symbol 
  # See http://tex.stackexchange.com/questions/228755/undefined-control-sequence-on-left-parenthesis-after-midrule-in-longtable
  col.names[1] <- stringr::str_replace(col.names[1], "^\\(", "\\\\relax(")
  col.names[1] <- stringr::str_replace(col.names[1], "^\\[", "\\\\relax[")
  res <- paste(c(
    "", "", knitr::kable(x, ..., row.names = row.names, col.names = col.names, 
                         escape = escape, digits = digits)),
    collapse = "\n")
  knitr::asis_output(res)
}

## ------------------------------------------------------------------------
load("Diagnostics/Rossi.RData")

# We only keep some of the columns
Rossi <- Rossi[
  , c("id","start","stop","event", "fin", "age", "prio", "employed.cumsum")]

## ------------------------------------------------------------------------
# Number of unique individauls
length(unique(Rossi$id)) 

# Show example for a given individual
Rossi[Rossi$id == 2, ]

## ------------------------------------------------------------------------
# See the varying and non-varying covariates
# The output shows the mean number of unique values for each individual
tmp <- 
  by(Rossi[, ], Rossi$id, function(dat) 
    apply(dat, 2, function(x) sum(!duplicated(x))))
colMeans(do.call(rbind, as.list(tmp)))

## ---- echo = FALSE-------------------------------------------------------
# Indivuals are observed for at most one year
stopifnot(setequal(unique(Rossi$stop), 1:52))

## ---- echo = FALSE-------------------------------------------------------
# Individuals have event through out the period
plot(xtabs(~ Rossi$stop[Rossi$event == 1]), xlab = "Week number", 
     ylab = "Number of recidivism")

## ---- echo = FALSE-------------------------------------------------------
# All individuals have gabs of one week 
stopifnot(all(
  unlist(tapply(Rossi$stop, Rossi$id, diff)) == 1))

# People have at most one event
is_event <- Rossi[Rossi$id %in% Rossi$id[Rossi$event == 1], ]
stopifnot(all(tapply(is_event$event, is_event$id, sum) == 1))

# People alwas get arrested on the last record
stopifnot(all(
  tapply(1:nrow(is_event), is_event$id, function(is){
    rows <- is_event[is, ]
    max(rows$stop) == rows$stop[rows$event == 1]
  })))

## ------------------------------------------------------------------------
library(dynamichazard)
dd_rossi <- ddhazard(
  Surv(start, stop, event) ~ fin + age + prio + employed.cumsum, 
  data = Rossi, id = Rossi$id, by = 1, max_T = 52, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

## ------------------------------------------------------------------------
plot(dd_rossi)

## ------------------------------------------------------------------------
hats <- hatvalues(dd_rossi)

## ------------------------------------------------------------------------
str(hats[1:3]) # Print str of first three matrices
head(hats[[1]], 10) # Print the head of first matrix

## ----stack_hats, echo = FALSE--------------------------------------------
stack_hats <- function(hats){
  # Stack
  resids_hats <- data.frame(do.call(rbind, hats), row.names = NULL)
  
  # Add the interval number
  n_rows <- unlist(lapply(hats, nrow))
  interval_n <- unlist(sapply(1:length(n_rows), function(i) rep(i, n_rows[i])))
  
  resids_hats <- cbind(resids_hats, interval_n = interval_n)
  
  # Sort by id and interval number
  resids_hats <- 
    resids_hats[order(resids_hats$id, resids_hats$interval_n), ]
  
  resids_hats
}

## ------------------------------------------------------------------------
hats_stacked <- stack_hats(hats)

head(hats_stacked)

## ------------------------------------------------------------------------
# Draw histogram of hat values
hist(log10(hats_stacked$hat_value), main = "",
     xlab = "Histogram of log10 of Hat values")

## ------------------------------------------------------------------------
# Print the largest values 
max_hat <- tapply(hats_stacked$hat_value, hats_stacked$id, max)
head(sort(max_hat, decreasing = T), 5)

## ----rossi_hat_plot------------------------------------------------------
# We will highlight the individuals with the highest hatvalues
is_large <- 
  names(head(sort(max_hat, decreasing = T), 5))

# Plot hat values
plot(range(hats_stacked$interval_n), c(0, 0.22), type = "n",
     xlab = "Interval number", ylab = "Hat value")

invisible(
  by(hats_stacked, hats_stacked$id, function(rows){
    has_large <- rows$id[1] %in% is_large
    col <- if(has_large) "Red" else "Black"
    lines(rows$interval_n, rows$hat_value, lty = 2, 
          col = col)
    
    if(has_large){
      pch <- as.character(which(rows$id[1] == is_large))
      points(rows$interval_n, rows$hat_value, pch = pch, col = col)
    }
  }))

## ------------------------------------------------------------------------
# These are the individuals id
is_large

# We print the last record each of these
Rossi_subset <- Rossi[
  unlist(sapply(is_large, function(x) which(x == Rossi$id))), ]
Rossi_subset <- Rossi_subset[nrow(Rossi_subset):1, ]
Rossi_subset[!duplicated(Rossi_subset$id), ]

## ------------------------------------------------------------------------
tmp <- xtabs(~Rossi$prio[!duplicated(Rossi$id)]) 
plot(as.numeric(names(tmp)), c(tmp), ylab = "frequency", type = "h", 
     xlab = "Number of prior convictions")

## ------------------------------------------------------------------------
tmp <- xtabs(~log(Rossi$prio[!duplicated(Rossi$id)] + 1))
plot(as.numeric(names(tmp)), c(tmp), ylab = "frequency", type = "h", 
     xlab = "Log(Number of prior convictions + 1)")

## ------------------------------------------------------------------------
dd_rossi_trans <- ddhazard(
  Surv(start, stop, event) ~ fin + age + log(prio + 1) + employed.cumsum, 
  data = Rossi, id = Rossi$id, by = 1, max_T = 52, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

plot(dd_rossi_trans)

## ---- echo = FALSE-------------------------------------------------------
hats <- hatvalues(dd_rossi_trans)

hats_stacked <- stack_hats(hats)

## ----rossi_hat_plot, echo = FALSE, fig.path = "figure/replot-"-----------
# We will highlight the individuals with the highest hatvalues
is_large <- 
  names(head(sort(max_hat, decreasing = T), 5))

# Plot hat values
plot(range(hats_stacked$interval_n), c(0, 0.22), type = "n",
     xlab = "Interval number", ylab = "Hat value")

invisible(
  by(hats_stacked, hats_stacked$id, function(rows){
    has_large <- rows$id[1] %in% is_large
    col <- if(has_large) "Red" else "Black"
    lines(rows$interval_n, rows$hat_value, lty = 2, 
          col = col)
    
    if(has_large){
      pch <- as.character(which(rows$id[1] == is_large))
      points(rows$interval_n, rows$hat_value, pch = pch, col = col)
    }
  }))

## ------------------------------------------------------------------------
# Fit the two models
f1 <- ddhazard(
  Surv(start, stop, event) ~ fin + age + prio + employed.cumsum, 
  data = Rossi, id = Rossi$id, by = 1, max_T = 52, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

f2 <- ddhazard(
  Surv(start, stop, event) ~ fin + age + log(prio + 1) + employed.cumsum , 
  data = Rossi, id = Rossi$id, by = 1, max_T = 52, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

# Compute residuals
res1 <- residuals(f1, type = "pearson")
res2 <- residuals(f2, type = "pearson")

# Compute logistic loss 
log_error1 <- unlist(
  lapply(res1$residuals, function(x) 
    ifelse(x[, "Y"] == 1, log(x[, "p_est"]), log(1 - x[, "p_est"]))))
log_error2 <- unlist(
  lapply(res2$residuals, function(x) 
    ifelse(x[, "Y"] == 1, log(x[, "p_est"]), log(1 - x[, "p_est"]))))

# Compare mean
mean(log_error1) # Mean error for first model
mean(log_error2) # Mean error for second model

## ------------------------------------------------------------------------
summary(log_error1)
summary(log_error2)

## ------------------------------------------------------------------------
load("Diagnostics/whas500.RData")

hist(whas500$lenfol[whas500$fstat == 1], breaks = 20, 
     xlab = "Time of death", main = "")

## ------------------------------------------------------------------------
# We only keep some of the columns
whas500 <- whas500[
  , c("id", "lenfol", "fstat", "gender",  "age", "bmi", "hr", "cvd")]

# First rows
head(whas500, 10)

# Summary stats
summary(whas500[, c("age", "bmi", "hr", "gender",  "cvd")])

## ------------------------------------------------------------------------
dd_whas <- ddhazard(
  Surv(lenfol, fstat) ~ gender + age + bmi + hr + cvd,
  data = whas500, by = 100, max_T = 2000, 
  Q_0 = diag(10000, 6), Q = diag(.1, 6))

plot(dd_whas)

## ------------------------------------------------------------------------
dd_whas <- ddhazard(
  Surv(lenfol, fstat) ~ age + bmi + hr + cvd,
  data = whas500, by = 100, max_T = 2000, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

plot(dd_whas)

## ------------------------------------------------------------------------
dd_whas_no_age <- ddhazard(
  Surv(lenfol, fstat) ~ bmi + hr + cvd,   # No age
  data = whas500, by = 100, max_T = 1700, 
  Q_0 = diag(10000, 4), Q = diag(.1, 4))

plot(dd_whas_no_age)

## ------------------------------------------------------------------------
obs_res <- residuals(dd_whas_no_age, type = "pearson")

## ------------------------------------------------------------------------
# We have matrix for each interval
length(obs_res$residuals)    

# Shows the structure of the matrices. We only print take the first 5 matrices
str(obs_res$residuals[1:5])  

# Print the first entries of the first interval
head(obs_res$residuals[[1]]) 

## ----stack_res, echo = FALSE---------------------------------------------
stack_residuals <- function(fit, resids){
  if(!inherits(resids, "fahrmeier_94_res"))
    stop("Residuals must have class 'fahrmeier_94_res'")
  if(!inherits(fit, "fahrmeier_94"))
    stop("fit must have class 'fahrmeier_94'")
  
  # Stack the residuals
  resids_stacked <- 
    data.frame(do.call(rbind, resids[[1]]), row.names = NULL)
  
  # Add the interval number and id
  n_rows <- unlist(lapply(resids$residuals, nrow))
  interval_n <- unlist(sapply(1:length(n_rows), function(i) rep(i, n_rows[i])))
  
  resids_stacked <- cbind(
    resids_stacked, 
    interval_n = interval_n,
    id = fit$id[resids_stacked$row_num])
  
  # Sort by id and interval number
  resids_stacked <- 
    resids_stacked[order(resids_stacked$id, resids_stacked$interval_n), ]
  
  resids_stacked
}

## ------------------------------------------------------------------------
resids_stacked <- stack_residuals(fit = dd_whas_no_age, resids = obs_res)

# print the first entries
head(resids_stacked, 10)

## ------------------------------------------------------------------------
# Add age variable
resids_stacked$age <- 
  whas500$age[resids_stacked$row_num]

# Stratify
age_levels <- quantile(whas500$age, seq(0, 1, by = .2))
age_levels
resids_stacked$age_cut <- cut(resids_stacked$age, age_levels)

# Compute the means 
cut_means <- 
  tapply(resids_stacked$residuals, 
         list(resids_stacked$interval_n, resids_stacked$age_cut), 
         mean)

head(cut_means)

# Plot against time
colfunc <- colorRampPalette(c("Black", "Blue"))
cols <- colfunc(ncol(cut_means))

matplot(dd_whas_no_age$times[-1], apply(cut_means, 2, cumsum), 
        type = "l", col = cols, xlab = "Time", lwd = 2,
        lty = 1, ylab = "Cumulated mean Pearson residuals")
abline(h = 0, lty = 2)

legend("topleft", bty = "n", 
       lty = rep(1, ncol(cut_means)),
       legend = colnames(cut_means), 
       col = cols, lwd = 2,
       cex = par()$cex * .8)

## ------------------------------------------------------------------------
stat_res <- residuals(dd_rossi, type = "std_space_error")

str(stat_res)

## ---- echo=FALSE---------------------------------------------------------
stopifnot(all(abs(stat_res$residuals) < 1.96))

## ------------------------------------------------------------------------
plot(stat_res, mod = dd_rossi, p_cex = .75, ylim = c(-2, 2))

## ------------------------------------------------------------------------
# Get non-standardized residuals
stat_res <- residuals(dd_rossi, type = "space_error")

# Standardize marginally
for(i in 1:nrow(stat_res$residuals))
  stat_res$residuals[i, ] <- stat_res$residuals[i, ] /
    sqrt(diag(stat_res$Covariances[, , i]))

# Plot
par(mfcol = c(2, 3))
for(i in 1:ncol(stat_res$residuals))
  plot(stat_res, mod = dd_rossi, p_cex = .75, cov_index = i, 
       ylab = colnames(stat_res$residuals)[i], 
       ylim = c(-2, 2))

## ------------------------------------------------------------------------
stat_res <- residuals(dd_whas, type = "std_space_error")

plot(stat_res, mod = dd_whas, ylim = c(-4, 4), p_cex = .8)

## ---- echo=FALSE---------------------------------------------------------
if(requireNamespace("dichromat", quietly = T, warn.conflicts = F)){
  cols <- c("#000000", dichromat::colorschemes$Categorical.12[c(6, 8, 10)])
} else
  cols <- c("#000000", "#009E73", "#e79f00", "#9ad0f3")

## ------------------------------------------------------------------------
sim_func <- with(environment(ddhazard), test_sim_func_logit)

## ------------------------------------------------------------------------
# Simulate the coefficients
set.seed(556189)
betas <- matrix(rnorm(21 * 4), ncol = 4)
betas[, 1] <- betas[, 1] * 0.25  # reduce the variance of the intercept
betas <- apply(betas, 2, cumsum) # accumulate the innovations
betas[, 1] <- betas[, 1] - 4     # we reduce the intercept

# Plot the simulated coefficients
matplot(betas, col = cols, lty = 1, type = "l")

## ------------------------------------------------------------------------
# Simulate series
sim_dat <- sim_func(
  n_series = 500,      # number of individuals
  t_max = 20,          # the last stop time
  x_range = 2,         # what is the uniform range to draw from
  x_mean = 0,          # the mean of the uniform covariates
  n_vars = 3,          # 4 - 1 for the intercept
  lambda = 1/10,       # lambda in the exponential distribution for time 
                       # between updates of covariate vectors
  betas = betas)

## ------------------------------------------------------------------------
head(sim_dat$res, 10)

## ------------------------------------------------------------------------
f1 <- ddhazard(
  Surv(tstart, tstop, event) ~ x1 + x2 + x3,
  sim_dat$res, by = 1, max_T = 20, id = sim_dat$res$id, 
  Q_0 = diag(10000, 4), Q = diag(.1, 4))

matplot(betas, col = cols, lty = 1, type = "l")
matplot(f1$state_vecs, col = cols, lty = 2, type = "l", add = T)

## ------------------------------------------------------------------------
stat_res <- residuals(f1, type = "std_space_error")

plot(stat_res, mod = f1, p_cex = .8)

## ------------------------------------------------------------------------
qqnorm(c(stat_res$residuals), pch = 16, cex = .8)
qqline(c(stat_res$residuals))

## ---- echo=FALSE, fig.height=8-------------------------------------------
# CHECK: arguments below match those above
par(mfcol = c(3, 3))

set.seed(189780)
for(i in 1:3){
  # Simulate the coefficients
  betas <- matrix(rnorm(21 * 4), ncol = 4)
  betas[, 1] <- betas[, 1] * 0.25  
  betas <- apply(betas, 2, cumsum) 
  betas[, 1] <- betas[, 1] - 4
  
  # Simulate series
  sim_dat <- sim_func(
    n_series = 500,      
    t_max = 20,          
    x_range = 2,         
    x_mean = 0,          
    n_vars = 3,          
    lambda = 1/10,        
                         
    betas = betas)
  
  # Fit
  f1 <- ddhazard(
    Surv(tstart, tstop, event) ~ x1 + x2 + x3,
    sim_dat$res, by = 1, max_T = 20, id = sim_dat$res$id, 
    Q_0 = diag(10000, 4), Q = diag(.1, 4))
  
  # Plot coefficients
  matplot(betas, col = cols, lty = 1, type = "l",
          ylim = range(betas, f1$state_vecs))
  matplot(f1$state_vecs, col = cols, lty = 2, type = "l", add = T)
  
  # Plot errors
  stat_res <- residuals(f1, type = "std_space_error")
  plot(stat_res, mod = f1, p_cex = .8)
  
  # QQ-plos
  qqnorm(c(stat_res$residuals), pch = 16, cex = .8)
  qqline(c(stat_res$residuals))
}


## ------------------------------------------------------------------------
obs_res <- residuals(dd_rossi, type = "pearson")

## ------------------------------------------------------------------------
# Stack residuals
resids_stacked <- stack_residuals(fit = dd_rossi, resids = obs_res)

# Compute cummulated residuals
resids_stacked$residuals_cum <- unlist(tapply(
  resids_stacked$residuals, resids_stacked$id, cumsum))

# Plot the cumulated residuals for each individual
plot(c(1, 52), range(resids_stacked$residuals_cum), type = "n", 
     xlab = "Interval number", ylab = "Cumulated Person residuals")
invisible(
  tapply(resids_stacked$residuals_cum, resids_stacked$id, lines, 
         col = gray(0, alpha = .2)))

## ------------------------------------------------------------------------
max_cum <- tapply(resids_stacked$residuals_cum, resids_stacked$id, max)

is_max <- names(max_cum)[order(max_cum, decreasing = T)[1:5]]
is_max # Id of the those with the largest values

## ------------------------------------------------------------------------
Rossi_max_subset <- Rossi[Rossi$id %in% is_max, ]
Rossi_max_subset <- Rossi_max_subset[nrow(Rossi_max_subset):1, ]
Rossi_max_subset[!duplicated(Rossi_max_subset$id), ]

## ------------------------------------------------------------------------
age <- Rossi$age[!duplicated(Rossi$id)]
hist(age, breaks = 20)

## ------------------------------------------------------------------------
summary(Rossi$prio[!duplicated(Rossi$id)])

## ------------------------------------------------------------------------
min_cum <- tapply(resids_stacked$residuals_cum, resids_stacked$id, min)

is_min <- names(min_cum)[order(min_cum)[1:5]]
is_min

# We print the last record each of these
Rossi_min_subset <- Rossi[Rossi$id %in% is_min, ]
Rossi_min_subset <- Rossi_min_subset[nrow(Rossi_min_subset):1, ]
Rossi_min_subset[!duplicated(Rossi_min_subset$id), ]

## ------------------------------------------------------------------------
Rossi_sub <- Rossi[!Rossi$id %in% c(is_max, is_min), ]

dd_rossi_tmp <- ddhazard(
  Surv(start, stop, event) ~ fin + age + prio + employed.cumsum, 
  data = Rossi_sub, id = Rossi_sub$id, by = 1, max_T = 52, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

par(mfcol = c(2, 3))
for(i in 1:5){
  plot(dd_rossi, cov_index = i)
  plot(dd_rossi_tmp, cov_index = i, col = "red", add = T)
}

## ------------------------------------------------------------------------
hats <- hatvalues(dd_whas)
hats_stacked <- stack_hats(hats)

## ------------------------------------------------------------------------
# Compute cumulated hat values
hats_stacked$hats_cum <- unlist(tapply(
  hats_stacked$hat_value, hats_stacked$id, cumsum))

# Plot the cumulated residuals for each individual
plot(c(1, 20), range(hats_stacked$hats_cum), type = "n", 
     xlab = "Interval number", ylab = "Cumulated hat values")
invisible(
  tapply(hats_stacked$hats_cum, hats_stacked$id, lines, 
         col = gray(0, alpha = .2)))

## ------------------------------------------------------------------------
max_cum <- tapply(hats_stacked$hats_cum, hats_stacked$id, max)
is_max <- order(max_cum, decreasing = T)[1:3]
is_max

# The records for these
whas500[is_max, ]

## ------------------------------------------------------------------------
# Averages of hat values
hats_stacked$hats_avg <- hats_stacked$hats_cum / hats_stacked$interval_n
  
# Plot against time
plot(c(1, 20), range(hats_stacked$hats_avg), type = "n", 
     xlab = "Interval number", ylab = "Avg. hat values")
invisible(
  tapply(hats_stacked$hats_avg, hats_stacked$id, lines, 
         col = gray(0, alpha = .2)))

## ------------------------------------------------------------------------
max_avg <- tapply(hats_stacked$hats_avg, hats_stacked$id, max)
is_max_avg <- order(max_avg, decreasing = T)[1:5]
is_max_avg

# The records for these
whas500[is_max_avg, ]

## ------------------------------------------------------------------------
# Setup parameters for the plot
cols <- rep("Black", 500)
cols[1:500 %in% is_max_avg] <- "Blue"
cols[1:500 %in% is_max] <- "Red"

cexs <- ifelse(cols != "Black", par()$cex * 1.25, par()$cex * .75)
pchs <- ifelse(whas500$fstat == 1 & whas500$lenfol <= 2000, 16, 1)

plot(whas500[, c("age", "hr", "bmi")], pch = pchs, cex = cexs, col = cols)

## ------------------------------------------------------------------------
plot(whas500$lenfol, whas500$hr, col = cols, pch = pchs, 
     xlab = "Follow-up time", ylab = "Heart rate")
plot(whas500$lenfol, whas500$age, col = cols, pch = pchs, 
     xlab = "Follow-up time", ylab = "Age")

## ------------------------------------------------------------------------
dd_whas <- ddhazard(
  Surv(lenfol, fstat) ~ age + bmi + hr + cvd,
  data = whas500, by = 100, max_T = 1700, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

dd_whas_no_extreme <- ddhazard(
  Surv(lenfol, fstat) ~ age + bmi + hr + cvd,
  
  data = whas500[-c(is_max, is_max_avg), ], # we exclude the "extreme" persons
  
  by = 100, max_T = 1700, 
  Q_0 = diag(10000, 5), Q = diag(.1, 5))

## ------------------------------------------------------------------------
par(mfcol = c(2,3))
for(i in 1:5){
  plot(dd_whas, cov_index = i)
  plot(dd_whas_no_extreme, cov_index = i, add = T, col = "DarkBlue")
}

## ----stack_hats, eval = FALSE--------------------------------------------
#  stack_hats <- function(hats){
#    # Stack
#    resids_hats <- data.frame(do.call(rbind, hats), row.names = NULL)
#  
#    # Add the interval number
#    n_rows <- unlist(lapply(hats, nrow))
#    interval_n <- unlist(sapply(1:length(n_rows), function(i) rep(i, n_rows[i])))
#  
#    resids_hats <- cbind(resids_hats, interval_n = interval_n)
#  
#    # Sort by id and interval number
#    resids_hats <-
#      resids_hats[order(resids_hats$id, resids_hats$interval_n), ]
#  
#    resids_hats
#  }

## ----stack_res, eval = FALSE---------------------------------------------
#  stack_residuals <- function(fit, resids){
#    if(!inherits(resids, "fahrmeier_94_res"))
#      stop("Residuals must have class 'fahrmeier_94_res'")
#    if(!inherits(fit, "fahrmeier_94"))
#      stop("fit must have class 'fahrmeier_94'")
#  
#    # Stack the residuals
#    resids_stacked <-
#      data.frame(do.call(rbind, resids[[1]]), row.names = NULL)
#  
#    # Add the interval number and id
#    n_rows <- unlist(lapply(resids$residuals, nrow))
#    interval_n <- unlist(sapply(1:length(n_rows), function(i) rep(i, n_rows[i])))
#  
#    resids_stacked <- cbind(
#      resids_stacked,
#      interval_n = interval_n,
#      id = fit$id[resids_stacked$row_num])
#  
#    # Sort by id and interval number
#    resids_stacked <-
#      resids_stacked[order(resids_stacked$id, resids_stacked$interval_n), ]
#  
#    resids_stacked
#  }

