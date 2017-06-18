## ----setup, include=FALSE------------------------------------------------
knitr::knit_hooks$set(default_opts = function(before, options, envir) {
    if (before){
      options(digist = 4)
      par(
        mar = c(5, 4, 1, 1),
        bty = "n",
        xaxs = "i",
        pch=16,
        cex= (cex <- .4),
        cex.axis = .8 / cex,
        cex.lab = .8 / cex,
        lwd= 1)
    }
})

options(digist = 4)
knitr::opts_chunk$set(
  echo = TRUE, warning = F, message = F, dpi = 126, fig.height=3.5, fig.width = 6)
knitr::opts_knit$set(warning = F, message = F,  default_opts = T)

## ----echo=FALSE----------------------------------------------------------
tryCatch({
  current_sha <- paste0("@", httr::content(
    httr::GET("https://api.github.com/repos/boennecd/dynamichazard/git/refs/heads/master")
    )$object$sha)
}, error = function(...){ current_sha <<- "" })

stopifnot(length(current_sha) > 0 && class(current_sha) == "character")

current_version <- paste0("boennecd/dynamichazard", current_sha)

## ------------------------------------------------------------------------
current_version # the string you need to pass to devtools::install_github

## ----eval=FALSE----------------------------------------------------------
#  devtools::install_github(current_version)

## ----eval=FALSE----------------------------------------------------------
#  install.packages("dynamichazard")

## ---- echo=FALSE---------------------------------------------------------
library(survival)
source("../R/test_utils.R")

start_fun <- function(t_0 = t_0, t_max = t_max) max(0, runif(1, t_0 - t_max, t_max - 1 - 1e-8))

# set.seed(print(round(runif(1, max = 1e6))))
set.seed(126265)
simple_ex <- test_sim_func_logit(
    n_series = 2e3, 
    beta_start = c(1.5, -1),
    intercept_start = - 3, 
    sds = rep(.5, 3),
    t_max = 28,
    n_vars = 2,
    lambda = 1/5,
    x_range = 1,
    x_mean = 0,
    tstart_sampl_func = start_fun)$res
# sum(simple_ex$event)

stopifnot(any(simple_ex$event[simple_ex$id == 1] == 1) && 
              all(simple_ex$event[simple_ex$id == 2] == 0))

## ---- echo = FALSE-------------------------------------------------------
knitr::kable(head(simple_ex, 10), digits = 4)

## ------------------------------------------------------------------------
library(dynamichazard)
library(survival)
dd_fit_short <- ddhazard(
  Surv(tstart, tstop, event) ~ x1 + x2, # Formula like for coxph from survival
       data = simple_ex,
       by = 1,                          # Length of time intervals
       Q = diag(0.1, 3),                # Covariance matrix in state eqn 
       Q_0 = diag(10000, 3),            # Covariance matrix for initial state
                                        # vector
       max_T = 28,                      # Last time we observe
       id = simple_ex$id                # id of individuals
  )                      

# Print diagonal of covariance matrix
diag(dd_fit_short$Q) 

## ------------------------------------------------------------------------
library(dynamichazard)
library(survival)
dd_fit_wide <- ddhazard(
  Surv(tstart, tstop, event) ~ x1 + x2, 
       data = simple_ex,
       by = 2,                          # increased
       Q = diag(0.1, 3),               
       Q_0 = diag(10000, 3), 
       max_T = 28,                      
       id = simple_ex$id)            

# Print relative differences between diagonal of covariance matrices
Q_short <- dd_fit_short$Q
Q_wide <- dd_fit_wide$Q
diag((Q_wide - Q_short) / Q_short)   

## ---- fig.height= 5------------------------------------------------------
par(mfcol = c(2, 2), mar = c(5, 4, 1, 1))

for(i in 1:3){
  plot(dd_fit_short, cov_index = i, col = "Black")
  plot(dd_fit_wide, cov_index = i, col = "Red", add = T)
}

## ---- eval=FALSE---------------------------------------------------------
#  dynamichazard::ddhazard_app()

## ----binning_fig, echo=FALSE, results="hide", fig.cap = "Illustration of going from event times to binary variables. Each horizontal line represents an individual. A cross indicates that new covariates are observed while a filled circle indicates that the individual have died. Open circles indicates that the individual is right censored. Vertical dashed lines are time interval borders", fig.height=3----
par(cex = .8, mar = c(1, 4, 1, 2))
plot(c(0, 4), c(0, 1), type="n", xlab="", ylab="", axes = F)

abline(v = c(0.5, 1.5, 2.5), lty = 2)

text(1, 0.01, "1st interval", adj = .5)
text(2, 0.01, "2nd interval", adj = .5)

n_series = 7
y_pos = seq(0, 1, length.out = n_series + 2)[-c(1, n_series +2)]

captioner::captioner()("binning_fig")

set.seed(1992)
x_vals_and_point_codes = list(
  cbind(c(0, .8, 2.2, 3, 3.7) + c(.1, rep(0, 4)),
        c(rep(4, 4), 1)),
  cbind(c(0.1, 1, 1.5 + runif(1)),
        c(4, 4, 16)),
  cbind(c(0.1, .8, 1.9, 2 + runif(1, min = 0, max = .25)),
        c(4, 4, 4, 16)),
  cbind(c(0.1, runif(1) + .33), c(4, 16)),
  cbind(c(0.1, .6, 2.1, 3.1 + runif(1)),
        c(4, 4, 4, 16)),
  cbind(c(2, 2.33),
        c(4, 16)), 
  cbind(c(0.1, 1.3),
        c(4, 1)))

x_vals_and_point_codes = x_vals_and_point_codes[
  c(n_series, sample(n_series - 1, n_series - 1, replace = F))]

for(i in seq_along(x_vals_and_point_codes)){
  vals = x_vals_and_point_codes[[i]]
  y = y_pos[i]
  
  xs = vals[, 1]
  n_xs = length(xs)
  
  # add lines
  segments(xs[-n_xs], rep(y, n_xs - 1),
           xs[-1], rep(y, n_xs - 1))
  
  # Add point
  points(xs, rep(y, n_xs), pch = vals[, 2], 
         cex = ifelse(vals[, 2] == 1, par()$cex * 2.5, par()$cex))
  
  # Add numbers
  text(xs, rep(y + .05, n_xs), as.character(0:(n_xs -1)))
}

# add letters
text(rep(0, n_series), rev(y_pos), letters[1:n_series], cex = par()$cex * 1.5)

## ---- echo=FALSE, fig.height=3-------------------------------------------
f <- function(l) exp(-l) * (exp(l) - 1) * (1 - l) / l

par(mfcol = c(1,2), 
    cex = .67, cex.axis = 1, cex.lab = 1,
    lwd = .67)

tmp_cex <- .33
plot(f, xlim = c(1e-2, 1e1), frame = FALSE, yaxt='n',
     xlab = expression(lambda), ylab = "Mean", cex = tmp_cex)
axis(2, at = c(-1, 0, 1))
abline(h = 0, lty = 2)

plot(f, xlim = c(1e-3, 1e3), log = "x", frame = FALSE,
     yaxt='n',
     xlab = expression(lambda), ylab = "Mean", cex = tmp_cex)
axis(2, at = c(-1, 0, 1))
abline(h = 0, lty = 2)

## ---- echo=FALSE---------------------------------------------------------
set.seed(1010101012)
data_frame <- 
  data.frame(id = c(1, 1, 1, 2, 2),
             tstop = c(0, 2, 3, 0, 2), tstart = c(2, 3, 4, 2, 4), y = c(0, 0, 1, 0, 0), 
             x1 = round(rnorm(5), 2), x2 = round(rnorm(5), 2))

## ------------------------------------------------------------------------
# The data we use
head(data_frame)

## ---- eval=FALSE---------------------------------------------------------
#  # The fit
#  fit <- ddhazard(Surv(tstart, tstop, y) ~ x1 + ddFixed(x2), data_frame,
#                  by = 1, # time interval lengths are 1
#                  id = data_frame$id, model = "exponential")

## ---- echo=FALSE---------------------------------------------------------
data_frame_new <- data_frame
data_frame_new[1, 2:3] <- c(.5, 2)
data_frame_new <- rbind(c(1, 0, .5, 0, .43, .33),
                        data_frame_new)

## ------------------------------------------------------------------------
head(data_frame_new)

