## ----setup, include=FALSE------------------------------------------------
knitr::knit_hooks$set(
  mySettings  = function(before, options, envir){
    if (before && options$mySettings){ 
      par(
        mar = c(10, 10, 4, 4),
        bty = "n",
        xaxs = "i",
        pch=16,
        cex= (cex <- .4),
        cex.axis = .8 / cex,
        cex.lab = .8 / cex,
        lwd= 1)
      options(digits = 3, width = 80, warn = -1)
    }})

knitr::opts_chunk$set(echo = TRUE, mySettings=TRUE, fig.height=3.5, fig.width = 6,
                      warning = F, message = F)
knitr::opts_knit$set(warning = F, message = F)

## ----echo=FALSE---------------------------------------------------------------
tryCatch({
  current_sha <- paste0("@", httr::content(
    httr::GET("https://api.github.com/repos/boennecd/dynamichazard/git/refs/heads/master")
    )$object$sha)
}, error = function(...){ current_sha <<- "" })

stopifnot(length(current_sha) > 0 && class(current_sha) == "character")

current_version <- paste0("boennecd/dynamichazard@", current_sha)

## -----------------------------------------------------------------------------
current_version # The string you need to pass devtools::install_github

## ----eval=FALSE---------------------------------------------------------------
#  devtools::install_github(current_version)

## -----------------------------------------------------------------------------
# PBC data set from survival with time variying covariates
# Details of tmerge are not important in this scope. The code is included
# to make you able to reproduce the results
# See: https://cran.r-project.org/web/packages/survival/vignettes/timedep.pdf
library(survival)
temp <- subset(pbc, id <= 312, select=c(id, sex, time, status, edema, age))
pbc2 <- tmerge(temp, temp, id=id, death = event(time, status))
pbc2 <- tmerge(pbc2, pbcseq, id=id, albumin = tdc(day, albumin),
               protime = tdc(day, protime), bili = tdc(day, bili))
pbc2 <- pbc2[, c("id", "tstart", "tstop", "death", "sex", "edema", 
                 "age", "albumin", "protime", "bili")]

## -----------------------------------------------------------------------------
head(pbc2)

## -----------------------------------------------------------------------------
(ex <- pbc2[pbc2$id == 282, ])

## -----------------------------------------------------------------------------
glm_simple <- glm(death == 2 ~ age + edema + log(albumin) + log(protime) + 
                    log(bili), binomial, pbc2)

glm_simple$coefficients

## -----------------------------------------------------------------------------
library(dynamichazard)
pbc2_big_frame <- get_survival_case_weigths_and_data(
  Surv(tstart, tstop, death == 2) ~ age + edema + log(albumin) + log(protime) +
    log(bili), data = pbc2, id = pbc2$id, by = 100, max_T = 3600, 
  use_weights = F)
pbc2_big_frame <- pbc2_big_frame$X

## -----------------------------------------------------------------------------
pbc2_big_frame[pbc2_big_frame$id == 282, ]

## -----------------------------------------------------------------------------
pbc2_small_frame <- get_survival_case_weigths_and_data(
  Surv(tstart, tstop, death == 2) ~ age + edema + log(albumin) + log(protime) +
    log(bili), data = pbc2, id = pbc2$id, by = 100, max_T = 3600, 
  use_weights = T)
pbc2_small_frame <- pbc2_small_frame$X

## -----------------------------------------------------------------------------
pbc2_small_frame[pbc2_small_frame$id == 282, ]

## -----------------------------------------------------------------------------
pbc2[pbc2$id == 268, ] # the orginal data
pbc2_small_frame[pbc2_small_frame$id == 268, ] # new data

## -----------------------------------------------------------------------------
glm_fit_big <- glm(Y ~ age + edema + log(albumin) + log(protime) + 
                    log(bili), binomial, pbc2_big_frame)
glm_fit_small <- glm(Y ~ age + edema + log(albumin) + log(protime) + 
                      log(bili), binomial, pbc2_small_frame, 
                     weights = pbc2_small_frame$weights) # <-- we use weights

## -----------------------------------------------------------------------------
all.equal(glm_fit_big$coefficients, glm_fit_small$coefficients)

## -----------------------------------------------------------------------------
knitr::kable(
  rbind("glm with bins" = glm_fit_big$coefficients, 
        "glm without bins" = glm_simple$coefficients, 
        "Sds with bins" = 
          summary(glm_fit_big)[["coefficients"]][, "Std. Error"],
        "Sds without bins" = 
          summary(glm_simple)[["coefficients"]][, "Std. Error"]),
  row.names = T)

## -----------------------------------------------------------------------------
static_glm_fit <- static_glm(
  Surv(tstart, tstop, death == 2) ~ age + edema + log(albumin) + log(protime) +
    log(bili), data = pbc2, id = pbc2$id, by = 100, max_T = 3600)

all.equal(static_glm_fit$coefficients, glm_fit_big$coefficients)

## ----message=FALSE------------------------------------------------------------
library(mgcv, quietly = T)
spline_fit <- gam(
  formula = Y ~ 
    # cr is cubic spline with k knots
    s(t, bs = "cr", k = 3, by = age) + 
    s(t, bs = "cr", k = 3, by = edema) + 
    s(t, bs = "cr", k = 5, by = log(albumin)) + 
    s(t, bs = "cr", k = 5, by = log(protime)) + 
    s(t, bs = "cr", k = 5, by = log(bili)),
  family = binomial, data = pbc2_big_frame,
  method = "GCV.Cp") # estimate smoothing parameters with generalized cross 
                     # validation  
                  

## ---- fig.height = 4.5, mySettings=FALSE, fig.cap="Plots of estimated effects in the GAM model. The effective degree of freedom is noted in the parentheses on the y axis and is computed given the number knots and final tuning parameter for spline function. For instance, `s(t,2.43):age` means that the effective degrees of freedom for `age` is `2.43`"----
plot(spline_fit, scale = 0, page = 1, rug = F)

## -----------------------------------------------------------------------------
glm_fit_big$coefficients

## -----------------------------------------------------------------------------
spline_fit$coefficients["(Intercept)"]

## ---- message=FALSE, fig.height = 7, mySettings=FALSE-------------------------
library(timereg)
cox_fit<- timecox(Surv(tstart / 365, tstop / 365, death == 2) ~ age + edema +
                        log(albumin) + const(log(protime)) + log(bili), pbc2,
                  start.time=0, 
                  max.time = 3000 / 365, # <-- decreased
                  id = pbc2$id, bandwidth = 0.35)

par(mfcol = c(3, 2))
plot(cox_fit)

## ---- message=FALSE, fig.height = 6, mySettings=FALSE-------------------------
summary(cox_fit)

## ---- fig.height = 4.5, mySettings=FALSE--------------------------------------
library(dynamichazard)
dd_fit <- ddhazard(Surv(tstart, tstop, death == 2) ~ age + edema +
                        log(albumin) + log(protime) + log(bili), pbc2,
                   id = pbc2$id, by = 100, max_T = 3600, 
                   Q_0 = diag(rep(10, 6)), Q = diag(rep(0.001, 6)))

plot(dd_fit)

## ---- fig.height = 4.5, mySettings=FALSE--------------------------------------
dd_fit <- ddhazard(Surv(tstart, tstop, death == 2) ~ age + edema +
                        log(albumin) + log(protime) + log(bili), pbc2,
                   id = pbc2$id, by = 100, max_T = 3600,
                   a_0 = glm_fit_big$coefficients,
                   control = list(NR_eps = 0.1), # <-- tolerance in correction step
                   Q_0 = diag(rep(1, 6)),        # <-- decreased Q_0
                   Q = diag(rep(0.001, 6)))

# Plot result
plot(dd_fit)

## ---- fig.height = 4.5, mySettings=FALSE--------------------------------------
dd_fit_UKF <- ddhazard(Surv(tstart, tstop, death == 2) ~ age +
                         edema + log(albumin) + log(protime) + log(bili), pbc2,
                   id = pbc2$id, by = 100, max_T = 3600, 
                   Q_0 = diag(rep(1, 6)), Q = diag(rep(0.01, 6)),
                   control = list(method = "UKF", beta = 0, alpha = 1,
                                  eps = 0.1, n_max = 1e4))

plot(dd_fit_UKF)

## ---- eval=FALSE--------------------------------------------------------------
#  # Not run
#  tmp_file <- file("pick_some_file_name.txt")
#  sink(tmp_file)
#  dd_fit_UKF <- ddhazard(Surv(tstart, tstop, death == 2) ~ age +
#                           edema + log(albumin) + log(protime) + log(bili), pbc2,
#                     id = pbc2$id, by = 100, max_T = 3600,
#                     Q_0 = diag(rep(1, 6)), Q = diag(rep(0.01, 6)),
#                     control =
#                       list(method = "UKF", beta = 0, alpha = 1,
#                            debug = T)) # <-- prints information in each iteration
#  sink()
#  close(tmp_file)

## -----------------------------------------------------------------------------
diag(dd_fit$Q)

## ---- fig.height = 4.5, mySettings=FALSE--------------------------------------
dd_fit_UKF <- ddhazard(
  Surv(tstart, tstop, death == 2) ~ age +
    edema + log(albumin) + log(protime) + log(bili), pbc2,
  id = pbc2$id, by = 100, max_T = 3600, 
  Q_0 = diag(c(0.001, 0.00001, rep(0.001, 4))) * 100, # <-- decreased
  Q = diag(rep(0.0001, 6)),                           # <-- decreased
  control = 
    list(method = "UKF", beta = 0, alpha = 1, eps = 0.01))

plot(dd_fit_UKF)

## ---- warning=FALSE-----------------------------------------------------------
dd_fit <- ddhazard(
  Surv(tstart, tstop, death == 2) ~ ddFixed(1) + 
    ddFixed(age) + ddFixed(log(albumin)) + edema +ddFixed(log(protime)) + log(bili), 
  pbc2, id = pbc2$id, by = 100, max_T = 3600, 
  Q_0 = diag(rep(1, 2)), Q = diag(rep(0.0001, 2)), control = list(eps = 0.02))

## ---- fig.height = 3.5, mySettings=FALSE--------------------------------------
plot(dd_fit)

## ---- fig.height = 6----------------------------------------------------------
dd_fit$fixed_effects

## ---- fig.height = 3.5, mySettings=FALSE--------------------------------------
dd_fit <- ddhazard(
  Surv(tstart, tstop, death == 2) ~ ddFixed(1) + 
    ddFixed(age) + ddFixed(log(albumin)) + edema +ddFixed(log(protime)) + log(bili), 
  pbc2, id = pbc2$id, by = 100, max_T = 3600, 
  Q_0 = diag(rep(1, 2)), Q = diag(rep(0.0001, 2)), 
  control = list(eps = 0.02, NR_eps = 0.1))

plot(dd_fit)

## -----------------------------------------------------------------------------
dd_fit <- ddhazard(Surv(tstart, tstop, death == 2) ~ 
                     ddFixed(1) + ddFixed(age) + 
                     ddFixed(edema) + ddFixed(log(albumin)) + 
                     ddFixed(log(protime)) + log(bili), pbc2,
                   id = pbc2$id, by = 100, max_T = 3600,
                   order = 2,             # <-- second order
                   Q_0 = diag(c(10, 10)), # <-- needs more elements
                   Q = matrix(0.001))

plot(dd_fit)

dd_fit$fixed_effects

