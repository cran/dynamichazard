context("Testing ddhazard w/ generic things and w/ the the EKF ")

# Test on data set that is one of Farhmiers papers
result <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1,
  control = list(est_Q_0 = F,
                 save_data = F, save_risk_set = F),
  a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
  max_T = 45,
  id = head_neck_cancer$id, order = 1)

result <- result[c("state_vecs", "state_vars","lag_one_cov", "model")]

test_that("get previous results with head_neck", {
  expect_equal(result, read_to_test("ddhazard_head_neck"))
})

test_that("Invalid penalty terms throw error", {
  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = list(denom_term = 0)),
    regexp = "Method not implemented with penalty term \\(control\\$denom_term\\) equal to 0")

  expect_error(
    ddhazard(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, # Use by month intervals
      control = list(denom_term = -1)),
    regexp = "Method not implemented with penalty term \\(control\\$denom_term\\) equal to -1")
})

test_that("Get expected warning when no Q or Q_0 is passed", {
  args <- list(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  tmp <- args
  tmp$Q_0 <- NULL
  expect_warning(
    do.call(ddhazard, tmp),
    paste(sQuote("Q_0"), "not supplied. It has been set to a diagonal matrix with diagonal entries equal to"))

  tmp <- args
  tmp$Q <- NULL
  expect_warning(
    do.call(ddhazard, tmp),
    paste(sQuote("Q"), "not supplied. It has been set to a diagonal matrix with diagonal entries equal to"))
})



test_that("Changing convergence criteria change output",{
  arg_list <- list(
    formula = survival::Surv(stop, event) ~ group,
    data = head_neck_cancer,
    by = 1, # Use by month intervals
    id = head_neck_cancer$id,
    Q_0 = diag(1e5, 2), Q = diag(.1, 2),
    control = list(criteria = "delta_coef", eps = .002))

  suppressMessages(res1 <- do.call(ddhazard, arg_list))
  arg_list$control$criteria <- "delta_likeli"
  suppressMessages(res2 <- do.call(ddhazard, arg_list))

  expect_true(res1$n_iter != res2$n_iter)
})

result_exp <- ddhazard(
  formula = survival::Surv(start, stop, event) ~ group,
  data = head_neck_cancer,
  by = 1, Q_0 = diag(10000, 2),
  Q = diag(1e-3, 2), a_0 = c(0, 0),
  max_T = 30,
  id = head_neck_cancer$id, order = 1,
  model = "exponential")

test_that("exponential model and logit moels hazzard functions differs", {
  result = ddhazard(
    formula = survival::Surv(start, stop, event) ~ group,
    data = head_neck_cancer,
    by = 1,
    control = list(est_Q_0 = F,
                   save_data = F, save_risk_set = F),
    a_0 = rep(0, 2), Q_0 = diag(100000, 2), Q = diag(0.01, 2),
    max_T = 45,
    id = head_neck_cancer$id, order = 1)

  expect_true(result_exp$model != result$model)
  expect_true(result_exp$family$name() != result$family$name())
})

test_that("Testing names of output from ddhazard on head and neck cancer dataset", {
  expect_equal(colnames(result_exp$state_vecs), c("(Intercept)", "group1"))
  expect_equal(unlist(dimnames(result_exp$state_vars)), unlist(list(c("(Intercept)", "group1"), c("(Intercept)", "group1"), NULL)))
  expect_equal(unlist(dimnames(result_exp$Q)), rep(c("(Intercept)", "group1"), 2))
  expect_equal(unlist(dimnames(result_exp$Q_0)), rep(c("(Intercept)", "group1"), 2))
})

# plot(result_exp)
result_exp <- result_exp[
  names(result_exp) %in%
    c("state_vecs", "state_vars", "fixed_effects", "Q")]
test_that("Result of exponential model on head_neck_data match previous results", {
  expect_known_value(
    result_exp, file = "head_neck_exp.RDS", update = FALSE,
    check.attributes = FALSE)
})
rm(result_exp)


test_that("Unmacthed control variable throw error",
          expect_error({
            result = ddhazard(
              formula = survival::Surv(start, stop, event) ~ group,
              data = head_neck_cancer,
              by = 1, # Use by month intervals
              a_0 = rep(0, 2), Q_0 = diag(1, 2), # Initial value
              max_T = 45,
              id = head_neck_cancer$id, order = 1,
              control = list(None_existing_parem = 1)
            )}, regexp = "These control parameters are not recognized"))

test_that("Different non-integer time_scales gives the same result with ddhazard", {
  skip_on_cran()
  scales <- exp(seq(-2, 2, .05))

  fit_exp <- expression(
    ddhazard(
      formula = survival::Surv(start * .by, stop * .by, event) ~ group,
      data = head_neck_cancer,
      by = .by,
      control = list(est_Q_0 = F, save_data = F),
      a_0 = rep(0, 2), Q_0 = diag(1e2, 2), Q = diag(1e-2 / .by, 2),
      id = head_neck_cancer$id, order = 1))

  .by <- scales[1]
  f1 <- eval(fit_exp)
  for(.by in scales[-1]){
    f2 <- eval(fit_exp)
    info <- paste("by =", .by)
    expect_equal(f1$risk_set$risk_sets, f2$risk_set$risk_sets, info = info)
    expect_equal(f1$risk_set$is_event_in, f2$risk_set$is_event_in, info = info)
    expect_equal(f1$state_vecs, f2$state_vecs, info = info)
  }
})

test_that(
  "Old expoential models gives the same results and yields expected message", {
    args <- list(
      formula = survival::Surv(start, stop, event) ~ group,
      data = head_neck_cancer,
      by = 1, Q_0 = diag(10000, 2),
      Q = diag(1e-3, 2), a_0 = c(0, 0),
      max_T = 30,
      id = head_neck_cancer$id, order = 1, control = list(eps = .1))

    args$model <- "exponential"
    expect_silent(f1 <- do.call(ddhazard, args))

    for(m in c("exp_bin", "exp_clip_time", "exp_clip_time_w_jump")){
      eval(bquote({
        args$model <- .(m)
        expect_message(
          f2 <- do.call(ddhazard, args),
          .(paste0(".", m, ". is not used after version 0\\.5\\.0\\.")))
        expect_equal(f1[c("state_vars", "state_vecs")],
                     f2[c("state_vars", "state_vecs")],
                     info = .(m))
      }))
    }

  })

########
# Test on simulated data

test_that("Result of exponential model gives previous results w/ simulated data", {
  args <- list(
    formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
    data = exp_sim_500$res,
    by = 1,
    Q_0 = diag(1e5, 11),
    Q = diag(1e-3, 11),
    control = list(
      save_data = F, save_risk_set = F,
      method = "EKF"),
    max_T = 10,
    id = exp_sim_500$res$id, order = 1,
    model = "exponential")

  result_exp <- do.call(ddhazard, args)

  # matplot(exp_sim_500$betas, type = "l", lty = 1)
  # matplot(result_exp$state_vecs, lty = 2, type = "l", add = T)
  result_exp <- result_exp[c("state_vars", "state_vecs", "Q")]
  expect_known_value(
    result_exp, "sim_exp.RDS", tolerance = 1e-5, update = FALSE)
})


test_that("Permutating data does not change the results", {
  args <- list(
    Surv(stop, event) ~ group, head_neck_cancer,
    by = 1, max_T = 40,
    Q_0 = diag(rep(10000, 2)), Q = diag(rep(0.1, 2)))

  r1 <- do.call(ddhazard, args)

  args$control <- c(args$control, list(permu = T))

  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

  #####
  # With fixed effects
  args <- list(
    Surv(tstart, tstop, death == 2) ~ age + ddFixed(edema) +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3600,
    Q_0 = diag(rep(10000, 5)), Q = diag(rep(0.001, 5)),
    control = list(n_threads = 1))

  r1 <- do.call(ddhazard, args)
  args$control <- c(args$control, list(permu = T))
  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))

  #####
  # With weigths
  set.seed(94884214)
  w <- sample(1:3, nrow(pbc2), replace = T)

  args <- list(
    Surv(tstart, tstop, death == 2) ~ age + edema +
      log(albumin) + log(protime) + log(bili), pbc2,
    id = pbc2$id, by = 100, max_T = 3000,
    weights = w,
    Q_0 = diag(rep(10000, 6)), Q = diag(rep(0.001, 6)),
    control = list(n_threads = 1, LR = .6))

  r1 <- do.call(ddhazard, args)
  args$control <- c(args$control, list(permu = T))
  r2 <- do.call(ddhazard, args)

  # plot(r1)
  # plot(r2)

  expect_equal(r1$state_vecs, r2$state_vecs)
  expect_true(is.character(
    all.equal(r1$state_vecs, r2$state_vecs,tolerance = 1e-16)))
})
