context("Testing test_utils")

set.seed(101010)

test_that("Testing util functions to sim for test", {
  test_pairs <- cbind("new funcs" = c(get_exp_draw(), get_unif_draw(), get_norm_draw()),
                     "funcs" = c(rexp, runif, rnorm))

  for(i in seq_len(nrow(test_pairs))){
    new_func <- test_pairs[i, 1][[1]]
    func <- test_pairs[i, 2][[1]]

    expect_equal(length(new_func(1)), 1)
    expect_equal(length(new_func(10)), 10)
    expect_equal(length(new_func(100)), 100)

    seed <- 101
    set.seed(seed)
    new_func(re_draw = T)
    res1 <- new_func(10)
    set.seed(seed)
    res2 <- func(10)
    expect_equal(res1, res2)
  }
})

test_that("Testing util functions to sim series for tests", {
  skip_on_cran()

  set.seed(4321)
  n_series <- 1e4
  t_max <- 10

  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    tmp = func(n_series, t_max = t_max)$res

    expect_equal(unique(tmp[, "id"]), seq_len(n_series))

    expect_true(all(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts > 10)) < 2))
    expect_true(any(tapply(tmp[, "tstop"], tmp[, "id"], function(ts) sum(ts >= 10)) > 0))

    expect_true(all(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) < 2))
    expect_true(any(tapply(tmp[, "event"], tmp[, "id"], function(ev) sum(ev)) > 0))

    expect_true(all(tmp[, "tstart"] < tmp[, "tstop"]))

    expect_true(all(tmp[-1, "id"] != tmp[-nrow(tmp), "id"] |
                      tmp[-1, "tstart"] == tmp[-nrow(tmp), "tstop"]))

    expect_true(all(tmp[-1, "id"] == tmp[-nrow(tmp), "id"] |
                      tmp[-nrow(tmp), "tstop"] > 4 |
                      tmp[-nrow(tmp), "event"]))

    expect_true(all(!(tmp$res$event == 1 &
                        tmp$res$tstop > t_max)))
  }
})


test_that("Whether covs are fixed or not depends on is fixed argument", {
  set.seed(1999293)
  inter_start <- rnorm(1)
  beta_start <- rnorm(10)
  for(func in c(test_sim_func_logit, test_sim_func_exp)){
    sims <- func(n_series = 1e1, n_vars = 10, beta_start = beta_start,
                 intercept_start = inter_start, t_max = 10)

    expect_true(all(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0)))

    sims <- test_sim_func_exp(n_series = 1e2, n_vars = 10, beta_start = beta_start,
                              intercept_start = inter_start, t_max = 10, is_fixed = c(1, 4))
    expect_equal(apply(sims$betas, 2, function(x) sum(duplicated(x)) == 0),
                 !seq_len(11) %in% c(1, 4))
  }
})

test_that("Sim functions gives previous results", {
  set.seed(897730)
  sim_logit <- test_sim_func_logit(1e2, n_vars = 3)
  sim_exp <- test_sim_func_exp(1e2, n_vars = 3)

  # save_to_test(sim_logit, "util_sim_logit")
  # save_to_test(sim_exp, "util_sim_exp")

  expect_equal(sim_logit, read_to_test("util_sim_logit"), tolerance = 1.490116e-08)
  expect_equal(sim_exp, read_to_test("util_sim_exp"), tolerance = 1.490116e-08)
})