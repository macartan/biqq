context("BIQQ clue analysis")

fit_biqq <-
  init_biqq(stanmodel = stan_biqq,
            data = data_biqq_init)


test_fit <- biqq(fit_biqq, alpha_prior = c(5,1,1,1), iter = 5000)

test_that("Basic BIQQ runs", {
  expect_equal(mean(summary(test_fit)$summary[,10]), 1, tolerance = .0001, scale = 1)
})

