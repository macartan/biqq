context("BIQQ plot posterior")

nr_fit <- init_biqq(model_code = stan_nr, data = data_nr_init)
nr_example <- biqq_nr(fit = nr_fit, chains = 4, cores = 1)

es_fit <- init_biqq(model_code = stan_es, data = data_es_init)
es_example <- biqq_es(fit = es_fit, chains = 4, cores = 1)

test_that("plot of posterior works",{
  # Natural resources example
  expect_equal(class(plot_biqq(fit = nr_example, title = "Posteriors | Both clues weak")),
               c("gg", "ggplot"))
  # Electoral systems example
  expect_equal(class(plot_biqq(fit = es_example, title = "")),
               c("gg", "ggplot"))
})

