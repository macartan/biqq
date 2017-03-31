context("Initiate BIQQ")

# Create model fit for Origins of Electoral Systems example using init_biqq function
es_fit <- init_biqq(model_code = stan_es, data = data_es_init)

# Create model fit for Natural Resource Curse example using init_biqq function
nr_fit <- init_biqq(model_code = stan_nr, data = data_nr_init)

test_that("initiation of BIQQ paper examples work", {
  # Electoral system example
  expect_equivalent(class(es_fit), "stanfit")
  # Natural resources example
  expect_equivalent(class(nr_fit), "stanfit")
})

test_that("initiation of custom BIQQ program works", {

})
