context("BIQQ clue analysis")

es_fit <-
  init_biqq(model_code = stan_es,
            data = data_es_init)

sims <- 1 # Number of simulations
N <- nrow(data_es) # Number of cases

# Generate hyperparameters for Stan
q0_alpha_beta <- mapply(m = c(q.a0 = 0.1, # Assumptions on mean of q0's
                              q.b0 = 0.1,
                              q.c0 = 0.05,
                              q.d0 = 0.3),
                        sd = rep(.01, times = 4),
                        FUN = beta_prior)

q1_alpha_beta <- mapply(m = c(q.a1 = 0.95, # Assumptions on mean of q1's
                              q.b1 = 0.9,
                              q.c1 = 0.475,
                              q.d1 = 0.5),
                        sd = rep(.01, times = 4),
                        FUN = beta_prior)
# Rstan setup options
rstan_options(auto_write = TRUE)
options(mc.cores = 1)

# Run the analysis
betas <-
  parallel::mclapply(X          = rep(0:N, each = sims),
                     FUN        = biqq::biqq_es_k,
                     fit        = es_fit,
                     clue_data  = data_es,
                     q0_alpha   = q0_alpha_beta["alpha",],
                     q0_beta    = q0_alpha_beta["beta",],
                     q1_alpha   = q1_alpha_beta["alpha",],
                     q1_beta    = q1_alpha_beta["beta",],
                     chains     = 1,
                     cores      = 1,
                     extract_pars = "abcd",
                     out_fun = function(x) { mean(x$abcd[,2] - x$abcd[,1]) })

betas <- matrix(unlist(betas), ncol = sims, byrow = TRUE)

test_that("Clue analysis for ES example works", {
  expect_equal(dim(betas), c(21,1))
})
