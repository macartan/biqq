#' BIQQ Model for Natural Resource Curse Example
#'
#' Fit BIQQ Model for Natural Resource Curse Example
#'
#' @param fit Stan model fit object
#' @param XYKK XYKK data for model fitting
#' @param alpha_prior Numeric vector of length 4. Dirichlet distribution parameters for proportions of 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param q1bd1_alpha Numeric vector of length 2. Alpha shape parameters for Beta distribution of probabilities of observing first clue given that type is (b)eneficial or (d)estinied. Defaults to \code{c(1,1)}
#' @param q1bd1_beta Numeric vector of length 2. Beta shape parameters for Beta distribution of probabilities of observing first clue given that type is (b)eneficial or (d)estinied. Defaults to \code{c(1,1)}
#' @param q2bd1_alpha Numeric vector of length 2. Alpha shape parameters for Beta distribution of probabilities of observing second clue given that type is (b)eneficial or (d)estinied. Defaults to \code{c(1,1)}
#' @param q2bd1_beta Numeric vector of length 2. Beta shape parameters for Beta distribution of probabilities of observing second clue given that type is (b)eneficial or (d)estinied. Defaults to \code{c(1,1)}
#' @param chains Integer. Number of MC chains. Defaults to 4
#' @param iter Integer. Total number of iterations in each chain. Defaults to 20000,
#' @param warmup Integer. Number of warm-up iterations in each chain. Defaults to 5000
#' @param seed Integer. Number for random number generation. Defaults to 1000
#' @param cores Integer. Number of cores to use for parallel computation. Defaults to 1
#'
#' @keywords stan biqq example natural resources
#'
#' @examples
#'
#' \dontrun{
#'
#' library(biqq)
#'
#' # Create model fit for Natural Resource Curse example using init_biqq function
#' nr_fit <-
#'   init_biqq(model_code = stan_nr,
#'             data = data_nr_init)
#'
#' # Fit the model for Natural Resource Curse example
#' biqq_nr(fit = nr_fit)
#' }
#'
#' @import rstan parallel
#'
#' @export

biqq_nr <- function(fit         = NULL,
                    XYKK        = c(61,6,64,1,1,3),
                    alpha_prior = rep(1, times = 4),
                    pi_alpha    = rep(1, times = 4),
                    pi_beta     = rep(1, times = 4),
                    q1bd1_alpha = c(1,1),
                    q1bd1_beta  = c(1,1),
                    q2bd1_alpha = c(1,1),
                    q2bd1_beta  = c(1,1),
                    iter        = 20000,
                    chains      = 4,
                    warmup      = 5000,
                    seed        = 1000,
                    cores       = 1) {

  # Load namespaces of required packages
  requireNamespace("rstan", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)

  # Rstan setup options
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  # Check that the model fit is provided
  if (is.null(fit))
    stop("You must provide a model fit object created by the call of rstan::stan(). The model fit for Natural Resource Curse example can be obtained by calling init_biqq(model_code = stan_nr, data = data_nr_init).")

  # Generate the data list
  data <- list(
    XYKK         = XYKK,
    alpha_prior  = alpha_prior,
    pi_alpha     = pi_alpha,
    pi_beta      = pi_beta,
    q1bd1_alpha  = q1bd1_alpha,
    q1bd1_beta   = q1bd1_beta,
    q2bd1_alpha  = q2bd1_alpha,
    q2bd1_beta   = q2bd1_beta
  )

  # Fit the Stan model
  posterior <- suppressMessages(
    rstan::stan(fit     = fit,
                data    = data,
                iter    = iter,
                chains  = chains,
                warmup  = warmup,
                verbose = FALSE,
                seed    = seed,
                cores   = cores,
                refresh = 0))

  return(posterior)
}
