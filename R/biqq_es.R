#' BIQQ Model for Origins of Electoral Systems Example
#'
#' Fit BIQQ Model for Origins of Electoral Systems Example
#'
#' @param fit Stan model fit object
#' @param XY XY data for model fitting (clues not sought)
#' @param XYK XYK data for model fitting (clues sought)
#' @param alpha_prior Numeric vector of length 4. Dirichlet distribution parameters for proportions of 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param q0_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of not observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q0_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of not observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q1_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q1_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param chains Integer. Number of MC chains. Defaults to 2
#' @param iter Integer. Total number of iterations in each chain. Defaults to 20000,
#' @param warmup Integer. Number of warm-up iterations in each chain. Defaults to 5000
#' @param cores Integer. Number of cores to use for parallel computation. Defaults to 1
#'
#' @keywords stan biqq example electoral systems
#'
#' @examples
#'
#' \dontrun{
#'
#' library(biqq)
#'
#' # Create model fit for Origins of Electoral Systems example using init_biqq function
#' es_fit <-
#'   init_biqq(model_code = stan_es,
#'             data = data_es_init)
#'
#' # Fit the model for Origins of Electoral Systems example
#' biqq_es(fit = es_fit)
#' }
#'
#' @import rstan parallel
#'
#' @export

biqq_es <- function(fit         = NULL,
                    XY          = rep(1, times = 4),
                    XYK         = rep(1, times = 8),
                    alpha_prior = rep(1, times = 4),
                    pi_alpha    = rep(1, times = 4),
                    pi_beta     = rep(1, times = 4),
                    q0_alpha    = rep(1, times = 4),
                    q0_beta     = rep(1, times = 4),
                    q1_alpha    = rep(1, times = 4),
                    q1_beta     = rep(1, times = 4),
                    iter        = 20000,
                    chains      = 2,
                    cores       = 1,
                    warmup      = 5000) {

  # Load namespaces of required packages
  requireNamespace("rstan", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)

  # Rstan setup options
  rstan_options(auto_write = TRUE)
  options(mc.cores = parallel::detectCores())

  # Check that the model fit is provided
  if (is.null(fit))
    stop("You must provide a model fit object created by the call of rstan::stan(). The model fit for Electoral Systems example can be obtained by calling init_biqq(model_code = stan_es, data = data_es_init).")

  # Generate the data list
  data <-
    list(XY          = XY,
         XYK         = XYK,
         alpha_prior = alpha_prior,
         pi_alpha    = pi_alpha,
         pi_beta     = pi_beta,
         q0_alpha    = q0_alpha,
         q0_beta     = q0_beta,
         q1_alpha    = q1_alpha,
         q1_beta     = q1_beta
    )

  # Fit the Stan model
  posterior <-
    suppressMessages(
      stan(fit     = fit,
           data    = data,
           iter    = iter,
           chains  = chains,
           warmup  = warmup,
           cores   = cores,
           verbose = FALSE,
           refresh = 0)
    )

  return(posterior)
}
