#' BIQQ Model for k Clues for Origins of Electoral Systems Example
#'
#' Fit BIQQ Model for Origins of Electoral Systems example assuming that only k clues are sought
#'
#' @param k Integer. Number of clues sought
#' @param clue_data Data frame. All known XYK cases from which only k cases are randomly drawn for BIQQ model
#' @param fit Stan model fit object
#' @param alpha_prior Numeric vector of length 4. Dirichlet distribution parameters for proportions of 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param pi_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of assignment for 4 possible types in the population. Defaults to \code{c(1,1,1,1)}
#' @param q0_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of not observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q0_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of not observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q1_alpha Numeric vector of length 4. Alpha shape parameters for Beta distribution of probabilities of observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param q1_beta Numeric vector of length 4. Beta shape parameters for Beta distribution of probabilities of observing clue given that it was sought and any of four possible types. Defaults to \code{c(1,1,1,1)}
#' @param extract_pars Character vector. Names of posterior parameters to extract as taken by \code{\link{extract}} function. These should be sufficient to compute output of interest for BIQQ program using \code{out_fun}. Defaults to \code{extract_pars = "abcd"}, which implies a call \code{rstan::extract(posterior, pars = "abcd")}
#' @param out_fun Function. The function, which takes posterior samples extracted using \code{extract_pars} as an argument and returns quantity/ies of interest. Defaults to \code{out_fun = function(x) { mean(x$abcd[,2] - x$abcd[,1]) } }, which returns mean of posterior difference between proportions of beneficial and adverse types in the population. This function is used in the replication of Kreuzer study presented in Humphreys, Jacobs (2016)
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
#' es_fit <-
#'   init_biqq(model_code = stan_es,
#'             data = data_es_init)
#'
#' sims <- 10 # Number of simulations
#' N <- nrow(data_es) # Number of cases
#'
#' # Generate hyperparameters for Stan
#' q0_alpha_beta <- mapply(m = c(q.a0 = 0.1, # Assumptions on mean of q0's
#'                               q.b0 = 0.1,
#'                               q.c0 = 0.05,
#'                               q.d0 = 0.3),
#'                         sd = rep(.01, times = 4),
#'                         FUN = beta_prior)
#'
#' q1_alpha_beta <- mapply(m = c(q.a1 = 0.95, # Assumptions on mean of q1's
#'                               q.b1 = 0.9,
#'                               q.c1 = 0.475,
#'                               q.d1 = 0.5),
#'                         sd = rep(.01, times = 4),
#'                         FUN = beta_prior)
#' # Rstan setup options
#' rstan_options(auto_write = TRUE)
#' options(mc.cores = parallel::detectCores())
#'
#' # Run the analysis
#' betas <-
#'   parallel::mclapply(X          = rep(0:N, each = sims),
#'                      FUN        = biqq::biqq_es_k,
#'                      fit        = es_fit,
#'                      clue_data  = data_es,
#'                      q0_alpha   = q0_alpha_beta["alpha",],
#'                      q0_beta    = q0_alpha_beta["beta",],
#'                      q1_alpha   = q1_alpha_beta["alpha",],
#'                      q1_beta    = q1_alpha_beta["beta",],
#'                      chains     = 2,
#'                      cores      = 1,
#'                      extract_pars = "abcd",
#'                      out_fun = function(x) { mean(x$abcd[,2] - x$abcd[,1]) })
#'
#' betas <- matrix(unlist(betas), ncol = sims, byrow = TRUE)
#' }
#'
#' @import rstan parallel
#'
#' @export


biqq_es_k <-
  function(k, clue_data,
           fit          = NULL,
           alpha_prior  = rep(1, times = 4),
           pi_alpha     = rep(1, times = 4),
           pi_beta      = rep(1, times = 4),
           q0_alpha     = rep(1, times = 4),
           q0_beta      = rep(1, times = 4),
           q1_alpha     = rep(1, times = 4),
           q1_beta      = rep(1, times = 4),
           extract_pars = "abcd",
           out_fun      = function(x) { mean(x$abcd[,2] - x$abcd[,1]) },
           iter         = 8000,
           chains       = 2,
           warmup       = 1000,
           cores        = 1
           ) {

  # Load namespaces of required packages
  requireNamespace("rstan", quietly = TRUE)
  requireNamespace("parallel", quietly = TRUE)

  # Check that the model fit is provided
  if (is.null(fit))
    stop("You must provide a model fit object created by the call of rstan::stan(). The model fit for Electoral Systems example can be obtained by calling init_biqq(model_code = stan_es, data = data_es_init).")

  N <- nrow(clue_data)

  if (k > N) stop("Number of cases studied (k) should be less or equal to total number of cases avaliable in clue_data")

  clues <- sample(x = N, size = k) # sample clues

  if (k == 0) {
    XYK   <- rep(0, times = 8)
  } else {
    XYK   <- unname( biqq_binary_classes(clue_data[clues,],
                                         ordering = c(1,5,3,7,2,6,4,8)) )
  }
  if (k == 20) {
    XY <- rep(0, times = 4)
  } else {
    XY <- unname( biqq_binary_classes(clue_data[!(1:N %in% clues), c("X","Y")],
                                      ordering = c(1,3,2,4)) )
  }

  posterior <-
    biqq_es(fit         = fit,
            XY          = XY,
            XYK         = XYK,
            alpha_prior = alpha_prior,
            pi_alpha    = pi_alpha,
            pi_beta     = pi_beta,
            q0_alpha    = q0_alpha,
            q0_beta     = q0_beta,
            q1_alpha    = q1_alpha,
            q1_beta     = q1_beta,
            iter        = iter,
            chains      = chains,
            warmup      = warmup,
            cores       = cores)

  posterior <- rstan::extract(posterior, pars = extract_pars)
  return( out_fun(posterior) )
}
