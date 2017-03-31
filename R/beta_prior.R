#' Re-parametrize Mean and Standard Deviation parameters for Beta Distribution into Shape Parameters
#'
#' Function takes mean and standard deviation as an argument and returns conventional alpha and beta shape parameters
#'
#' @param m Numeric. Mean of the Beta distribution
#' @param sd Positive numeric. Standard deviation of the Beta distribution
#'
#' @keywords beta distribution
#'
#' @examples
#'
#' \dontrun{
#' # Generate hyperparameters for Stan model for Origins of Electoral Systems example
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
#'
#' }
#'
#' @export

# Function for translating mean and sd into shape hyperparameters (alpha and beta)
beta_prior <- function(m, sd) {
  # Check assumption on variance
  if ( (sd^2 < m * (1 - m)) == F ) stop("The variance must be less than (m * (1 - m)).")

  alpha <- -m * (sd^2 + m^2 - m)/sd^2
  beta  <- (m - 1) * (sd^2 + m^2 - m)/sd^2

  return(c(alpha = alpha, beta = beta))
}
