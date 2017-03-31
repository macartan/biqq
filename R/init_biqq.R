#' Initiation of Stan Program for BIQQ
#'
#' The function takes the model code and data as an argument and returns stan_fit object to use in BIQQ
#'
#' @param model_code A character string either containing the model definition or the name of a character string object in the workspace. This parameter is used only if parameter file is not specified. When fit is specified, the model compiled previously is used so specifying model_code is ignored
#' @param stanmodel An object of class \code{stanmodel}
#' @param data A named list or environment providing the data for the model, or a character vector for all the names of objects used as data. See the Note section below
#' @keywords stan init example
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
#' # Print model fit for Natural Resource Curse example
#' nr_fit
#'
#' # Create model fit for Origins of Electoral Systems example using init_biqq function
#' es_fit <-
#'   init_biqq(model_code = stan_es,
#'             data = data_es_init)
#'
#' # Print model fit for Origins of Electoral Systems example
#' es_fit
#' }
#'
#' @import rstan
#'
#' @export

init_biqq <- function(model_code = NULL,
                      stanmodel = NULL,
                      data) {

  requireNamespace("rstan", quietly = TRUE)

  if (!is.null(model_code)) {
  suppressMessages(
    rstan::stan(model_code = model_code,
                data       = data,
                iter       = 1,
                chains     = 1,
                cores      = 1,
                verbose    = FALSE)
  )
  } else if (!is.null(stanmodel)) {
    suppressMessages(
      rstan::sampling(object     = stanmodel,
                      data       = data,
                      iter       = 1,
                      chains     = 1,
                      cores      = 1,
                      verbose    = FALSE)
    )
  }
}

