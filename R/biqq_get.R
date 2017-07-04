#' biqq_get
#'
#' Asks for a quantity from a model
#' @param model a fitted model
#' @param Q  quantity of interest expressed as function. Defaults to ATE. Can be any function of data. Alternatively some key words predefined for some quantities of interest, eg ATE, phi.
#' @param digits digits returned
#' @export
#'
#' @keywords biqq
#' @examples
#' model <- biqq()
#' biqq_get(model)
#' biqq_get(model, Q = function(d) c(lambda_b  = mean(d$abcd.2)))

biqq_get <- function(model         = fit_biqq,
                     Q = NULL,
                     digits = 2,
                     ...) {
  if(is.function(Q)) {f <- Q}
  else {
  if(is.null(Q))      Q <- "ATE"
  if(Q == "ATE")      f <- function(d) with(d, c(ate = mean(abcd.2 - abcd.1), var_ate = var(abcd.2 - abcd.1)))
  if(Q == "phi_b1d1") f <- function(d) with(d, c(phi_b_1  = mean(d$phi_b_1),  phi_d_1  = mean(d$phi_d_1)))
  if(Q == "phi")      f <- function(d) with(d, apply(dplyr::select(d, starts_with("q")), 2, mean))
  if(Q == "f")        f <- function(d) with(d, apply(dplyr::select(d, starts_with("f")), 2, mean))}

  out <- f(as.data.frame(extract(model)))
  round(out, digits = digits)
  }
