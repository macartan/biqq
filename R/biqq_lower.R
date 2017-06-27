#' biqq from lower level model, taking data as arguments
#'
#' Prior information is provided for Y_u_prior (16 elements) and K_u_prior (4 elements)
#'
#' @param fit a fitted biqq_lower model
#' @param XY  XY data -- 4-vector  00 01 10 11
#' @param XYK XYK data -- 12-vector 000 001 010 011...
#' @param XYK_data XYK data provided as a dataframe
#' @param Y_u_prior Dirichlet hyperparameters on 16 u_Y types
#' @param K_u_prior Dirichlet hyperparameters on 4 u_K types
#' @param pi_alpha 16 alpha Beta pramaters for pis
#' @param pi_alpha 16 beta Beta pramaters for pis
#' @param doubly_decisive if TRUE pobative value defaults to clue traces  for a and b types (not extreme 0/1)
#' @return A stan object
#' @export
#'
#' @keywords biqq dag
#' @examples
#' print(biqq_lower(XYK = c(0, 1, 1, 0, 1, 0, 0, 1)))
#' print(biqq_lower(XYK = c(0, 1, 1, 0, 1, 0, 0, 1), doubly_decisive = TRUE))
#'

biqq_lower <- function(fit         = fit_biqq_lower,
                 XY          = c(  0,  0,  0, 0),
                 XYK         = c(  0,  0,
                                   0,  0,
                                   0,  0,
                                   0,  0),
                 XYK_data    = NULL, # option to provide XYK as data frame
                 Y_u_prior       = rep(.5,  16),
                 K_u_prior       = rep(.5,4),
                 pi_alpha    = matrix(.5, 16, 4),
                 pi_beta     = matrix(.5, 16, 4),
                 doubly_decisive = FALSE,
                 iter        = 1000,
                 chains      = 2,
                  ...
)  {

  if(!is.null(XYK_data)) {XY  <- biqq_cases(XYK_data)[1:4]
                          XYX <- biqq_cases(XYK_data)[5:12]}
  if(doubly_decisive){
    Y_u_prior       = rep(.1,  16)
    Y_u_prior[c(2,5,12,15)] <- 2
    K_u_prior       =  c(2,.1,.1,2)
    }

  data <-   list(
    XY         = XY,
    XYK         = XYK,
    XYK_data    = NULL, # option to provide XYK as data frame
    Y_u_prior   = Y_u_prior,
    K_u_prior   = K_u_prior,
    pi_alpha    = pi_alpha,
    pi_beta     = pi_beta
  )
  #  suppressWarnings can also be added to reduce print out
  posterior <- suppressMessages(
    stan(fit     = fit,
         data    = data,
         iter    = iter,
         chains  = chains)
  )

  return(posterior)
}
