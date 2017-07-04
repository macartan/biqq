#' biqq function, taking data as arguments
#' Basic model; used to define fast biqq function
#'
#' @param fit a fitted model
#' @param XY  XY data -- 4-vector  00 01 10 11
#' @param XYK XYK data -- 12-vector 000 001 010 011...
#' @param XYK_data XYK data provided as a dataframe
#' @param doubly_decisive if TRUE pobative value defaults to clue traces for a and b types only
#' @return A stan object
#' @export
#'
#' @keywords biqq
#' @examples
#' print(biqq(alpha_prior = c(5,1,1,1)))
#'

biqq <- function(fit         = fit_biqq,

                 XY          = c(  0,0,0,0),
                 XYK         = c(  0,  0,
                                   0,  0,
                                   0,  0,
                                   0,  0),
                 XYK_data    = NULL, # option to provide XYK as data frame
                 alpha_prior = c(  1,  1,  1,  1),
                 pi_alpha    = c(  1,  1,  1,  1),
                 pi_beta     = c(  1,  1,  1,  1),
                 q0_alpha    = c(  1,  1,  1,  1),
                 q0_beta     = c(  1,  1,  1,  1),
                 q1_alpha    = c(  1,  1,  1,  1),
                 q1_beta     = c(  1,  1,  1,  1),
                 doubly_decisive = FALSE,
                 iter        = 1000,
                 chains      = 2,
                  ...
)  {

  if(!is.null(XYK_data)) {XY  <- biqq_cases(XYK_data)[1:4]
                          XYX <- biqq_cases(XYK_data)[5:12]}
  if(doubly_decisive){
    q0_alpha  <- c(  1000,  1000,  1,  1)
    q0_beta   <- c(  1,  1,  1000,  1000)
    q1_alpha  <- c(  1000,  1000,  1,  1)
    q1_beta   <- c(  1,  1,  1000,  1000)
  }

  data <- list(
    XY         = XY,
    XYK         = XYK,
    alpha_prior = alpha_prior,
    pi_alpha    = pi_alpha,
    pi_beta     = pi_beta,
    q0_alpha    = q0_alpha,
    q0_beta     = q0_beta,
    q1_alpha    = q1_alpha,
    q1_beta     = q1_beta
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
