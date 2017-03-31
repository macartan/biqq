#' BIQQ Model for Natural Resource Curse Example
#'
#' Fit BIQQ Model for Natural Resource Curse Example
#'
#' @param fit Stan model fit object
#' @param title Character string of graph name. Defaults to \code{""}
#' @param x_fun Function which takes extracted posterior as an argument and returns values plotted on x axis
#' @param x_name Character string of x axis statistic name
#' @param y_fun Function which takes extracted posterior as an argument and returns values plotted on y axis
#' @param y_name Character string of y axis statistic name
#' @param object_colors Character vector of length 2. Colors to be used for plotting of posterior draws and scaled density respectively. Defaults to \code{object_colors = c("darkcyan", "darkred")}
#' @param adjust_value Scalar. Adjustment of the smoothing bandwidth, the bandwidth used is \code{adjust*bw}, where \code{bw} is a result of "nrd0" method
#'
#' @keywords stan biqq plot
#'
#' @examples
#'
#' \dontrun{
#'  nr_fit <- init_biqq(model_code = stan_nr, data = data_nr_init)
#'  nr_example <- biqq_nr(fit = nr_fit, chains = 4, cores = 1)
#'
#'  plot_biqq(fit = nr_example, title = "Posteriors | Both clues weak")
#' }
#'
#' @import rstan ggplot2
#'
#' @export


plot_biqq <- function(fit,
                      title = "",
                      x_fun = function(x) { x$abcd[,2] - x$abcd[,1] },
                      x_name = "ATE",
                      y_fun = function(x) { x$pi_c },
                      y_name = "Assignment propensity for c types",
                      object_colors = c("darkcyan", "darkred"),
                      adjust_value = .8) {

  requireNamespace("rstan", quietly = TRUE)
  requireNamespace("ggplot2", quietly = TRUE)

  if (class(fit) != "stanfit")
    stop("Fit should be an object of class stanfit, as produced by the call of stan(), biqq_es() or biqq_nr().")

  posterior <- rstan::extract(fit)

  df <- data.frame(x_post = x_fun(posterior), y_post = y_fun(posterior))

  p <-
    ggplot(df) +
    geom_point(aes_string(x = "x_post", y = "y_post", colour = '"Posterior draws"'), alpha = 0.1, size = .05) +
    geom_line(aes_string(x = "x_post", y = ".5 * ..scaled..", col = '"Scaled density"'), stat = "density",
              adjust = adjust_value, alpha = 1/2, size = .6) +
    theme_linedraw(base_size = 12, base_family = "Helvetica") +
    theme(legend.position = "top") +
    labs(title = title,
         x = paste0(x_name, " (posterior mean = ", round(mean(x_fun(posterior)),3), ")"),
         y = y_name) +
    scale_color_manual(values = object_colors) +
    guides(colour =
             guide_legend(title = NULL,
                          override.aes = list(shape = c(16,NA),
                                              linetype = c(0,1),
                                              size = c(2,1),
                                              alpha = c(1,1))))

  return(p)

}
