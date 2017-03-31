
#' Function to count the set in each permutation of X,Y,K (Data in that order) and output summary
#' Consisent with this naming ordering from Stan Code:
#'  int<lower=0> XY[4];         // summary data, providing number of XY: order: 00, 01, 10, 11
#'  int<lower=0> XYK[8];        // summary data, providing number of XYK: order 000, 001, 010,...
#' @param D   a dataframe or matrix with three columns, X, Y, K
#' @keywords Conversion biqq
#'
#' @examples biqq_cases(D = data.frame(X = c(0,1,1), Y = c(1,1,0), K = c(NA, NA, 1)))

biqq_cases = function(D){
  n <- nrow(D)
  D[is.na(D)] <- 3
  count.class <- function(s)  sum(colSums(sapply(1:n, function(i) D[i,]==c(s)))==3)
  class.types <- cbind(X=c(rep(0:1, each = 2),rep(0:1, each = 4)),
                       Y=c(rep(0:1,2), rep(rep(0:1, each = 2),2)),
                       K=c(rep(3, 4), rep(0:1,4)))
  classes = apply(class.types,1,count.class)
  names(classes) = c("00x", "01x", "10x","11x", "000", "001", "010", "011", "100", "101", "110", "111")
  return(classes)
}

