
#' Generate permutations
#' @param v a vector
#' @export
#' @examples
#' perm(v = c(2,2,2))
#' perm(v = 3:1)

perm <- function (v) {
  sapply(1:length(v), function(x) {
    rep(rep(1:v[x], each = prod(v[x:length(v)])/v[x]), length.out = prod(v))
  }) - 1
  }


#' Returns a matrix of XY, XYK data
#'
#' Given XY_base and number of clues sought.
#' Standard ordering should be : XY: 00, 01, 10, 11
#' NOTE TO SELF: Currently there are dupicate rows to be removed
#' @param  XY_base is the data before considering K data
#' @param  k number of clues sought
#' @export
#' @examples
#' possible_XYK()
#' possible_XYK(XY_base  = c(1,1,1,0), k=2)
#'
possible_XYK <- function(XY_base  = c(1,1,1,1), k=1){
  n  <- sum(XY_base)
  cs <- combn(1:n,k)
  results <- perm(rep(2, k))
  Ks <-  sapply(1:ncol(cs), function(i) {
    sapply(1:nrow(results), function(j){
      K_sought <- rep(NA,n)
      K_sought[cs[,i]] <- results[j,]
      K_sought
    })}
  )
  Ks <-  matrix(Ks, nrow = n)
  X <- rep(c(0,0,1,1), times = XY_base)
  Y <- rep(c(0,1,0,1), times = XY_base)
  out <- sapply(1:ncol(Ks), function(j) biqq_cases(data.frame(X, Y, K=Ks[,j])))
  out[, !duplicated(t(out))]
}

#' Calculate losses over multiple data possibilities
#' @param  XY_base XY data prior to considering K data
#' @param  k number of clues sought
#' @param  my_loss A loss function
#' @export
#' @examples
#' my_loss <- function(model){
#' D <-extract(model, pars='abcd')
#' return(var(D$abcd[,2]-D$abcd[,1]))
#' }
#' losses(my_loss = my_loss)

losses <- function(XY_base  = c(1,1,1,1), k=1, my_loss, ...){
  D <- possible_XYK(XY_base  = XY_base, k = k)
  XYs  <- D[1:4,]
  XYKs <- D[5:12,]

  out <- sapply(1:ncol(D), function(j) {
    c(XYs[,j], XYKs[,j],
      my_loss(biqq(XY = XYs[,j], XYK = XYKs[,j], ...)))
  })
  out <- t(out)
  colnames(out) <- c("00","01","10","11","000","001","010","011","100","101","110","111", "post var")
  return(out)
  }



#' Assessing likelihood of different clue search results given prior
#' @param model Prior (or posterior, given XY)
#' @param k A 4-vector saying how many clues are to be sought in each case.
#' @return Dataframe showing XYK combinations and the likelihood of seeing so many positive clues in each XY case
#' @export
#' @examples
#' model <- biqq()
#' prob_k_seen(model = model, k = 0:3)

prob_k_seen <- function(model, k = c(1,1,1,1)){
  # In order 00 01 10 11
  # Note, probability of K conditional on being in cell...
  posterior <- rstan::extract(model, permuted = TRUE)
  ps <- with(posterior, {
    sapply(1:4, function(j) {                     # loop through 4 XY cells
      apply(
        sapply(0:k[j], function(i) dbinom(i, k[j],  # binomial weights for k possibilities
                                          w_XYK[,2*j]/w_XY[,j])), 2, mean)            # Expected probability over the runs -- need to check that expectation can be taken here: effectively sum_i p_data var(data)
    })})

  data.frame(
    X = rep(c(0, 0, 1, 1), k+1),
    Y = rep(c(0, 1, 0, 1), k+1),
    K = as.vector(unlist(sapply(k, function(i) 0:i))),
    p = as.vector(unlist(ps)))
}



#' Generate Case Frequencies from Data
#'
#' Function takes data frame and ordering as arguments and returns ordered named vector with frequencies of all unique rows in the data frame
#'
#' @param df Data frame. Consists of all observed binary data
#' @param ordering  Numeric vector. Vectors of ordering for all unique combinations in the output frequencies vector
#' @export

biqq_binary_classes <- function(df,
                                ordering = c(1,5,3,7,2,6,4,8) ) {

  class_types <- expand.grid(as.data.frame(matrix(rep(0:1, ncol(df)), nrow = length(0:1))))
  names(class_types) <- names(df)

  if (length(ordering) != nrow(class_types))
    stop( paste0("Ordering should be provided for all possible combinations of ",
                 paste(names(df), collapse = ","),
                 ". This implies that the length of ordering vector should be ",
                 nrow(class_types), ".") )

  if (!all(sort(ordering) == 1:nrow(class_types)))
    stop( paste0("Ordering include the following numbers (each one time): ",
                 paste(1:nrow(class_types), collapse = ","), ".") )

  counts <- stats::aggregate(freq ~ ., data = base::transform(df, freq = 1), FUN = length)
  classes <- base::merge(class_types, counts, by = colnames(class_types), all.x = TRUE)
  classes$freq[is.na(classes$freq)] <- 0

  rownames(classes) <-
    apply(classes, MARGIN = 1,
          FUN = function(x) paste0("N", paste0(x[colnames(class_types)], collapse = "") ) )

  return(apply(classes[,"freq", drop = FALSE],
               MARGIN = 1, FUN = as.numeric)[ordering])
}
