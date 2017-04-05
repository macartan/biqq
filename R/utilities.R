
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


#' Calculate a statistic givrn some model and a function f.
#'
#' Defaults to the posterior variance on lambda_b - lambda_a for simple biqq model
#' @param model a fitted biqq model
#' @param L a loss function, defined over the dataframe extracted from model L. Defaults to posterior variance of ATE
#' @export
#' @examples
#' model <- biqq()
#' my_loss(model, f = function() {var(abcd[,2])})

my_loss <- function(model,
                    L = function(D) {var(D$abcd[,2]-D$abcd[,1])}
            ){
              D <-extract(model, pars='abcd')
              out <- L(D)
              return(out)
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
#' @param  f A loss function
#' @param  strategies if only a subset of strategies are examined strategies can be indicated  of the form "0-0-1-1" meaning seek one clue for 10 cases and one for 11 cases, or more generally c("1-1-0-0", "0-0-1-1"). Use this option with caution to be sure that naming of strategies is correct and unique
#' @export
#' @examples
#' my_loss <- function(model){
#' D <-extract(model, pars='abcd')
#' return(var(D$abcd[,2]-D$abcd[,1]))
#' }
#' losses(f = my_loss)
#' losses(strategies = c("1-0-0-0","0-0-0-1"))
#' losses(strategies = c("1-0-0-1"))
#' losses(strategies = c("1-0-0-1"), k = 2)
#'
losses <- function(XY_base  = c(1,1,1,1), k=1,  f = my_loss, strategies = NULL, ...){

  D <- possible_XYK(XY_base  = XY_base, k = k)
  x <- t(D)
  s <- sapply(1:nrow(x), function(j) {
    +     c(x[j,5]+x[j,6], x[j,7]+x[j,8], x[j,9]+x[j,10], x[j,11]+x[j,12])})
  s <- apply(s, 2, paste, collapse = "-")
  colnames(D) <- s

  if(!is.null(strategies)){
  if(sum(s %in% strategies)==0)  stop("Strategy not found in list of possible strategies. Check assumptions on base XY data or value of k.")
  D <- D[, s %in% strategies]
  }

  XYs  <- D[1:4,]
  XYKs <- D[5:12,]

  out <- sapply(1:ncol(D), function(j) {
    c(XYs[,j], XYKs[,j],
      f(biqq(XY = XYs[,j], XYK = XYKs[,j], ...)))
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




#' Function to identify possible case selection strategies and assess expected gains
#'
#' Assesses the set of strategies that could be implemented by looking for clues in  k cases given existing X,Y data, recorded in a 4-vector, `XY_base`. The output lists the possible strategies, assesses the probability of different data realizations that coudl result from each strategy, the posterior variance that arises in each case, and the overall expected posterior variance from each strategy under consideration.
#' @param k A scalar saying how many clues are to be sought in total
#' @param  XY_base is the data before considering K data
#' @param  strategies optional: if only a subset of strategies are examined strategies can be indicated  of the form "0-0-1-1" meaning seek one clue for 10 cases and one for 11 cases, or more generally c("1-1-0-0", "0-0-1-1"). Use this option with caution to be sure that naming of strategies is correct and unique
#' @return A list
#' @export
#' @examples
#' Expected_Var(k = 1, XY_base = rep(1,4), chains = 2)
#' Expected_Var(k = 1, XY_base = rep(1,4), chains = 2, strategies = c("1-0-0-1"))
#' Expected_Var(k = 2, XY_base = rep(1,4), chains = 2, strategies = c("1-0-0-1"))

Expected_Var <- function(k = 1, XY_base = rep(1,4), strategies = NULL, f = my_loss, name = "Results",  ...){

  model <- biqq(XY = XY_base, ...)
  D <-extract(model, pars='abcd')
  base_ate <-  {mean(D$abcd[,2]-D$abcd[,1])}
  base_loss <- f(model)


  x <- losses(k=k, XY_base = XY_base, strategies = strategies, f=f, ...)

  # Probability of a given outcome given the strategy
  g <- function(j,
                m = model,
                strategy =  c(x[j,5]+x[j,6], x[j,7]+x[j,8], x[j,9]+x[j,10], x[j,11]+x[j,12])) {
    p <- prob_k_seen(model = m, k = strategy)
    prod(
      p[p$X==0 & p$Y==0,4][(x[j,6]+1)],
      p[p$X==0 & p$Y==1,4][(x[j,8]+1)],
      p[p$X==1 & p$Y==0,4][(x[j,10]+1)],
      p[p$X==1 & p$Y==1,4][(x[j,12]+1)])}


    strat_names <- sapply(1:nrow(x), function(j) {
      c(x[j,5]+x[j,6], x[j,7]+x[j,8], x[j,9]+x[j,10], x[j,11]+x[j,12])})
    strat_names <- apply(strat_names, 2, paste, collapse = "-")

  # Probability of each outcome for each strategy
  p_each <- sapply(unique(strat_names), function(ss) {
    sapply((1:nrow(x))[strat_names==ss], g)
  })

  rnames <- paste("When", 0:(nrow(p_each)-1), "clues found: ")
  if(length(rnames)>1) rnames[2] <- "When 1 clue found: "
  rownames(p_each) <- rnames

  # Expected posterior variance associate with each strategy
  var_each <- sapply(unique(strat_names), function(s) {x[strat_names==s, 13]})
  rownames(var_each) <- rnames

  # Expected posterior variance associate with each strategy
  Ex_Post_var <- sapply(unique(strat_names), function(s) {
    sapply((1:nrow(x))[strat_names==s], g)%*%x[strat_names==s, 13]
  })
  return(list(name, beliefs_from_XY_only =
                c(base_ate = base_ate, base_loss = base_loss),
              strategies_examined = unique(strat_names),
              probability_each_pattern = p_each, loss_for_each_pattern = var_each, Expected_loss_from_strategy = Ex_Post_var))
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
