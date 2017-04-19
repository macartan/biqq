
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


#' Declare a structural model
#'
#' @param var_functions A list of structural equations
#' @param var_names An optional list of node names
#' @param P An optional list of node names
#' @export
#' @examples
#' biqq_model()

biqq_model <- function(
  var_functions =
    list( f_S = function(U,V) U[1] < .5,
          f_X = function(U,V) U[2] < .5,
          f_C = function(U,V) (1- (V[1]*V[2]))*(U[3]>.2) +(U[3]<.1),
          f_R = function(U,V) V[2]*V[3]*(U[4]>.2) + (U[4]<.1),
          f_Y = function(U,V) V[3]*V[4]*(U[5]>.2)+ (U[5]<.1)
    ),
  var_names = 1:length(var_functions),
  P =  function() runif(length(var_functions))){
  model <- list(var_names = var_names, var_functions = var_functions, P =  P)
  return(model)
}


#' Function to generate a world from a structural model
#'
#' @param model A model, made using biqq_model
#' @param U  A context. Generally generated using model$P()
#' @param do an optional  list of do operations on nodes
#' @export
#' @examples
#' M <- biqq_model()
#' biqq_world(M, M$P())

biqq_world <- function(model, U, do = NULL){
  n <- length(model$var_names)
  if(!is.null(do)) for(i in 1:n){ if(!is.na(do[i])) {model$var_functions[[i]] <- function(U,V)  do[i]}}
  V <- rep(NA,n)
  for(i in 1:n) V[i] <- model$var_functions[[i]](U,V)
  names(V) = model$var_names
  unlist(V)
}



#' A function to figure out which worlds satisfy some query
#'
#' @param model A model, made using biqq_model
#' @param operations  A set of operations
#' @param query_v A query defined on outcomes on observables, V
#' @param sims optional number of simulations for draws
#' @param U optional matrix of contexts (worlds)
#' @param plotit Plot the results of the investigation as matrix of 2 way plots (pairs). Defaults to FALSE
#' @export
#'
#' @examples
#' M <- biqq_model()
#' biqq_world(M, M$P())
#' biqq_which(
#'    model,
#'    operations = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA)),
#'    query_v = function(x) (x[[1]][5] ==1) & (x[[2]][5] == 0),
#'    sims = 500,
#'    plotit = TRUE)
#'

biqq_which <- function(
  model,
  operations,
  query_v,
  sims = NULL,
  U = NULL, # Optionally provide a matrix of worlds: Useful for keeping U fixed under different experiments. Each U should be a single number, not a vector.
  plotit = FALSE,
  col1 = "red",
  col2 = "blue") {

  nU <- length(model$P())
  if(is.null(U)) U <- (replicate(sims, model$P()))
  sims <- ncol(U)
  V <-  matrix(NA, sims, nU)
  A  <- rep(NA, sims)
  k  <- length(operations)

  for(j in 1:sims){
    u  <- U[,j]
    v  <- biqq_world(model, u)
    x  <- list()
    a  <- rep(NA, k)
    for(i in 1:k){
      x[[i]]  <- biqq_world(model, u, do = operations[[i]])
    }
    A[j]  <- query_v(x)
    V[j,] <- v
  }

  V <- as.data.frame(V)
  colnames(V) <- model$var_names

  result <- list(U=U, V=V, A=A)

  if(plotit) { if(dim(result[[1]])[1]>2){stop("graphing currently only possible for U shocks on 2 vars")
    } else {
      pairs(t(result[[1]]), col = ifelse(result[[3]], col1, col2), labels = model$var_names, pch = ".")
 }}
  return(result)
}



#' Calculate what is learned from a set of observations, Ks, given existing data, W
#'
#' Defaults to the posterior variance on lambda_b - lambda_a for simple biqq model
#' @param model a fitted biqq model
#' @param f any function; defaults to f=var, posterior variance. For posterior mean set f = mean
#' @param Ks  A matrix with one column per node indicating whther to be sought or not. Defaults to NULL
#' @export
#' @examples
#' model <- biqq_model()
#' operations = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA))
#' query_v = function(x) (x[[1]][5] ==1) & (x[[2]][5] == 0)
#' biqq_learning(model,
#'               operations,
#'               query_v,
#'               Ks = rbind(
#'                      c(TRUE, FALSE, FALSE, FALSE, FALSE),
#'                      c(FALSE, TRUE, FALSE, FALSE, FALSE),
#'                      c(TRUE, TRUE, FALSE, FALSE, FALSE),
#'                      c(TRUE, TRUE, FALSE, TRUE, FALSE),
#'                      c(TRUE, TRUE, TRUE, TRUE, TRUE)
#'          )
#' )

biqq_learning <- function(
  model,
  operations,
  query_v,
  sims=100,
  U=NULL,
  W  = rep(NA, length(model$var_names)),
  Ks = NULL,  # Or a matrix with one column per node indicating whther to be sought or not
  f=var){
    if(!is.null(Ks) & is.null(dim(Ks))) Ks <- matrix(Ks, 1)
    if(!is.null(Ks) & !is.matrix(Ks)) stop("Ks should be provided as a matrix with as many columns as elements in V")
    if(is.null(U)) U <- (replicate(sims, model$P()))
    result <- biqq_which(
      model,
      operations = operations,
      query_v = query_v,
      U=U,
      plotit = FALSE)

    inW <- sapply(1:sims, function(i) !(0 %in% (1*(result$V[i,] == W))))
    # Two cases: when W is not possible and when it is:
    if(sum(inW) ==0 ) {      print("W is 0 probability event")
      out <- c(f(result$A))
      if(!is.null(Ks)) out <- c(out, rep(-1, 1+nrow(Ks)))

    } else {   # Possible W

      out <- c(f(result$A), f(result$A[inW]))

      epv <- sapply(1:nrow(Ks), function(k) {
        K <- Ks[k,]
        Kranges <- apply(as.matrix(result$V[,K]), 2, function(j) length(unique(j)))  # the set of clues

        if(!is.null(Ks) & length(Kranges)>0){
          Ks <- perm(Kranges)
          Kp <- rep(NA, length(K))
          Kp <- t(sapply(1: nrow(Ks), function(i)  {Kp[K] <- Ks[i,]; return(Kp)}))  # Pattern of Ks


          # Expected variance over all the combinations of K
          pE_k <- sapply(1:nrow(Kp), function(j) {
            Kv  <- Kp[j,]
            inK <- sapply(1:sims, function(i) !(0 %in% (1*(result$V[i,] == Kv))))
            pKv <- mean(inK[inW])   # probability of K pattern, given W

            V_given_WK <- f(result$A[inW & inK]) # 0 var for 0 prob event
            if(sum(inW & inK)<=1)   V_given_WK <- 0
            prior_V_given_W <- f(result$A)
            c(pKv, V_given_WK)
          })

          Epostv <- sum(pE_k[1,]*pE_k[2,])
        } else {
          Epostv <- out[2]
        }
      })
      out <- c(out, epv)
    }

    if(!is.null(Ks)) {names(out) <- c("Prior f", "f |W",
                                      sapply(1:nrow(Ks), function(r) {paste(1*Ks[r,], collapse = "")} ))}

    out
  }



#' Calculate a statistic given some model and a function f.
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


  x <- biqq:::losses(k=k, XY_base = XY_base, strategies = strategies, f=f, ...)

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
    pp <- sapply((1:nrow(x))[strat_names==ss], g)
    names(pp) <- apply(x[strat_names==ss, c(6,8, 10, 12)], 1, paste, collapse = "")
    pp  })



  # Expected posterior variance associate with each strategy
  var_each <- sapply(unique(strat_names), function(ss) {
    v <- x[strat_names==ss, 13]
    names(v) <- apply(x[strat_names==ss, c(6,8, 10, 12)], 1, paste, collapse = "")
    v
    })


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
