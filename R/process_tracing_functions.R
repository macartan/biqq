

#' Declare a structural model
#'
#' Declare a set of exogenous and endogenous variables and a set of functions,
#' \code{var_functions}, indicating how variables are affected by each other.
#'
#' @param var_functions a list of structural equations. Each element in the list must be a function specifying how each variable relates to its direct causes. Causes's functions should be specified before effects' functions i.e. if \code{X} causes \code{Y} then \code{f_X} should be specified before \code{f_Y}.
#' @param var_names an optional list of node names.
#' @param P a function to generate a context. The joint probability distribution over the exogenous variables, \code{U}.
#' @return a list containing the declared model.
#' @export
#' @examples
#' pt_model()

pt_model <- function(
  var_functions =
    list( f_S = function(U,V) U[1] < .5,
          f_X = function(U,V) U[2] < .5,
          f_C = function(U,V) (1- (V[1]*V[2]))*(U[3]>.2) +(U[3]<.1),
          f_R = function(U,V) V[2]*V[3]*(U[4]>.2) + (U[4]<.1),
          f_Y = function(U,V) V[3]*V[4]*(U[5]>.2)+ (U[5]<.1)
    ),
  var_names = c("S", "X","C", "R", "Y"),
  P =  function() runif(length(var_functions))){
  model <- list(var_names = var_names, var_functions = var_functions, P =  P)
  return(model)
}


#' Function to generate a world (values for V) from a structural model
#'
#'
#' @param model a model made using  \code{\link{pt_model}}.
#' @param U  a context. A realization of the exogenous variables in the structural model. Generally generated using \code{model$P()}.
#' @param do an optional list of do operations on nodes. An intervention.
#' @return a world.
#' @export
#' @examples
#' M <- pt_model()
#' pt_world(M)
#' pt_world(M, U = runif(5))

pt_world <- function(model, U = NULL, do = NULL){

  # Get U if needed
  if(is.null(U)){U <- model$P()}

  # Tidy U if needed
  if(!is.data.frame(U)) U <- data.frame(t(U))
  names(U) <- paste0("U_",model$var_names)

  # Add V
  n <- length(model$var_names)
  V <- rep(NA, n)

  if(!is.null(do)) for(i in 1:n){if(!is.na(do[i])) {
    model$var_functions[[i]] <- function(U,V)  do[[i]]}}

  for(i in 1:n) V[[i]] <- model$var_functions[[i]](U,V)

  V <- data.frame(t(V))
  names(V) <- model$var_names

  cbind(V,U)
  }

#' Function to generate a dataset of possible worlds from a structural model with V provided
#'
#'
#' @param model a model made using  \code{\link{pt_model}}.
#' @param N Number of simulations
#' @param U a sim*n_nodes dataframe .
#' @param do an optional list of do operations on nodes. An intervention.
#' @return a dataframe with worlds
#' @export
#' @examples
#' M <- pt_model()
#' pt_worlds(M, N = 4)
#' pt_worlds(M, N = 4)
#' U = data.frame(matrix(runif(10), 2, 5))
#' pt_worlds(M, U = U)

pt_worlds <- function(model, N = NULL, do = NULL, U = NULL){
  if(!is.null(U) & !is.null(N)) if(nrow(U)!=N)  stop("U should have N rows")
  if(is.null(U) & is.null(N))  stop("Please provide N or U")
  if(is.null(U) & !is.null(N)) U <- rep(NULL, N)
  if(is.null(N) & !is.null(U)) N <- nrow(U)

  x <- pt_world(model, do = do, U = U[1,])
  if(N>1) {for(j in 2:N) x <- rbind(x, pt_world(model, do = do, U = U[j,]))}
  x
  }


#' A function to figure out which worlds satisfy some query
#'
#' The output is a list that contains a matrix of contexts, U, the outcomes associated with them, V, and an indicator for whether the query is satisfied or not in each context, A.
#'
#' @param model A model made using \code{\link{pt_model}}.
#' @param U     A matrix of contexts (worlds). If missing this is generated based on sims.
#' @param ops  A set of operations.
#' @param query A query defined on outcomes on observables, V.
#' @param sims Optional number of simulations for draws; defaults to 500 if U not defined.
#'
#' @export
#'
#' @examples
#' pt_which(model = pt_model(),
#'          ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA)),
#'          query = function(x) (x[[1]][5] == 1) & (x[[2]][5] == 0),
#'          sims = 3)
#'

pt_which <- function(
  model,                        # the model
  U = NULL,                     # dimensionality: sims * n_nodes
  ops,
  query,
  sims = 500) {

  # U
  if(is.null(U)) {
    U <- data.frame(t(replicate(sims, model$P())))
    names(U) <- paste0("U_", model$var_names)
  }

  # V
  V  <- pt_worlds(model, U = U)[model$var_names]

  # A
  k <- length(ops)

  f <- function(u) {                 # For each context, apply operations, and take result
    x  <- list()
    for(i in 1:k){x[[i]]  <- pt_world(model, u, do = ops[[i]])}
    query(x)}

  list(U = U, A = apply(U, 1, f), V = V)

  }


#' A function to plot worlds that satisfy a query
#'
#' @param which_worlds A worlds evaluation object generated by `pt_which`.
#' @param col1 Color parameter for query = TRUE
#' @param col2 Color parameter for query = FALSE
#'
#' @export
#'
#' @examples
#' qw <- pt_which(pt_model(),
#'             ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA)),
#'             query = function(x) (x[[1]][5] == 1) & (x[[2]][5] == 0))
#'  pt_worlds_plot(qw)

pt_worlds_plot <- function(which_worlds, col1 = "black", col2 = "pink") {
  pairs(which_worlds$U,
        col    = ifelse(which_worlds$A, col1, col2),
        pch = ".")}

#'
#' #'
#' #' Calculate the Probability that a binary clue K = 1, given W.
#' #' K enters as a vector of the form c(TRUE, FALSE, ...)
#' #'
#' #' @param model A model, made using pt_model
#' #' @param ops  A set of operations
#' #' @param query A query defined on outcomes on observables, V
#' #' @param sims optional number of simulations for draws; defaults to 500 if U not defined
#' #' @param W Preexisting data
#' #' @param K Clues being looked for
#' #' @export
#' #' @examples
#' #' pKw(model = pt_model(),
#' #'     ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA)),
#' #'     query = function(x) (x[[1]][5] == 1) & (x[[2]][5] == 0),
#' #'     sims = 3,
#' #'     W = matrix(c(0,1,1,0,0), 1, 5),
#' #'     K = matrix(c(1,0,0,0,0), 1, 5))
#' pKw <- function(model, ops, query, sims=NULL, U=NULL, W, K){
#'
#'   if(is.null(dim(K))) K <- matrix(K, 1)
#'   if(!is.matrix(K))   stop("Ks should be provide as a matrix with as many columns as elements in V")
#'   if(sum(K)>1)        stop("Include only one clue please")
#'   if(is.null(U))      U <- pt_worlds(model, N = sims, include_V = FALSE)
#'
#'     pt_worlds(model, U = U)
#'   V <- pt_worlds[model$var_names]
#' #  sims  <- nrow(U)
#'
#' #  result <- pt_which(model, U=U[paste0("U_", model$var_names)], ops = ops, query = query)
#'
#'   in_W <- sapply(1:sims, function(i) !(0 %in% (1*(result$V[i,] == W))))
#'
#'   # Two cases: when W is not possible and when it is:
#'   if(sum(in_W) ==0 ) {out <- NA
#'   } else {   out <- mean(result$V[in_W,K] == 1) } # ie
#'   out
#' }


#' Calculate what is learned from a set of observations, Ks, given existing data, W
#'
#' Defaults to the posterior variance on lambda_b - lambda_a for simple pt model
#' @param model A model, made using pt_model
#' @param ops  A set of operations
#' @param query A query defined on outcomes on observables, V
#' @param W  Observed variables: Vector indicating what is known on each node. Values NA, 0, 1 etc.
#' @param Ks TRUE/FALSE matrix with one column per node indicating which clues to be sought.
#' @param sims optional number of simulations for draws; defaults to 500 if U not defined
#' @param f any function; defaults to f=var, posterior variance. For posterior mean set f = mean
#' @param detail defaults to FALSE. When TRUE output includes listing of possible patterns and probability of each
#' @export
#' @examples
#' model <- pt_model()
#' ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA))
#' query = function(x) (x[[1]][5] ==1) & (x[[2]][5] == 0)
#' Ks = rbind(c(TRUE, FALSE, FALSE, FALSE, FALSE),
#'            c(TRUE, TRUE, FALSE, FALSE, FALSE),
#'            c(TRUE, TRUE, FALSE, TRUE, FALSE),
#'            c(TRUE, TRUE, TRUE, TRUE, TRUE))
#' pt_learning(model, ops, query,
#'             sims = 100,
#'             Ks = Ks)
#' pt_learning(model, ops, query,
#'             sims = 100,
#'             W = c(0, NA, NA, NA, NA),
#'             Ks = Ks)
#' U  <- as.matrix(t(replicate(100, model$P())))
#' pt_learning(model, ops, query,
#'             U = U,
#'             W = c(0, NA, NA, NA, NA),
#'             Ks = Ks)
#'
pt_learning <- function(
  model,
  ops,
  query,
  sims     = 2000,
  U        = NULL,
  W        = NULL,
  Ks       = NULL,     # Or a matrix with one column per node indicating whether to be sought or not
  f        = var,
  add_means = FALSE,
  detail    = FALSE,  # provide detailed output on clue probabilities
  detailround = 5     # Digits on rounding for detailed output
  ){
    if(is.null(W))             W <- rep(NA, length(model$var_names))
    if(!is.null(Ks) & is.null(dim(Ks))) Ks <- matrix(Ks, 1)
    if(!is.null(Ks) & !is.matrix(Ks)) stop("Ks should be provided as a matrix with as many columns as elements in V")
    if(is.null(U))             U  <- as.matrix(t(replicate(sims, model$P())))

    sims  <- nrow(U)
    V     <- pt_worlds(model, U=U)[model$var_names]
    A     <- pt_which(model, U=U,   ops = ops, query = query)$A

    prior <- f(A)  # First result

    in_W  <- sapply(1:sims, function(i) !(0 %in% (1*(V[i,] == W))))

    # Two cases:
    # Impossible W:
    if(all(!in_W)) { print("W is 0 probability event")

    } else {

    # Possible W
     posterior_given_W <- f(A[in_W])   # Second result

     # Now go through each of the K patterns given in Ks, check all implied permutations
     # of possible findings, get prob for each and post loss function (f, var) for each
     # Expected probability and variance
     # Expected probability and variance
     epv <- function(K, detail = FALSE) {

       Kranges <- apply(as.matrix(V[,K]), 2, function(j) length(unique(j)))

       if(is.null(Ks) | !(length(Kranges)>0)) return(posterior_given_W)

       # All permutations of the clues that are sought
       perms         <- perm(Kranges)
       patterns      <- matrix(NA, nrow(perms), ncol(V))
       patterns[, K] <- perms
       n_patterns    <- nrow(patterns)


       # Expected prob and variance over all the combinations of K
       prob_var <- sapply(1:n_patterns, function(j) {

         in_K <- sapply(1:sims, function(i) !(0 %in% (1*(V[i,] == patterns[j,]))))
         p_Kv <- mean(in_K[in_W])   # probability of K pattern, given W

         V_given_WK <- f(A[in_W & in_K]) # 0 var for 0 prob event
         if(sum(in_W & in_K) <= 1)   V_given_WK <- 0

         c(prob = p_Kv, variance = V_given_WK)
       })

       # If requestd this segment takes and returns the detailed weights and variance
       if(detail) {
         pattern <- sapply(1:n_patterns, function(i)  paste(patterns[i,], collapse = "-"))
         x <- rbind(pattern, round(prob_var,detailround))
         rownames(x) <- c("K-pattern", "prob pattern", "var|pattern");
         return(x)  }

       return(sum(prob_var[1,]*prob_var[2,]))

     }

    # Expected variance for each strategy
    e_p_v <- apply(Ks, 1, epv)
    }

    # Output
    ##################################################################################

    if(all(in_W == 0)){ posterior_given_W <- NA
                        e_p_v             <- rep(NA, nrow(Ks))
                       }

    if(is.null(Ks)) return(c(prior = prior, v_w = posterior_given_W))

    # Main output
    out <- c(prior, posterior_given_W, e_p_v)
    names(out) <- c("Prior v", "v |W", sapply(1:nrow(Ks), function(r) {paste(1*Ks[r,], collapse = "")} ))

    # Detail option
    if(detail)    out <- list(out, apply(Ks, 1, epv, detail = TRUE))

    out
  }


#' Get matrix of possible Ws: the set of possible data for a case when data is sought for *up to* m out of n nodes
#'
#' @param model a fitted pt model
#' @param m Number of W variations to make
#' @export
#' @examples
#' possible_W(3, 2)
#' possible_W(3, 1)
possible_W <- function(n, m) {
  W <- perm(rep(3, n))
  W[W==2] <- NA
  W <- W[apply(!is.na(W), 1, sum) <= m,]
  W}

#' Get matrix of possible Ks observation possibilities ( is observed or it is not observed)  on up to m out of n nodes
#'
#' @param model a fitted pt model
#' @export
#' @examples
#' possible_K(3,1)
#'
possible_K <- function(n, m) {
  K <- perm(rep(2, n))
  K <- K[rowSums(K) <= m,]
  K==1}

#' Asssess expected learning from a host of strategies.
#'
#' @param model a fitted pt model
#' @param m_w Number of Ws to look at
#' @param m_k Number of Ks to look at
#' @export
#' @examples
#' model <- pt_model()
#' ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA))
#' query = function(x) (x[[1]][5] ==1) & (x[[2]][5] == 0)
#' pt_strategies(model = model, ops = ops, query = query, sims = 100)

pt_strategies <- function(model,ops, query,
                            m_W = 1,
                            m_K = 1,
                            n   = length(model$var_names),
                            possible_Ws = possible_W(n, m_W),
                            possible_Ks = possible_K(n, m_K),
                            sims  = 2000,
                            round = 2,
                            Rsq   = FALSE,
                            f     = var,
                            U     = NULL,
                            symbol_unknown ="?"){

  if(!is.matrix(possible_Ws)) stop("possible_Ws should be a matrix with a column for each V")
  if(!is.matrix(possible_Ks)) stop("possible_Ks should be a matrix with a column for each V")
  if(is.null(U))              U  <- as.matrix(t(replicate(sims, model$P())))

  nW <- nrow(possible_Ws)
  nK <- nrow(possible_Ks)

  x <- sapply(1:nW, function(i) {
       pt_learning(model, ops, query, U = U, W  = possible_Ws[i,], Ks = possible_Ks, f=f)})

  if(Rsq) {
    x2 <- 1+t((t(x) - x[2,])/x[2,])
    x[3:nrow(x)] <- x2[3:nrow(x)]}

  x <- round(t(x), round)

  # Add rownames that show the W data pattern
  wn             <- possible_Ws
  wn[is.na(wn)]  <- symbol_unknown
  rownames(x) <-  sapply(1:nrow(wn), function(r) {paste(wn[r,], collapse = "")} )
  x
}





#' Declare a process tracing design from a causal model
#'
#' Warning -- do not use redesign with design produced from this designer
#'
#' @param true_model A model: not available to applied researcher
#' @param estimation_model A model: as supposed by applied researcher
#' @param ops Operations that form part of the query
#' @param query The query
#' @param strategy Data gathering strateg -- a vector of TRUES or FALSES indicating which nodes are observed
#' @param N Population for which inferences might be made. Estimand averages over these.
#' @param sims Precision for posterior calculation
#' @export
#' @examples
#' model <- pt_model()
#' ops = list(c(0, NA, NA, NA, NA), c(1, NA, NA, NA, NA))
#' query = function(x) (x[[1]][5] ==1) & (x[[2]][5] == 0)
#' design <- pt_designer(model, ops = ops, query = query, sims = 500)
#' draw_estimates(design)


pt_designer <- function(true_model,
                        estimation_model = true_model,
                        ops,
                        query,
                        strategy = TRUE,   # Nodes that get observed
                        N = 1,             # Population for which inferences might be made
                        sims = 2000        # Precision for posterior calculation
){

  believed_distn   <- pt_which(estimation_model,  ops = ops, query = query, sims = sims)

  estimand_function <- function(data, label) {
    U <- data[paste0("U_", true_model$var_names)]
    x <- mean(pt_which(true_model, ops = ops, query  = query, U = U)$A)
    data.frame(estimand_label = label, estimand = x, stringsAsFactors = FALSE)}

  # function to get the mean posterior estimate over a set of cases
  pt_est = function(data, believed_distn, strategy){

    # Actual data
    Ks  <- data[true_model$var_names]

    # Observed data
    Ks[!strategy] <- NA

    # Possible data, consistent with observed Ks
    possible <- apply(believed_distn$V, 1, function(x) !(0 %in%  ((Ks == x))))

    # Posteriors condition on possible
    data.frame(estimate      = mean(believed_distn$A[possible]),
               prior_var     = var(believed_distn$A),
               prior         = mean(believed_distn$A),
               posterior_var = var(believed_distn$A[possible]))}


  diagnosands <- declare_diagnosands(
    bias = mean(estimate - estimand, na.rm = TRUE),
    rmse = sqrt(mean((estimate - estimand)^2, na.rm = TRUE)),
    mean_estimand = mean(estimand, na.rm = TRUE),
    mean_prior_var = mean(prior_var, na.rm = TRUE),
    mean_posterior_var = mean(posterior_var, na.rm = TRUE),
    failed_estimation = mean(is.na(estimate)),
    keep_defaults = FALSE)

  design <-
    declare_population(data = pt_worlds(true_model, N = N)) +
    declare_estimand(handler = estimand_function, label = "Q") +
    declare_sampling(n = 1) +
    declare_estimator(handler = tidy_estimator(pt_est),
                      believed_distn = believed_distn,
                      strategy = strategy,
                      estimand = "Q")

  set_diagnosands(design, diagnosands)

}
