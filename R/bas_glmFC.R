# Copyright (c) 2024 Merlise Clyde and contributors to BAS. All rights reserved.
# This work is licensed under a GNU GENERAL PUBLIC LICENSE Version 3.0
# License text is available at https://www.gnu.org/licenses/gpl-3.0.html
#
#' Bayesian Adaptive Sampling Without Replacement for Variable Selection in
#' Generalized Linear Models
#'
#' Sample with or without replacement from a posterior distribution on GLMs
#'
#' BAS provides several search algorithms to find high probability models for
#' use in Bayesian Model Averaging or Bayesian model selection. For p less than
#' 20-25, BAS can enumerate all models depending on memory availability, for
#' larger p, BAS samples without replacement using random or deterministic
#' sampling. The Bayesian Adaptive Sampling algorithm of Clyde, Ghosh, Littman
#' (2010) samples models without replacement using the initial sampling
#' probabilities, and will optionally update the sampling probabilities every
#' "update" models using the estimated marginal inclusion probabilities. BAS
#' uses different methods to obtain the \code{initprobs}, which may impact the
#' results in high-dimensional problems. The deterministic sampler provides a
#' list of the top models in order of an approximation of independence using
#' the provided \code{initprobs}.  This may be effective after running the
#' other algorithms to identify high probability models and works well if the
#' correlations of variables are small to modest.  The priors on coefficients
#' are mixtures of g-priors that provide approximations to the power prior.
#'
#' @param formula generalized linear model formula for the full model with all
#' predictors, Y ~ X.  All code assumes that an intercept will be included in
#' each model.
#' @param family a description of the error distribution and link function for
#' exponential family; currently only `binomial()` with the logistic link and
#' `poisson()` and `Gamma()`with the log link are available.
#' @param data data frame
#' @param weights optional vector of weights to be used in the fitting process.
#' May be missing in which case weights are 1.
#' @param subset subset of data used in fitting
#' @param contrasts an optional list. See the contrasts.arg of `model.matrix.default()`.
#' @param offset a priori known component to be included in the linear
#' predictor; by default 0.
#' @param na.action a function which indicates what should happen when the data
#' contain NAs. The default is "na.omit".
#' @param n.models number of unique models to keep. If NULL, BAS will attempt
#' to enumerate unless p > 35 or method="MCMC". For any of methods using MCMC
#' algorithms that sample with replacement, sampling will stop when the number
#' of iterations exceeds 'MCMC.iterations'. On exit
#' 'n.models' is updated to reflect the unique number of models that have been
#' sampled.
#' @param betaprior Prior on coefficients for model coefficients (except
#' intercept).  Options include 
#' \code{\link{g.prior}}, 
#' \code{\link{CCH}},
#' \code{\link{robust}}, 
#' \code{\link{intrinsic}}, 
#' \code{\link{beta.prime}},
#' \code{\link{EB.local}}, 
#' \code{\link{AIC}}, and 
#' \code{\link{BIC}}.
#' @param modelprior Family of prior distribution on the models.  Choices
#' include \code{\link{uniform}}, \code{\link{Bernoulli}},
#' \code{\link{beta.binomial}}, truncated Beta-Binomial,
#' \code{\link{tr.beta.binomial}}, and truncated power family
#' \code{\link{tr.power.prior}}.
#' @param initprobs vector of length p with the initial inclusion probabilities
#' used for sampling without replacement (the intercept will be included with
#' probability one and does not need to be added here) or a character string
#' giving the method used to construct the sampling probabilities if "Uniform"
#' each predictor variable is equally likely to be sampled (equivalent to
#' random sampling without replacement). If "eplogp", use the
#' \code{\link{eplogprob}} function to approximate the Bayes factor using
#' p-values to find initial marginal inclusion probabilities and sample
#' without replacement using these inclusion probabilities, which may be
#' updated using estimates of the marginal inclusion probabilities. "eplogp"
#' assumes that MLEs from the full model exist; for problems where that is not
#' the case or 'p' is large, initial sampling probabilities may be obtained
#' using \code{\link{eplogprob.marg}} which fits a model to each predictor
#' separately.  To run a Markov Chain to provide initial
#' estimates of marginal inclusion probabilities, use method="MCMC+BAS" below.
#' While the initprobs are not used in sampling for method="MCMC", this
#' determines the order of the variables in the lookup table and affects memory
#' allocation in large problems where enumeration is not feasible.  For
#' variables that should always be included set the corresponding initprobs to
#' 1, to override the `modelprior` or use `include.always` to force these variables
#' to always be included in the model.
#' @param include.always A formula with terms that should always be included
#' in the model with probability one.  By default this is `~ 1` meaning that the
#' intercept is always included.  
#' This will also override any of the values in `initprobs`
#' above by setting them to 1.
#' @param method A character variable indicating which sampling method to use:
#' method="BAS" uses Bayesian Adaptive Sampling (without replacement) using the
#' sampling probabilities given in initprobs and updates using the marginal
#' inclusion probabilities to direct the search/sample; method="MCMC" combines
#' a random walk Metropolis Hastings (as in MC3 of Raftery et al 1997) with a
#' random swap of a variable included with a variable that is currently
#' excluded (see Clyde, Ghosh, and Littman (2010) for details);
#' method="MCMC+BAS" runs an initial MCMC as above to calculate marginal
#' inclusion probabilities and then samples without replacement as in BAS;
#' method = "deterministic" runs a deterministic sampler using the initial
#' probabilities (no updating); this is recommended for fast enumeration or if a
#' model of independence is a good approximation to the joint posterior
#' distribution of the model indicators.  For BAS, the sampling probabilities
#' can be updated as more models are sampled. (see 'update' below).  We
#' recommend "MCMC+BAS" or "MCMC" for high dimensional problems.
#' @param update number of iterations between potential updates of the sampling
#' probabilities in the "BAS" method. If NULL do not update, otherwise the
#' algorithm will update using the marginal inclusion probabilities as they
#' change while sampling takes place.  For large model spaces, updating is
#' recommended. If the model space will be enumerated, leave at the default.
#' @param bestmodel optional binary vector representing a model to initialize
#' the sampling. If NULL sampling starts with the null model
#' @param prob.rw For any of the MCMC methods, probability of using the
#' random-walk proposal; otherwise use a random "flip" move to propose a new
#' model.
#' @param burnin.iterations Number of iterations to discard as part of burnin
#' when using any of the MCMC
#' options; should be greater than 'n.models'. By default 10*p.
#' @param MCMC.iterations Number of MCMC iterations for sampling using any of the MCMC
#' options; should be greater than 'n.models'. By default 1000*p.
#' @param thin For "MCMC", thin the MCMC chain every "thin" iterations; default 
#' is no thinning.  For large p, thinning can be used to significantly reduce memory
#' requirements as models and associated summaries are saved only every thin 
#' iterations.  For thin = p, the  model and associated output are recorded 
#' every p iterations,similar to the Gibbs sampler in SSVS.
#' @param control a list of parameters that control convergence in the fitting
#' process.  See the documentation for \code{glm.control()}
#' @param laplace logical variable for whether to use a Laplace approximate for
#' integration with respect to g to obtain the marginal likelihood.  If FALSE
#' the Cephes library is used which may be inaccurate for large n or large
#' values of the Wald Chisquared statistic.
#' @param renormalize logical variable for whether posterior probabilities
#' should be based on renormalizing marginal likelihoods times prior
#' probabilities or use Monte Carlo frequencies. Applies only to MCMC sampling.
#' @param force.heredity Logical variable to force all levels of a factor to be
#' included together and to include higher order interactions only if lower
#' order terms are included.  Currently only supported with `method='MCMC'`
#' and `method='BAS'` (experimental).
#' Default is FALSE.
#' @param GROW Logical variable to indicate that the output vectors in MCMC are growable.  
#' Rather than allocate space based on `n.models`, the vectors will grow as needed
#' if the number of unique models sampled exceeds the initial size of the allocated output vectors controlled by 
#' `n.models.init`.  This is useful when `n.models` is unknown
#' before reaching 'MCMC.iterations'.  Default is TRUE.
#' @param expand variable to control how much to grow vectors with GROW = TRUE 
#' if number of unique models exceeds the current size of the vectors. 
#' The default is 1.25 times, which allows vectors to grow by 25 percent.
#' @param n.models.init Initial size of output vectors if GROW = TRUE. The 
#' default is `n.models = 2500`.
#' @param bigmem Logical variable to indicate that there is access to
#' large amounts of memory (physical or virtual) for enumeration
#' with large model spaces, e.g. > 2^25.
#'
#' @return \code{bas.glm} returns an object of class \code{basglm}
#'
#' An object of class \code{basglm} is a list containing at least the following
#' components:
#'
#' \item{postprobs}{the posterior probabilities of the models selected}
#' \item{priorprobs}{the prior probabilities of the models selected}
#' \item{logmarg}{values of the log of the marginal likelihood for the models}
#' \item{n.vars}{total number of independent variables in the full model,
#' including the intercept} 
#' \item{size}{the number of independent variables in
#' each of the models, includes the intercept} 
#' \item{which}{a list of lists
#' with one list per model with variables that are included in the model}
#' \item{probne0}{the posterior probability that each variable is non-zero}
#' \item{mle}{list of lists with one list per model giving the GLM
#' estimate of each (nonzero) coefficient for each model.} 
#' \item{mle.se}{list of
#' lists with one list per model giving the GLM standard error of each
#' coefficient for each model} 
#' \item{deviance}{the GLM deviance for each model}
#' \item{modelprior}{the prior distribution on models that created the BMA 
#' object} 
#' \item{Q}{the Q statistic for each model used in the marginal
#' likelihood approximation} 
#' \item{Y}{response} 
#' \item{X}{matrix of predictors}
#' \item{family}{family object from the original call} 
#' \item{betaprior}{family object for prior on coefficients, including
#'  hyperparameters}
#' \item{modelprior}{family object for prior on the models}
#' \item{include.always}{indices of variables that are forced into the model}
#' @author Merlise Clyde (\email{clyde@@duke.edu}), Quanli Wang and Yingbo
#' Li
#' @references Li, Y. and Clyde, M. (2018) Mixtures of g-priors in Generalized
#' Linear Models. 
#' Journal of the American Statistical Association. 113:1828-1845 \cr
#' \doi{10.1080/01621459.2018.1469992} \cr
#' Clyde, M. Ghosh, J. and Littman, M. (2010) Bayesian Adaptive Sampling for
#' Variable Selection and Model Averaging. Journal of Computational Graphics
#' and Statistics.  20:80-101 \cr
#' \doi{10.1198/jcgs.2010.09049} \cr
#' Raftery, A.E, Madigan, D. and Hoeting, J.A. (1997) Bayesian Model Averaging
#' for Linear Regression Models. Journal of the American Statistical
#' Association.
#' @keywords GLM regression
#' @examples
#'
#' library(MASS)
#' data(Pima.tr)
#'
#'
#' # enumeration  with default method="BAS"
#' pima.cch = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="BAS",
#'               betaprior=CCH(a=1, b=532/2, s=0), family=binomial(),
#'               modelprior=beta.binomial(1,1))
#'
#' summary(pima.cch)
#' image(pima.cch)
#'
#' # Note MCMC.iterations are set to 2500 for illustration purposes due to time
#' # limitations for running examples on CRAN servers.
#' # Please check convergence diagnostics and run longer in practice
#'
#' pima.robust = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="MCMC", MCMC.iterations=2500,
#'               betaprior=robust(), family=binomial(),
#'               modelprior=beta.binomial(1,1))
#'
#' pima.BIC = bas.glm(type ~ ., data=Pima.tr, n.models= 2^7,
#'               method="MCMC+BAS", MCMC.iterations=2500,
#'               betaprior=bic.prior(), family=binomial(),
#'               modelprior=uniform())

#' # Poisson example
#' if(requireNamespace("glmbb", quietly=TRUE)) {
#'   data(crabs, package='glmbb')
#'   #short run for illustration
#'   crabs.bas = bas.glm(satell ~ color*spine*width + weight, data=crabs,
#'                       family=poisson(),
#'                       betaprior=EB.local(), modelprior=uniform(),
#'                       method="MCMC", n.models=2^10, MCMC.iterations=2500,
#'                       prob.rw=.95)
#'   
#'  # Gamma example
#'  if(requireNamespace("faraway", quietly=TRUE)) {
#'     data(wafer, package='faraway')
#'                       
#'     wafer_bas = bas.glm(resist~ ., data=wafer,  include.always = ~ .,
#'                         betaprior = bic.prior() ,
#'                         family = Gamma(link = "log"))
#'   }
#' }
#' @concept BMA
#' @concept variable selection
#' @family BMA functions
#' @rdname bas.glm
#' @export
bas.glmFC <- function(formula, family = binomial(link = "logit"),
                    data, weights, subset, contrasts=NULL, offset, na.action = "na.omit",
                    n.models = NULL,
                    betaprior = CCH(alpha = .5, beta = as.numeric(nrow(data)), s = 0),
                    modelprior = beta.binomial(1, 1),
                    initprobs = "Uniform",
                    include.always = ~1,
                    # NEW: cost vector for variables (if NULL, no cost is considered / all costs equal to 1, the baseline)
                    # b is introduced in the modelprior argument
                    var.costs = NULL,  
                    method = "MCMC",
                    update = NULL,
                    bestmodel = NULL,
                    prob.rw = 0.5,
                    burnin.iterations = NULL, MCMC.iterations = NULL, thin = 1,
                    control = glm.control(), laplace = FALSE, renormalize = FALSE,
                    force.heredity = FALSE, GROW = TRUE, expand = 1.25, n.models.init = 2500,
                    bigmem = FALSE) {
  num.updates <- 10
  call <- match.call()

  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  
  if (is.function(family)) {
    family <- family()
  }

  if (!(family$family %in% c("binomial", "poisson", "Gamma"))) {
    stop(paste("family ", family$family, "not implemented"))
  }
  else {
    if (family$family == "binomial") {
      if (family$link != "logit") {
        stop("Only logit link is implemented for binomial family currently")
      }
    }
    if (family$family == "Gamma") {
      if (family$link != "log") {
        stop("Only log link is implemented for Gamma family currently")
      }
    }
    if (family$family == "poisson") {
      if (family$link != "log") {
        stop("Only log link is implemented for poisson family currently")
      }
  }
  }
  
  
  if (missing(data)) {
    data <- environment(formula)
  }
  
  if (!inherits(modelprior, "prior")) stop("modelprior should be an object of class prior,  uniform(),  beta.binomial(), etc")

  if (!(method %in% c("BAS", "deterministic", "MCMC", "MCMC+BAS", "AMCMC"))) {
    stop(paste("No available sampling method:", method))
  }
  

  mfall <- match.call(expand.dots = FALSE)
  m <- match(c(
    "formula", "data", "subset", "weights", "na.action",
    "etastart", "mustart", "offset"
  ), names(mfall), 0L)
  mf <- mfall[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  n.NA <- length(attr(mf, "na.action"))

  if (n.NA > 0) {
    warning(paste(
      "dropping ", as.character(n.NA),
      "rows due to missing data"
    ))
  }

  Y <- model.response(mf, type = "any")

  mt <- attr(mf, "terms")
  X <- model.matrix(mt, mf, contrasts)
  #    Y = glm.obj$y
  #    X = glm.obj$x
  # Building the overparameterized model matrix
  depvars <- attr (mt, "term.labels") # Names of the independent (competing) variables 
  # This is valid if no variable is forced to be included in the model (i.e., include.always = ~1)
  # Otherwise, we could fit the new "null model" as in BayesVarSel
  # Only assign contrasts to factor variables
  factor_vars <- names(Filter(is.factor, mf))
  contrast.list <- list()
  if (length(factor_vars) > 0) {
    contrast.list <- setNames(
      lapply(mf[factor_vars], function(col) {
        if (is.factor(col)) {
          contrasts(col, contrasts = FALSE)
        } else {
          NULL
        }
      }),
      factor_vars
    )
  }
  else { # No factor in the set of competing variables...
    contrast.list <- NULL
  }
  X <- model.matrix(mt, mf, contrasts.arg = contrast.list) # Create the design matrix X with all levels of factors
  
  namesx <- dimnames(X)[[2]]
  namesx[1] <- "Intercept"

  p <- dim(X)[2] # Now (overparameterized design matrix): p = k + sum_j Lj + 1, as it includes the intercept
  
  # positions is a matrix with number of rows equal to the number of regressors
	# (either factor or numeric) and number of columns the number of columns of X
	# Each row describes the position (0-1) in X of a regressor (several positions in case
	# this regressor is a factor)

  positions <- matrix (0L,
    ncol = p - 1, # Number of columns in X (after removing the intercept)
    nrow = length(depvars) # Number of competing variables (length(depvars))
  ) 

  cols <- colnames(X)[-1] # drop intercept

  # For each competing variable (row)
  for (i in seq_along(depvars)) {
    pat <- paste0("^", depvars[i], "$")
    ## pat <- paste0("^", depvars[i], "(\\d+|$)")  # digit(s) OR end-of-string
    positions[i, ] <- as.integer(grepl( # searchs for patterns in character strings and returns a logical vector (indicating whether the pattern was found)
      pat, # text or regular expression to search for
      cols, # character vector in which to search
      perl = TRUE))
  }

  # Vector of length equal to the number of competing variables, with the number of levels of each variable (1 for numeric, >1 for factors)
  var_levels <- colSums (positions %*% t(positions))
  # rownames (var_levels) <- depvars

  # positionsX is a vector of the same length as rows has X (overparameterized design matrix) = number of competing variables
  # with 1 in the position with a numeric variable:
  positionsx <- as.numeric (var_levels == 1) 
  rownames (positions) <- depvars
  
   
  if (is.null(var.costs)) {
   # A questão dos discounted predictors nunca vai ser ativada, porque estamos a colocar 0s na primeira coluna
    var.costs <- matrix (1.00, nrow = length(depvars), ncol = length(depvars) + 1)  # initialize with 0s
    var.costs [, 1] <- 0.00
    # This way the baseline cost will be 1
  } else {

    # Using drop, because I don´t want this matrix to reduce to a vector when we only have 1 competing variable
    var.costs <- var.costs [depvars, c("Groups", depvars), drop = FALSE]
    # Remove row and column names
    rownames(var.costs) <- NULL
    colnames(var.costs) <- NULL

  }

  diff_groups <- unique(var.costs[, 1]) 

  groups <- matrix (0L,
      nrow = length (diff_groups), # Number of groups
      ncol = nrow (var.costs), # Number of competing variables
      dimnames = list(as.character(diff_groups), depvars) # Row and Column Names, respectively
  ) 
  
  for (i in seq_along(diff_groups)) {
    groups[i, ] <- as.integer(var.costs[, 1] == diff_groups[i])
  }

  if (nrow (var.costs) != length (depvars)) {
    stop ("var.costs must be a matrix with the same rows as the number of competing variables")
  } 
  

  nobs <- dim(X)[1]
  if (nobs == 0) {stop("Sample size is zero; check data and subset arguments")}
  #   weights = as.vector(model.weights(mf))

  weights <- as.vector(model.weights(mf))
  if (is.null(weights)) {
    weights <- rep(1, nobs)
  }

  offset <- model.offset(mf)
  if (is.null(offset)) offset <- rep(0, nobs)

  null.model = glm(Y ~ 1,
                      offset = offset,
                      family = eval(call$family))

  null.deviance = null.model$null.deviance
  loglik_null <- as.numeric(-0.5 * null.deviance)
 



  if (!is.numeric(initprobs)) {
    if (nobs <= p && initprobs == "eplogp") {
      stop(
        "Full model is not full rank so cannot use the eplogp bound to create starting sampling probabilities, perhpas use 'marg-eplogp' for fiting marginal models\n"
      )
    }
    initprobs <- switch(
      initprobs,
      "eplogp" = eplogprob(glm(Y ~ X - 1,
                               family = family, weights = weights,
                               offset = offset)),
      "marg-eplogp" = eplogprob.marg(Y, X),
      "uniform" = c(1.0, rep(.5, p - 1)),
      "Uniform" = c(1.0, rep(.5, p - 1))
    )
  }
  if (length(initprobs) == (p - 1)) {
    initprobs <- c(1.0, initprobs)
  }

  # set up variables to always include
  keep <- 1
  if ("include.always" %in% names(mfall)) {
    minc <- match(c("include.always", "data", "subset"), names(mfall), 0L)
    mfinc <- mfall[c(1L, minc)]
    mfinc$drop.unused.levels <- TRUE
    names(mfinc)[2] <- "formula"
    mfinc[[1L]] <- quote(stats::model.frame)
    mfinc <- eval(mfinc, parent.frame())
    mtinc <- attr(mfinc, "terms")
    X.always <- model.matrix(mtinc, mfinc, contrasts)

    keep <- c(1L, match(colnames(X.always)[-1], colnames(X)))
    initprobs[keep] <- 1.0
    if (ncol(X.always) == ncol(X)) {
      # just one model with all variables forced in
      # use method='BAS" as deterministic and MCMC fail in this context
      method <- "BAS"
    }
  }

  parents <- matrix(1, 1, 1)
  force.heredity <- FALSE # As in García-Donato and Paulo (2022)

  # if (method == "deterministic" | method == "MCMC+BAS" )
  #  force.heredity <- FALSE
  # 
  # if (force.heredity) {
  #   parents <- make.parents.of.interactions(mf, data)
  #
  #   # check to see if really necessary
  #   if (sum(parents) == nrow(parents)) {
  #     parents <- matrix(1, 1, 1)
  #     force.heredity <- FALSE
  #   }
  # }

  prob <- normalize.initprobs.lm(initprobs, p)

  # bestmodel: model to initialize the sampling
    if (is.null(bestmodel)) { 
    #    bestmodel = as.integer(initprobs)
    bestmodel <- c(1, rep(0, p - 1))
  }
  bestmodel[keep] <- 1
  
  # if (force.heredity) {
  #   update <- NULL # do not update tree  FIXME LATER
  #   if (prob.heredity(bestmodel, parents) == 0) {
  #     warning("bestmodel violates heredity conditions; resetting to null model.  Please check  include.always and bestmodel")
  #     bestmodel <- c(1, rep(0, p - 1))
  #   }
  # # initprobs <- c(1, seq(.95, .55, length = (p - 1))) # keep same order
  # }

  bestmodel <- as.integer(bestmodel)


  if (!GROW & method == "MCMC") method <- "MCMC_OLD"
  if (!GROW & method == "BAS") method <- "BAS_OLD"
  
  # p = k + sum_j Lj + 1, as it includes the intercept
  # Shouldn´t we adjust n.models to accomodate the new model space 
  # (as some variables may be forced to be included in the model: e.g., the intercept)?
  if (is.null(n.models)) {
    # n.models <- min(2^p, 2^16) # # change to 2^p to force enumeration regardless p
    n.models <- min(2^(p-1), 2^16) # I´ll always include the intercept / change to 2^(p-1) to force enumeration regardless p
    if (method == "MCMC")  n.models = min(n.models, n.models.init) 
    # FIXME add n.models.init as argument rather than specify here
  }
  if (is.null(MCMC.iterations)) {
    #MCMC.iterations <- as.integer(p * 1000)
    MCMC.iterations <- as.integer(max(10000, p * 1000)) # David
  }
  if (is.null(burnin.iterations)){
    burnin.iterations <- as.integer(p * 25) # We always have a burn-in period
  }
  
  
  n.models <- as.integer(normalize.n.models(n.models, p, prob, method, bigmem))

  modelprior <- normalize.modelprior(modelprior, p)

  modeldim <- as.integer(rep(0, n.models))

  #print(MCMC.iterations)

  # For method = "BAS"
  if (is.null(update)) {
    # if (force.heredity) {  # do not update tree for BAS
    #  update <- n.models + 1}
    # else {

      if (n.models == 2^(p - 1)) {
        update <- n.models + 1
      } else {
        (update <- n.models / num.updates)
      }

    # }
  }

  #  check on priors
  
  if (!inherits(betaprior, "prior")) stop("prior on coeeficients must be an object of type 'prior'")
  

  betaprior$hyper.parameters$loglik_null <- loglik_null
  #  	browser()

  if (betaprior$family == "BIC" & is.null(betaprior$n)) {
    betaprior <- bic.prior(as.numeric(nobs))
  }


  if (betaprior$family == "hyper-g/n" & is.null(betaprior$n)) {
    betaprior$hyper.parameters$theta <- 1 / nobs
    betaprior$n <- nobs
  }
  if (betaprior$family == "robust" & is.null(betaprior$n)) betaprior <- robust(as.numeric(nobs))

  if (betaprior$family == "intrinsic" & is.null(betaprior$n)) {
      betaprior$hyper.parameters$n <- as.numeric(nobs)
  }

  if (betaprior$family == "betaprime" & is.null(betaprior$hyper.parameters$n)) {
    betaprior$hyper.parameters$n <- as.numeric(nobs)
  }

  # call this to coerce response as needed for glm
 
  y = Y
  eval(family$initialize)
  storage.mode(y) <- "double"
  
  # an R function that directly invokes a compiled C function.
  # Added methods = c("Gibbs","GibbsBVS")
  result <- switch(method,
    "MCMC_OLD" = .Call(C_glm_mcmc,
      RY = y, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw), # Probability of using the the random SWAP proposal
      BURNIN_Iterations = as.integer(burnin.iterations),
      MCMC_Iteration = as.integer(MCMC.iterations),
      Rthin = as.integer(thin),
      family = family, Rcontrol = control,
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "MCMC" = .Call(C_glm_mcmc_grow,
                   RY = y, X = X,
                   Roffset = as.numeric(offset),
                   Rweights = as.numeric(weights),
                   Rprobinit = prob,
                   RnModels = as.integer(n.models),
                   modelprior = modelprior,
                   betaprior = betaprior,
                   positions = positions,
                   levels = var_levels,
                   costs = var.costs,
                   Rbestmodel = bestmodel,
                   plocal = as.numeric(1.0 - prob.rw),
                   BURNIN_Iterations = as.integer(burnin.iterations),
                   MCMC_Iteration = as.integer(MCMC.iterations),
                   Rthin = as.integer(thin),
                   family = family, Rcontrol = control,
                   Rlaplace = as.integer(laplace),
                   Rparents = parents, Rexpand = as.numeric(expand)
    ),
    "Gibbs" = .Call(C_glm_gibbssampler,
      RY = y, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw), # Probability of using the the random SWAP proposal
      BURNIN_Iterations = as.integer(burnin.iterations),
      MCMC_Iteration = as.integer(MCMC.iterations),
      Rthin = as.integer(thin),
      family = family, Rcontrol = control,
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "GibbsBVS" = .Call(C_glm_gibbsBVS,
      RY = y, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw), # Probability of using the the random SWAP proposal
      BURNIN_Iterations = as.integer(burnin.iterations),
      MCMC_Iteration = as.integer(MCMC.iterations),
      Rthin = as.integer(thin),
      family = family, Rcontrol = control,
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "BAS" = .Call(C_glm_sampleworep_grow,
      RY = y, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      RnModels = as.integer(n.models),
      #Rmodeldim = modeldim
      modelprior = modelprior,
      betaprior = betaprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw),
      family = family, Rcontrol = control,
      Rupdate = as.integer(update),
      Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "BAS_OLD" = .Call(C_glm_sampleworep,
                  RY = y, X = X,
                  Roffset = as.numeric(offset),
                  Rweights = as.numeric(weights),
                  Rprobinit = prob,
                  RnModels = as.integer(n.models),
                  #Rmodeldim = modeldim
                  modelprior = modelprior,
                  betaprior = betaprior,
                  betaprior = betaprior,
                  positions = positions,
                  levels = var_levels,
                  costs = var.costs,
                  Rbestmodel = bestmodel,
                  plocal = as.numeric(1.0 - prob.rw),
                  family = family, Rcontrol = control,
                  Rupdate = as.integer(update),
                  Rlaplace = as.integer(laplace),
                  Rparents = parents
    ),
    "MCMC+BAS" = .Call(C_glm_mcmcbas,
      RY = y,
      X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      Rbestmodel = bestmodel,
      plocal = as.numeric(1.0 - prob.rw),
      BURNIN_Iterations = as.integer(burnin.iterations),
      MCMC_Iteration = as.integer(MCMC.iterations),
      family = family, Rcontrol = control,
      Rupdate = as.integer(update), Rlaplace = as.integer(laplace),
      Rparents = parents
    ),
    "deterministic" = .Call(C_glm_deterministic,
      RY = y, X = X,
      Roffset = as.numeric(offset),
      Rweights = as.numeric(weights),
      Rprobinit = prob,
      Rmodeldim = modeldim,
      modelprior = modelprior,
      betaprior = betaprior,
      positions = positions,
      levels = var_levels,
      costs = var.costs,
      family = family,
      Rcontrol = control,
      Rlaplace = as.integer(laplace)
    )
  )


  result$namesx <- namesx
  result$n <- nrow(X)
  result$modelprior <- modelprior
  result$probne0[keep]  <- 1.0
  result$probne0.RN <- result$probne0
  result$postprobs.RN <- result$postprobs
  result$family <- family
  result$betaprior <- betaprior
  result$modelprior <- modelprior

  result$n.models <- length(result$postprobs)
  result$include.always <- keep

  #-------
  # David
  #------
  result$costs <- var.costs
  result$depvars <- depvars
  result$positions <- positions
  result$positionsx <- positionsx
  result$groups <- groups

  df <- rep(nobs - 1, result$n.models)

  if (betaprior$class == "IC") df <- df - result$size + 1
  
  result$df <- df
  result$R2 <- 1.0 - result$deviance/null.deviance
  result$n.vars <- p
  result$Y <- y
  result$X <- X
  result$weights = weights
  result$call <- call
  result$terms <- mt
  result$contrasts <- attr(X, "contrasts")
  result$xlevels <- .getXlevels(mt, mf)
  result$model <- mf

  # drop null model if it is present
  if (betaprior$family == "Jeffreys" & (min(result$size) == 1)) result <- .drop.null.bas(result)
  
  # github issue #74. drop models with zero prior probability
  
  if (any(result$priorprobs == 0)) {
    drop.models = result$priorprobs != 0
    
    result$mle = result$mle[drop.models]
    result$mle.se = result$mle.se[drop.models]
    result$mse = result$mse[drop.models]
    result$which = result$which[drop.models]
    result$freq = result$freq[drop.models]
    result$shrinkage = result$shrinkage[drop.models]
    result$R2 = result$R2[drop.models]
    result$logmarg = result$logmarg[drop.models]
    result$df = result$df[drop.models]
    result$size = result$size[drop.models]
    result$Q = result$Q[drop.models]
    result$rank = result$rank[drop.models]
    result$sampleprobs = result$sampleprobs[drop.models]
    result$postprobs = result$postprobs[subset = drop.models]
    result$priorprobs = result$priorprobs[subset = drop.models]
    result$intercept = result$intercept[drop.models]
    result$deviance = result$deviance[drop.models]
    result$n.models = length(result$postprobs)
  }

  list_models <- result$which
  
  models_matrix <- t(
    as.matrix(
      sapply(list_models,
      function(x, dim_p) {
        xx <- rep(0, dim_p)
        xx[x + 1] <- 1
        xx
      },
      p, 
      simplify = "matrix"
      )
    )
  )
  colnames(models_matrix) <- namesx # Includes the intercept...

  if (method == "MCMC" | method == "MCMC_OLD" | method == "Gibbs" | method == "GibbsBVS") {
    
    # Resampling the models found in MCMC to correct for the real set of competing models...
    if (modelprior$family == "SBSB") {new_priorprobs <- priorSBSB2 (models_matrix [, -1, drop = FALSE], positions)}
    if (modelprior$family == "CC") {new_priorprobs <- priorConstConst2 (models_matrix [, -1, drop = FALSE], positions)}
    if (modelprior$family == "SBC") {new_priorprobs <- priorSBConst2 (models_matrix [, -1, drop = FALSE], positions)}
    if (modelprior$family == "SB") {new_priorprobs <- priorSB2 (models_matrix [, -1, drop = FALSE], positions)}
    if (modelprior$family == "Uniform") {new_priorprobs <- priorConst2 (models_matrix [, -1, drop = FALSE], positions)}
    
    # Ver se faz sentido usar FND como está...
    if (modelprior$family == "FND") {
      new_priorprobs <- priorFND2 (models_matrix [, -1, drop = FALSE],
                                   positions, var.costs [, -1], modelprior$hyper.parameters, nobs)
    }
    if (modelprior$family == "FNDConst") { 
      new_priorprobs <- priorFNDConst2 (models_matrix [, -1, drop = FALSE],
                                        positions, var.costs, modelprior$hyper.parameters, nobs)
    }
    if (modelprior$family == "FNDSB") { 
      new_priorprobs <- priorFNDSB2 (models_matrix [, -1, drop = FALSE],
                                     positions, var.costs, modelprior$hyper.parameters, nobs)
    }

    # Seria engraçado perceber quantos modelos é que já são saturados antes do resampling; para ver se o resampling é bom ou não    
    result$num_bad_models <- sum (new_priorprobs == 0)
    
    # Create bvs_df first
    bvs_df <- data.frame (
      models_matrix,
      freq = result$freq,
      logmarg = result$logmarg,
      new_priorprobs = new_priorprobs,
      priorprobs = result$priorprobs, 
      model_num = seq_len (nrow(models_matrix)) # To keep track of the models
    )
    
    # Expanding the dataframe for resampling 
    idx_expanded <- rep.int(seq_len(nrow(bvs_df)), times = bvs_df$freq)
    result$MCMC.iter <- sum (result$freq) # Should be equal to length (idx_expanded) 
    
    expanded_df <- bvs_df [idx_expanded, ] # Care needed when dealing with the freq column
    
    # A questão agora é mesmo: modelos que foram visitados mais vezes, deviam ser mais prováveis de serem resampled:
    # se não dissermos que eles são iguais, garantimos isso fazendo com que apareçam mais vezes como candidatos,
    # uma vez que expandi o dataframe?
    
    # Would it be equivalent multiplying the prob argument in the sample function by freq?
    # We would then change X for seq_len (nrow(bvs_df))? Size would still need to be the same
    # , so we could still define posterior model probabilities according to Monte Carlo
    
    resamp <- sample ( # vector with the indices of the models to be resampled
     x = seq_len (nrow(expanded_df)),
     # Number of models to resample (number of MCMC iterations)
     size = nrow(expanded_df), # Same as sum (bvs_df$freq), I think 
     # Allows duplicate models to be selected
     replace = TRUE,
     # Importance resampling of model indices
     prob = expanded_df$new_priorprobs / expanded_df$priorprobs 
     # Odds = Prior Probs Corrected / Prior Probs Wrongly Defined (Original)
    ) 
    
    resampled_df <- expanded_df [resamp, ]
    
    # Count occurrences of each index
    # (similar to an histogram of the resampling over all rows: indexes as the x-axis and counts as the y-axis)
    draw_counts <- tabulate (resamp, nbins = nrow(expanded_df))
    
    # resampled_df$model_num maps every expanded row (after resampling) back to the unique model
    unique_model_num <- unique (resampled_df$model_num)
    
    # new_freq <- as.numeric (rowsum (draw_counts, group = resampled_df$model_num, reorder = FALSE)) (before)
    # collapse the counts by the original model number
    freq_by_model <- rowsum (draw_counts, # rowsum() returns a 1-col matrix
                             group   = expanded_df$model_num,
                             reorder = FALSE)
    
    # Model numbers should be integer (not character)
    model_ids <- as.integer(rownames(freq_by_model))
    
    # keep only the models that actually appeared in the resample
    nonzero <- freq_by_model [, 1] > 0
    bvs_df_resamp <- bvs_df [model_ids [nonzero], ]
    
    # The old freq column doesn´t seem very useful, to be honest...
    result$freq <- bvs_df_resamp$freq # Number of times each model kept in the resampling process was visited by MCMC
    
    # Updating the freq column (safe indexing — aligned by construction)
    bvs_df_resamp$freq <- freq_by_model [nonzero, 1] # Number of times each unique model was drawn in the resampling process
    # Update result object with resampled values
    result$resampling_freq <- bvs_df_resamp$freq
     
    result$logmarg <- bvs_df_resamp$logmarg
    result$n.models <- nrow (bvs_df_resamp) # New number of unique models
    
    result$priorprobs <- bvs_df_resamp$new_priorprobs
    result$old_priorprobs <- bvs_df_resamp$priorprobs # Old prior probs used in the MCMC algorithm
     
    models_matrix <- models_matrix [unique_model_num, , drop = FALSE]
    

    # The denominator in the expression below is the number of MCMC iterations 
    result$postprobs.MCMC <- result$resampling_freq / sum(result$resampling_freq) # freq is the number of times each model was visited by MCMC
    # PIPs for the variables in each model (at the level of the levels)
    result$probne0.MCMC <- colMeans (resampled_df [, c(1:p), drop = FALSE]) # Must include the intercept
    

    result$old_ratio <- result$n.Unique / (2 ^ (p - 1)) # The denominator includes repeated models
    # We´re treating numerical variables as if they were categorical variables with two levels...
    levels <- pmax (rowSums (positions), 2)
    result$new_ratio <- nrow (bvs_df_resamp) / prod ((2 ^ levels) - levels) # The denominator portraits the real number of competing models
    
    result$postprobs.RN <- compute_posterior (result$logmarg, result$priorprobs)
    
    result$probne0.RN <- as.vector (result$postprobs.RN %*% models_matrix [, -1])     
    
    # freq: number of times each model was visited by MCMC
    # sum(result$freq) is just the number of MCMC iterations
    result$postprobs.MCMC <- result$freq / sum(result$freq) 
    
    # If true, we´re basing posterior probabilities on renormalizing marginal likelihoods times prior probabilities
    if (!renormalize) { # If false, we base posterior probabilities on relative frequencies from MCMC
      result$probne0 <- result$probne0.MCMC # empirical posterior inclusion probabilities
      result$postprobs <- result$postprobs.MCMC #  empirical posterior model probabilities, computed as relative frequencies.
    }
  }

  log_marg <- result$logmarg  # log(marginal likelihood) of each model (including repeated ones)
  models_w_logmarg <- cbind (models_matrix, log_marg)
  # Matrix with the binary representation of the models and the log marginal likelihood of each model
  colnames(models_w_logmarg) <- c(namesx, "log_marg") 
  
  # Trying to identify the Null Model index (sometimes it is not the first model to be enumerated)
  only_zero <- sapply(list_models, function(x) length(x) == 1 && x == 0)
  null_model_index <- which (only_zero)

  if (length(null_model_index) > 0L) {
    
    lBFi0 <- exp(log_marg - log_marg[null_model_index[1L]]) # log(BF) of each model (including repeated ones) to the null model
    models_w_lBF <- cbind(models_matrix, log_BFi0 = lBFi0)

    #colnames(models_w_lBF) <- c (namesx, "log_BFi0")

    result$lBFi0 <- lBFi0
    result$models_w_lBF <- models_w_lBF [, -1]
  }

  # Now we convert the model indicator matrix to specify only the active variables in each model (not looking at the levels of factors)
	if (dim(positions)[1] > 1) # if there is more than one competing variable
    models_active_vars <- t( # I need to transpose the matrix to have the variables back as columns and rows as models
      apply ( # Produces a matrix where each column corresponds to the result of applying the function to one row of models_matrix (results are stacked column-wise).
      models_matrix [, -1, drop = FALSE], # positions does not include the intercept, so we remove the first column of models_matrix
      MARGIN = 1, # Apply function to each row (model) of the matrix
      # The matrices product computes the number of levels (or variables) active per group or factor
      FUN = function(x,M){as.numeric(x %*% M>0)}, # Tells us what are the active variables in each model (row) (not at the level of the levels in the factors)
      M = t(positions)
      ))
	else models_active_vars <- matrix (
      as.numeric(models_matrix [, -1, drop = FALSE] %*% t(positions)>0),
      ncol = 1 #, nrow = dim(models_matrix)[1], i.e., number of models
      )
    
  colnames(models_active_vars) <- depvars
  # matrix product symbol, %*%, is sometimes confused with the pipe operator, %>%
  pip <- as.vector (result$postprobs %*% models_active_vars)
  names (pip) <- depvars # Posterior Inclusion Probabilities at the level of the variables 
    
	if (dim(groups)[1] > 1) # if there is more than one group
    models_active_groups <- t ( # I need to transpose the matrix to have the groups back as columns and rows as models
      apply ( # Produces a matrix where each column corresponds to the result of applying the function to one row of models_matrix (results are stacked column-wise).
      models_active_vars [, , drop = FALSE],
      MARGIN = 1, # Apply function to each row (model) of the matrix
      # The matrices product computes the number of variables active per group
      FUN = function(x,M){as.numeric(x %*% M>0)}, # Tells us the active groups in each model (row)
      M = t(groups)
      )
    )
	else models_active_groups <- matrix (
        as.numeric(models_active_vars [, , drop = FALSE] %*% t(groups) > 0),
        ncol = 1
       )

  colnames(models_active_groups) <- as.character(diff_groups)
  pip_groups <- as.vector (result$postprobs %*% models_active_groups)
  # Posterior Inclusion Probabilities at the level of the variables
  names (pip_groups) <- as.character (diff_groups)
  result$pip_groups  <- pip_groups

  if (method == "MCMC" | method == "MCMC_OLD" | method == "Gibbs" | method == "GibbsBVS") {

    var_pip.RN <- as.vector (result$postprobs.RN %*% models_active_vars)
    names (var_pip.RN) <- depvars
    result$var_pip.RN <- var_pip.RN
    
    group_pip.RN <- as.vector (result$postprobs.RN %*% models_active_groups)
    names (group_pip.RN) <- as.character (diff_groups)
    result$group_pip.RN <- group_pip.RN

    # Weights
    w <- result$resampling_freq

    var_pip.MCMC <- colSums (models_active_vars * w) / sum (w)
    names (var_pip.MCMC) <- depvars
    result$var_pip.MCMC <- var_pip.MCMC
    
    group_pip.MCMC <- colSums (models_active_groups * w) / sum (w)
    names (group_pip.MCMC) <- as.character (diff_groups)
    result$group_pip.MCMC <- group_pip.MCMC

  } 
  else {
    
    # Prior Inclusion Probabilities at the level of the variables
    var_priorprobs <- as.vector (result$priorprobs %*% models_active_vars)
    names (var_priorprobs) <- depvars
    result$var_priorprobs <- var_priorprobs

    group_priorprobs <- as.vector (result$priorprobs %*% models_active_groups)
    names (group_priorprobs) <- as.character (diff_groups)
    result$group_priorprobs <- group_priorprobs
   
  }

  model_dim <- rowSums (models_matrix) - 1 # To exclude the intercept
  result$model_dim <- model_dim

  num_active_vars <- rowSums(models_active_vars) # models_active_vars doesn´t include the intercept
  result$num_active_vars <- num_active_vars

  result$models_matrix <- models_matrix [, -1]
  result$models_active_vars <- models_active_vars
  result$models_active_groups <- models_active_groups
  result$models_w_logmarg <- models_w_logmarg[, -1] # Remove the first column (the intercept)
  result$pip <- pip

  tot.cost <- model_cost (models_active_vars, var.costs)
  result$modelcost <- tot.cost

  class(result) <- c("basglm", "bas")
  return(result)
}


# Drop the null model from Jeffrey's prior

.drop.null.bas <- function(object) {
  n.models <- object$n.models


    p <- object$size
    drop <- (1:n.models)[p == 1]
    logmarg <- object$logmarg[-drop]
    prior <- object$priorprobs[-drop]

    postprobs <- .renormalize.postprobs(logmarg, log(prior))
    which <- which.matrix(object$which[-drop], object$n.var)

    object$probne0 <- as.vector(postprobs %*% which)
    object$postprobs <- postprobs

    method <- eval(object$call$method)
    if (method == "MCMC+BAS" | method == "MCMC" | method == "MCMC_OLD") {
      object$freq <- object$freq[-drop]
      object$probne0.MCMC <- as.vector(object$freq %*% which)/sum(object$freq)
    }

    object$priorprobs <- prior
    if (!is.null(object$sampleprobs)) object$sampleprobs <- object$sampleprobs[-drop]
    object$which <- object$which[-drop]
    object$logmarg <- logmarg
    object$deviance <- object$deviance[-drop]
    object$intercept <- object$intercept[-drop]
    object$size <- object$size[-drop]
    object$Q <- object$Q[-drop]
    object$R2 <- object$R2[-drop]
    object$mle <- object$mle[-drop]
    object$mle.se <- object$mle.se[-drop]
    object$shrinkage <- object$shrinkage[-drop]
    object$n.models <- n.models - 1
    object$df <- object$df[-drop]
  return(object)
}

.renormalize.postprobs <- function(logmarg, logprior) {
  probs <- logmarg + logprior # posterior log = (marg * prior) log
  probs <- exp(probs - max(probs)) # exp(posterior log / max (posterior log))
  probs <- probs / sum(probs)
  return(probs)
}

