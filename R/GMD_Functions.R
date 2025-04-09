#####GMD Functions#####
require(mclust)
require(mvtnorm)
require(parallel)
require(doParallel)

#Note that all functions herein that take functional data as input assume the shape (number of time points) x (number of functions)

#' The following function will be used to calculate the L2 inner product of a function
#' 
#' @param input A numerical vector of increasing value where the function(s) were evaluated
#' @param values A numerical vector consisting of the function evaluations at the corresponding input
#' @return A numerical value of the inner product

trap_integration <- function(input, values){
  
  # Sanity check that confirms the length of input == length values
  stopifnot(length(input) == length(values))
  # Sanity check to confirm that input is an increasing vector
  stopifnot( all( diff(input) > 0 ) )
  
  # Use the trapezoidal rule to approximate the integral at each point in the domain
  return( sum( 0.5 * (values[-1] + values[-length(values)]) * diff(input) ) )
  
}

#' The following function calculates the inner product between each function and each basis function
#' without explicitly looping over the number of functions
#' @param input a numeric vector of increasing time points; it is where each observed function is evaluated
#' @param x a numeric matrix where each column is an observed function
#' @param phi a numeric matrix where each column is one of the orthogonal basis functions
#' @return a vector consisting of the n*p inner product values; the order will go (x[,1],phi[,1]), (x[,2],phi[,1]), ..., (x[,n],phi[,1]), ... ..., (x[,1],phi[,p]), ..., (x[,n],phi[,p])

xiFind <- function(x, phi, input){
  
  # Sanity check to confirm the dimensions match
  stopifnot( length(input) == dim(x)[1] & length(input) == dim(phi)[1] )
  # Sanity check to confirm that input is an increasing vector
  stopifnot( all( diff(input) > 0 ) )
  
  # Universal constants
  p <- dim(phi)[2];    # p = number of basis functions
  n <- dim(x)[2];      # n = number of functions observed
  tr <- length(input); # tr = number of time points in the domain
  
  # Create a matrix of dimension tr-by-(p*n); each column will be integrated
  temp <- matrix(data = rep(x, times = p), nrow = tr, ncol = n*p) *
    matrix(data = rep(t(phi), each = n), nrow = tr, ncol = n*p, byrow = TRUE)
  
  # Integrate each column of the matrix
  temp <- apply(X = temp, MARGIN = 2,
                FUN = function(x, input){trap_integration(input = input, values = x)},
                input = input)
  
  return(temp)
}

# Since the methods rely upon estimating the parameters of a Gaussian mixture
# model, the code will use mclust::Mclust to estimate those parameters. 
# However, errors and warnings can occur using Mclust, so we need to catch them 
# and let the user know to adjust the number of clusters or the number of 
# parameters (by enforcing a common covariance structure)

myTryCatch <- function(expr){
  warn <- err <- NULL
  value <- withCallingHandlers(
    tryCatch(expr, error = function(e) {
      err <<- e
      NULL
    }), warning = function(w) {
      warn <<- w
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = warn, error = err)
}

#'
#' Create a function which will automate conformal prediction on functional data using
#' the methods presend in Lei, Rinaldo, Wasserman (2015) "A conformal prediction approach to explore functional data"
#' This first function will only produce the NonConformity Scores for the calibration data
#' A second function will take the output from the first and create the "band" using the testing data.
#' A third function will take the output from the first and create the p-values for the testing data.
#' 
#' @param input a numeric vector of increasing points; this vector represents the domain upon which each device is evaluated
#' @param train a matrix of size length(input)-by-n_train; each column represents one device in the training set
#' @param cal a matrix of size length(input)-by-n_cal; each column represents one device in the calibration set
#' @param p a non-negative integer that represents the number of orthogonal basis functions to project the device into
#' @param K a non-negative integer that represents the number of clusters (i.e. the number of Gaussian processes used to construct each device)
#' @return a list
#' 

functionalCP <- function(input, train, cal, p, K){
  
  # Sanity checks
  stopifnot( length(input) == dim(train)[1] & length(input) == dim(cal)[1] ) # Confirms that each function is measured at the same number of points
  stopifnot( all(diff(input) >= 0 ) ) # Confirms that the independent variable is increasing
  stopifnot( length(p) == 1L & is.numeric(p) & p >= 1L ) # Confirms that the number of principal components is at least 1
  stopifnot( length(K) == 1L & is.numeric(K) & K >= 1L ) # Confirms that the number of clusters is at least one
  
  # Universal constants
  n_train <- dim(train)[2] # n_train = number of devices in training set
  n_cal <- dim(cal)[2] # n_cal = number of devices in calibration set
  tp <- length(input) # tp = number of points in the domain
  
  # More sanity checks
  stopifnot(p < tp)
  
  # Step 1: Calculate the orthonormal basis functions using the training data
  Sigma <- train %*% t(train) - n_train * rowMeans(train) %*% t( rowMeans(train) )
  basisFns <- eigen(Sigma)$vectors
  
  # Normalize the basis functions
  basisFns <- apply(X = basisFns, MARGIN = 2,
                    FUN = function(x, input){x/sqrt(trap_integration(input = input, values = x^2))},
                    input = input)
  
  # Step 2: Keep only the first p orthonormal basis functions
  phi <- basisFns[,1:p]
  
  # Step 3: Find the principal component scores, xi, for the training and calibration devices
  xi_train <- matrix(data = xiFind(x = train, phi = phi, input = input),
                     nrow = p, ncol = n_train, byrow = TRUE)
  xi_cal <- matrix(data = xiFind(x = cal, phi = phi, input = input),
                   nrow = p, ncol = n_cal, byrow = TRUE)
  
  # Step 4: Find the parameter estimates for each of the K clusters
  estimates <- myTryCatch( 
    mclust::densityMclust( data = t(xi_train), G = K, modelNames = "VVV", warn = TRUE, verbose = FALSE, plot = FALSE) )
  
  # Step 4b: If K == 1 then reconfigure the data to be the correct dimensions
  if(K == 1){
    estimates$value$parameters$variance$cholsigma <- array(data = estimates$value$parameters$variance$cholsigma, dim = c(p,p,1))
    estimates$value$parameters$variance$sigma <- array(data = estimates$value$parameters$variance$sigma, dim = c(p, p, 1))
  }
  
  # Step 5: Find the nonconformity measure (the negative density) for each calibration device
  if(!is.null(estimates$warning)){stop("Error: Not enough data for fitting the multivariate Gaussians. Consider lowering the number of clusters, changing the covariance structure, or increasing the number of devices in the training set.")}
  nc_scores <- -1 * mclust::predict.densityMclust(object = estimates$value, newdata = t(xi_cal), what = "dens")
  
  # Step 6: Return the values
  output <- list(basisFns = phi,
                 ncm_scores = nc_scores,
                 densMclust = estimates$value,
                 input = input,
                 train = train,
                 cal = cal)
  class(output) <- append("GaussianMixture", class(output))
  
  return(output)
}

#'
#' This function takes the output from the functionalCP function and produces the p-values for devices
#' in the testing set.
#' @param Gmix an object of class GaussianMixture
#' @param test a matrix of size length(input)-by-n_test
#' @return a vector of p-values for something
#' 

functionalCPpvals <- function(Gmix, test){
  
  # Sanity Check
  stopifnot(is(Gmix, "GaussianMixture"))
  
  # Acquire useful data from Gmix
  train <- Gmix$train
  n_train <- dim(train)[2]
  cal <- Gmix$cal
  n_cal <- dim(cal)[2]
  input <- Gmix$input
  tp <- length(input)
  p <- Gmix$densMclust$d
  K <- Gmix$densMclust$G
  phi <- Gmix$basisFns
  n_test <- dim(test)[2]
  ncm_scores <- Gmix$ncm_scores
  
  # Additional Sanity checks
  stopifnot(dim(test)[1] == tp)
  
  # Step 1: Convert the test data into principal component scores
  xi_test <- matrix(data = xiFind(x = test, phi = phi, input = input),
                    nrow = p, ncol = n_test, byrow = TRUE)
  
  # Step 2: Calculate the nonconformity measure for each cluster for each device in test set
  nc_scores <- -1 * mclust::predict.densityMclust(object = Gmix$densMclust, newdata = t(xi_test), what = "dens")
  
  # Step 3: Calculate the p-values for each device in the test set
  pvalues <- lapply(X = seq_along(nc_scores),
                    FUN = function(x, ncm, nc){
                      return( ( sum(as.numeric(ncm > nc[x])) + runif(n = 1)*(1 + sum(as.numeric(ncm == nc[x]) ) ) )/( n_cal + 1 ) ) },
                    ncm = ncm_scores, nc = nc_scores)
  pvalues <- unlist(pvalues)
  
  # Return the p-values
  return(pvalues)
  
}


#'
#' Like the previous function, functionalCP, this function is the first step
#' in producing conformal prediction bounds or p-values for a set of functions
#' using the third method presented in Lei, Rinaldo, Wasserman (2015) "A conformal prediction approach to explore functional data"
#' This function will only produce the NonConformity Scores for the calibration data
#' using the maximum component in the mixture. An additional function will take
#' the output from this function to produce a "band" and a third function will take
#' this output and a set of test functions to produce a p-value for each test function.
#' 
#' @param input a numeric vector of increasing points; this vector represents the domain upon which each device is evaluated
#' @param train a matrix of size length(input)-by-n_train; each column represents one device in the training set
#' @param cal a matrix of size length(input)-by-n_cal; each column represents one device in the calibration set
#' @param p a non-negative integer that represents the number of orthogonal basis functions to project the device into
#' @param K a non-negative integer that represents the number of clusters (i.e. the number of Gaussian processes used to construct each device)
#' @return a list
#' 

functionalCP_m <- function(input, train, cal, p, K){
  
  # Sanity checks
  stopifnot( length(input) == dim(train)[1] & length(input) == dim(cal)[1] ) # Confirms that each function is measured at the same number of points
  stopifnot( all(diff(input) >= 0 ) ) # Confirms that the independent variable is increasing
  stopifnot( length(p) == 1L & is.numeric(p) & p >= 1L ) # Confirms that the number of principal components is at least 1
  stopifnot( length(K) == 1L & is.numeric(K) & K >= 1L ) # Confirms that the number of clusters is at least one
  
  # Universal constants
  n_train <- dim(train)[2] # n_train = number of devices in training set
  n_cal <- dim(cal)[2] # n_cal = number of devices in calibration set
  tp <- length(input) # tp = number of points in the domain
  
  # More sanity checks
  stopifnot(p < tp)
  
  # Step 1: Calculate the orthonormal basis functions using the training data
  Sigma <- train %*% t(train) - n_train * rowMeans(train) %*% t( rowMeans(train) )
  basisFns <- eigen(Sigma)$vectors
  
  # Normalize the basis functions
  basisFns <- apply(X = basisFns, MARGIN = 2,
                    FUN = function(x, input){x/sqrt(trap_integration(input = input, values = x^2))},
                    input = input)
  
  # Step 2: Keep only the first p orthonormal basis functions
  phi <- basisFns[,1:p]
  
  # Step 3: Find the principal component scores, xi, for the training and calibration devices
  xi_train <- matrix(data = xiFind(x = train, phi = phi, input = input),
                     nrow = p, ncol = n_train, byrow = TRUE)
  xi_cal <- matrix(data = xiFind(x = cal, phi = phi, input = input),
                   nrow = p, ncol = n_cal, byrow = TRUE)
  
  # Step 4: Find the parameter estimates for each of the K clusters
  estimates <- myTryCatch( 
    mclust::densityMclust( data = t(xi_train), G = K, modelNames = "VVV", warn = TRUE, verbose = FALSE, plot = FALSE) )

  # Step 4b: If K == 1 then reconfigure the data to be the correct dimensions
  if(K == 1){
    estimates$value$parameters$variance$cholsigma <- array(data = estimates$value$parameters$variance$cholsigma, dim = c(p,p,1))
    estimates$value$parameters$variance$sigma <- array(data = estimates$value$parameters$variance$sigma, dim = c(p, p, 1))
  }
  
  # Step 5: Estimate the density of the each function in the calibration set for each of the K clusters
  if(!is.null(estimates$warning)){stop("Error: Not enough data for fitting the multivariate Gaussians. Consider lowering the number of clusters, changing the covariance structure, or increasing the number of devices in the training set.")}
  ncm_scores <- mclust::predict.densityMclust(object = estimates$value, newdata = t(xi_cal), what = 'cdens')
  ncm_scores <- ncm_scores * matrix(data = rep(estimates$value$parameters$pro, n_cal), nrow = n_cal, ncol = K, byrow = TRUE)
  ncm_scores <- -1 * as.numeric( apply(ncm_scores, 1, max) )
    
  # Step 6: Return the values
  output <- list(basisFns = phi,
                 ncm_scores = ncm_scores,
                 densMclust = estimates$value,
                 input = input,
                 train = train,
                 cal = cal)
  class(output) <- append("GaussianMixture", class(output))
  
  return(output)
}

#'
#' This function takes the output from the functionalCP_m function and produces the p-values for devices
#' in the testing set.
#' @param Gmix an object of class GaussianMixture
#' @param test a matrix of size length(input)-by-n_test
#' @return a vector of p-values for something
#' 

functionalCPpvals_m <- function(Gmix, test){
  
  # Sanity Check
  stopifnot(is(Gmix, "GaussianMixture"))
  
  # Acquire useful data from Gmix
  train <- Gmix$train
  n_train <- dim(train)[2]
  cal <- Gmix$cal
  n_cal <- dim(cal)[2]
  input <- Gmix$input
  tp <- length(input)
  p <- Gmix$densMclust$d
  K <- Gmix$densMclust$G
  phi <- Gmix$basisFns
  n_test <- dim(test)[2]
  ncm_scores <- Gmix$ncm_scores
  
  # Additional Sanity checks
  stopifnot(dim(test)[1] == tp)
  
  # Step 1: Convert the test data into principal component scores
  xi_test <- matrix(data = xiFind(x = test, phi = phi, input = input),
                    nrow = p, ncol = n_test, byrow = TRUE)
  
  # Step 2: Calculate the nonconformity measure for each cluster and for each device (function)
  # in the test set
  nc_scores <- mclust::predict.densityMclust(object = Gmix$densMclust, newdata = t(xi_test), what = 'cdens')
  nc_scores <- nc_scores * matrix(data = rep(Gmix$densMclust$parameters$pro, n_test), nrow = n_test, ncol = K, byrow = T )
  nc_scores <- -1 * as.numeric( apply(nc_scores, 1, max) )
  
  # Step 3: Calculate the p-values for each device in the test set
  pvalues <- lapply(X = seq_along(nc_scores),
                    FUN = function(x, ncm, nc){
                      return( ( sum(as.numeric(ncm > nc[x])) + runif(n = 1)*(1 + sum(as.numeric(ncm == nc[x]) ) ) )/( n_cal + 1 ) ) },
                    ncm = ncm_scores, nc = nc_scores)
  pvalues <- unlist(pvalues)
  
  # Return the p-values
  return(pvalues)
  
}


