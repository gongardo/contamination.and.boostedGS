make_hyper_par <- function(y, X, g = n, h = c(1, 1), Z = NULL, prior_type = 1, eta = NULL, realindex = NULL, p = NULL){
  

  h<- c(1,1); cat("--\n Info: Argument h is forced to be c(1,1): prior over the model space is Scott-Berger.\n")
  g<- n; cat("--\n Info: Argument g is not used. In all cases g-Zellner with g=n is used.\n")
  cat("--\n")
	  
	
  X <- as.matrix(X)
  y <- as.vector(y)
  colnames(X) <- NULL
  colnames(y) <- NULL
  
  n <- nrow(X)
  
  if (is.null(p)) p <- ncol(X)
  
  if (is.null(Z)) {
    X <- scale(X)
    y <- scale(y, scale = FALSE)
  }
  else {
    # remove the fixed linear effects
    Z <- as.matrix(cbind(1, Z))
    y <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% y
    X <- (diag(n) - Z %*% solve(t(Z) %*% Z) %*% t(Z)) %*% X
  }
  
  if (prior_type != 2 & prior_type != 3 & prior_type != 11 & prior_type != 12 & prior_type != 13) stop("prior not supported.\n")
	  
  #prior types for PARNI and ASI
  #original:
  if (prior_type == 2) {
    g_prior_type <- "g"
    diag_V <- colSums(X^2)
  }
  #added for contamination
  if (prior_type == 3) {
    g_prior_type <- "gcont"
    diag_V <- colSums(X^2)
  }   

  #prior types for TGS and wTGS
  #fast implementation:
  if (prior_type == 12) {
    g_prior_type <- "g2TGS"
    diag_V <- colSums(X^2)
  }
  #added for contamination
  if (prior_type == 13) {
    g_prior_type <- "gcontTGS"
    diag_V <- colSums(X^2)
  }
  
  
  
  
  hyper_par <- list(n = n,
                    p = p,
                    g = g,
                    h = h,
                    X = X,
                    y = y,
                    yty = sum(y^2), 
                    ytX = t(y) %*% X,
                    diag_V = diag_V,
                    g_prior_type = g_prior_type,
					realindex = realindex,
					eta = eta)
  
  return(hyper_par)
}
