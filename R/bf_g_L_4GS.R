bf_g_L_4GS <- function(gamma, hyper_par, p.ini = 1){
  #this function is the adaptation of bf_g_L by Griffin et al () to be used in Gibbs
  #the new p.ini argument is from which component we want to compute the ratio
  #the result is a p-p.ini+1 vector with component 1:
  #BF(gamma_1=gamma_1,..,gamma_{p.ini}=1,gamma_p=gamma_p)/BF(gamma_1=gamma_1,..,gamma_{p.ini}=0,gamma_p=gamma_p)
  #component 2:       #BF(gamma_1=gamma_1,..,gamma_{p.ini}=gamma_{p.ini},gamma_{p.ini+1}=1,gamma_p=gamma_p)/BF(gamma_1=gamma_1,..,gamma_{p.ini}=gamma_{p.ini},gamma_{p.ini+1}=0,gamma_p=gamma_p)
  #and so on...
  	
  includes <- which(gamma == 1)
  #includes <- includes[includes >= p.ini]
    
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  # XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  X <- hyper_par$X
  diag_V <- hyper_par$diag_V
  # y <- hyper_par$y
  # h <- hyper_par$h
  
  g_ratio <- g/(g+1)
  inv_sqrtg1 <- 1/sqrt(g+1)
  
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)
  
  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in p.ini:p){
      tilda_A <- yty - g_ratio * ytX[j]^2/diag_V[j]
      
      A_ratio <- A/tilda_A
      
      BF[j] <- (A_ratio)^n_power * inv_sqrtg1
    }
  }
  else{
    
	zj <- sum(includes < p.ini)
    
    Xg <- X[, includes]
    XgtX <- t(Xg) %*% X
    XgtXg <- XgtX[, includes]  # XtX[includes, includes]
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(XgtXg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((L_F %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty * g_ratio
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio
    
    # print(( d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio)[1])
    
    BF <- (A / tilda_A_vec)^n_power * inv_sqrtg1
    
    for (j in includes[includes >= p.ini]) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
      }
      else{
        
        zj <- zj+1
        
        # dj <- 1/F[zj,zj]
        # tilde_A <- A + (ytXg %*% F[,zj] )^2 * g_ratio * dj
        # A_ratio <- tilde_A/A
        
        A_ratio <- 1 + (ytXg %*% F[,zj] )^2 * g_ratio / (A * F[zj,zj])
        
      }
      
      
      BF[j] <- A_ratio^n_power * inv_sqrtg1
      
    }
  }
  
  return(BF[p.ini:p])
}

bf_g_L_4GS_log <- function(gamma, hyper_par, p.ini = 1){
	#same as above but returning log
  #this function is the adaptation of bf_g_L by Griffin et al () to be used in Gibbs
  #the new p.ini argument is from which component we want to compute the ratio
  #the result is a p-p.ini+1 vector with component 1:
  #BF(gamma_1=gamma_1,..,gamma_{p.ini}=1,gamma_p=gamma_p)/BF(gamma_1=gamma_1,..,gamma_{p.ini}=0,gamma_p=gamma_p)
  #component 2:       #BF(gamma_1=gamma_1,..,gamma_{p.ini}=gamma_{p.ini},gamma_{p.ini+1}=1,gamma_p=gamma_p)/BF(gamma_1=gamma_1,..,gamma_{p.ini}=gamma_{p.ini},gamma_{p.ini+1}=0,gamma_p=gamma_p)
  #and so on...
  	
  includes <- which(gamma == 1)
  #includes <- includes[includes >= p.ini]
    
  n <- hyper_par$n
  p <- hyper_par$p
  yty <- hyper_par$yty
  # XtX <- hyper_par$XtX
  ytX <- hyper_par$ytX
  g <- hyper_par$g
  X <- hyper_par$X
  diag_V <- hyper_par$diag_V
  # y <- hyper_par$y
  # h <- hyper_par$h
  
  g_ratio <- g/(g+1)
  inv_sqrtg1 <- 1/sqrt(g+1)
  
  n_power <- n/2
  
  BF <- rep(NA, p)
  p_gam <- sum(gamma)
  
  # empty model
  if(p_gam == 0){
    A <- yty
    for(j in p.ini:p){
      tilda_A <- yty - g_ratio * ytX[j]^2/diag_V[j]
      
      A_ratio <- A/tilda_A
      
      BF[j] <- n_power*log(A_ratio) + log(inv_sqrtg1) #(A_ratio)^n_power * inv_sqrtg1
    }
  }
  else{
    
	zj <- sum(includes < p.ini)
    
    Xg <- X[, includes]
    XgtX <- t(Xg) %*% X
    XgtXg <- XgtX[, includes]  # XtX[includes, includes]
    
    ytXg <- matrix(ytX[includes],nrow = 1,ncol = p_gam)
    
    F <- solve(XgtXg)
    L_F <- chol(F)
    
    ytXgFXgty <- sum((L_F %*% t(ytXg))^2)
    
    A <- yty - ytXgFXgty * g_ratio
    
    # XgtX <- matrix(XtX[includes,], nrow = p_gam, ncol = p)
    # XgtX <- XtX[includes,]
    
    d_vec <- 1 / (diag_V - rowSums((t(XgtX) %*% t(L_F))^2))
    
    ytXFXtxj_vec <- ytXg %*% F %*% XgtX
    tilda_A_vec <- A - d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio
    
    # print(( d_vec * (ytXFXtxj_vec - ytX)^2 * g_ratio)[1])
    these.indexes<- which((1:p) >= p.ini & gamma!=1)
    BF[these.indexes] <- n_power*(log(A) - log(tilda_A_vec[these.indexes])) + log(inv_sqrtg1) #(A / tilda_A_vec)^n_power * inv_sqrtg1
    
    for (j in includes[includes >= p.ini]) {
      
      if(p_gam == 1){
        
        tilde_A <- yty
        A_ratio <- tilde_A/A
        
      }
      else{
        
        zj <- zj+1
        
        # dj <- 1/F[zj,zj]
        # tilde_A <- A + (ytXg %*% F[,zj] )^2 * g_ratio * dj
        # A_ratio <- tilde_A/A
        
        A_ratio <- 1 + (ytXg %*% F[,zj] )^2 * g_ratio / (A * F[zj,zj])
        
      }
      
      
      BF[j] <- n_power*log(A_ratio) + log(inv_sqrtg1) #A_ratio^n_power * inv_sqrtg1
      
    }
  }
  
  return(BF[p.ini:p])
}

