#### MAIN FUNCTIONS FOR SAMPLERS ######
###
GS<-function(p,gamma_start=NULL,
             T,burn_in=0,thin=1,
             hyper_par=NULL,vars_selected=c(1,2)){ 
  ## random scan Gibbs Sampler for Bayesian variable selection problems
  
  stop("not yet implemented\n")    	  
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform random scan Gibbs Sampling step
    for (iter in 1:thin){
      j<-sample.int(n = p,size = 1)
      fc_j<-single_full_cond(j=j,gamma=gamma,hyper_par=hyper_par)$fc_j
      # ap<-1-fc_j## not metropolized
      ap<-(1-fc_j)/fc_j## metropolized
      if(runif(1) < ap){
        gamma[j]<-1-gamma[j]
        indices_sequence[t]<-j
      }
    }
    ## update inclusion estimators
    if(t>burn_in){
      est_inclusion_probs<-est_inclusion_probs+gamma/T
    }
    ## store indices sequence for output analysis
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma))
}

TGS<-function(p,gamma_start=NULL,
              T,burn_in=0,
              hyper_par=NULL,vars_selected=c(1,2)){
  ## TGS algorithm for Bayesian variable selection problems
  
  #code added to allow for the fast implementation and the contaminated implementation
  if (hyper_par$g_prior_type != "g2TGS" & hyper_par$g_prior_type != "gcontTGS") stop("Prior type not supported")
  if (hyper_par$g_prior_type == "g2TGS") {single_full_cond<- single_full_cond2; full_cond<- full_cond2}
  if (hyper_par$g_prior_type == "gcontTGS") {single_full_cond<- single_full_cond_cont; full_cond<- full_cond_cont}
  #end of code added
  
  #code added to save visited models
  visited.models<- list()
  #end of code added
  
    	    
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  est_inclusion_probs_1<-rep(NA,T)
  est_inclusion_probs_2<-rep(NA,T)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  sample_weights<-rep(NA,T)
  ## compute full conditionals
  output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform Tempered Gibbs Sampling step
    if(length(which(output_fc$fc==0))>0){
      j<-sample(which(output_fc$fc==0),size = 1)
    }else{
      j<-which(cumsum(1/output_fc$fc)>(runif(1)*sum(1/output_fc$fc)))[1]
    }
    gamma[j]<-1-gamma[j]
    # gamma[j]<-rbinom(1,1,0.5)
    ## compute full conditionals and update Rao-Blackwelied estimator
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=output_fc$stored)
    if(t>burn_in){
      sample_weights[t-burn_in]<-(p/2)/sum(1/output_fc$fc)
      est_inclusion_probs<-est_inclusion_probs+sample_weights[t-burn_in]*
        (gamma*output_fc$fc+(1-gamma)*(1-output_fc$fc))
      est_inclusion_probs_1[t-burn_in]<-est_inclusion_probs[vars_selected[1]]
      est_inclusion_probs_2[t-burn_in]<-est_inclusion_probs[vars_selected[2]]
    }
    ## store indices sequence for output analysis
    indices_sequence[t]<-j
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
	
    #code added to save visited models
	visited.models[[t]]<- which(gamma==1)
	#end of code added
	
  }
  est_inclusion_probs<-est_inclusion_probs/sum(sample_weights)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              sample_weights=sample_weights,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2,
			  visited.models=visited.models))
}

wTGS<-function(p,gamma_start=NULL,
              T,burn_in=0,
              hyper_par=NULL,vars_selected=c(1,2)){
  ## wTGS algorithm for Bayesian variable selection problems
  
  #code added to allow for the fast implementation and the contaminated implementation
  if (hyper_par$g_prior_type != "g2TGS" & hyper_par$g_prior_type != "gcontTGS") stop("Prior type not supported")
  if (hyper_par$g_prior_type == "g2TGS") {single_full_cond<- single_full_cond2; full_cond<- full_cond2}
  if (hyper_par$g_prior_type == "gcontTGS") {single_full_cond<- single_full_cond_cont; full_cond<- full_cond_cont}
  #end of code added
  
  #code added to save visited models
  visited.models<- list()
  #end of code added
  
    	      
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  est_inclusion_probs_1<-rep(NA,T)
  est_inclusion_probs_2<-rep(NA,T)
  indices_sequence<-rep(NA,burn_in+T)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  sample_weights<-rep(NA,T)
  ## compute full conditionals
  output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
  flip_rates<-compute_flip_rates(output_fc=output_fc,gamma=gamma)
 
  ## MCMC iteration
  for (t in 1:(T+burn_in)){ #
   ## perform Markov chain transition step
   if(length(which(output_fc$fc==0))>0){
     if(length(which(output_fc$fc==0))==1){ j<-which(output_fc$fc==0)}
      else{j<-sample(as.vector(which(output_fc$fc==0)),size = 1)} #Seems that with a unique value it samples from 1 to this value
    }else{
      j<-which(cumsum(flip_rates)>(runif(1)*sum(flip_rates)))[1]
    }
    gamma[j]<-1-gamma[j]
    ## compute full conditionals and update Rao-Blackwelied estimator
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=output_fc$stored)
    flip_rates<-compute_flip_rates(output_fc=output_fc,gamma=gamma)
    # flip_rates[which(gamma==0)]<-(1-output_fc$fc[which(gamma==0)])/output_fc$fc[which(gamma==0)]
    # flip_rates[which(gamma==1)]<-1
    if(t>burn_in){
      sample_weights[t-burn_in]<-(p/2)/sum(flip_rates)
      est_inclusion_probs<-est_inclusion_probs+sample_weights[t-burn_in]*
        (gamma*output_fc$fc+(1-gamma)*(1-output_fc$fc))
      est_inclusion_probs_1[t-burn_in]<-est_inclusion_probs[vars_selected[1]]
      est_inclusion_probs_2[t-burn_in]<-est_inclusion_probs[vars_selected[2]]
    }
    ## store indices sequence for output analysis
    indices_sequence[t]<-j
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
	
    #code added to save visited models
	visited.models[[t]]<- which(gamma==1)
	#end of code added
	
  }
  est_inclusion_probs<-est_inclusion_probs/sum(sample_weights)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              sample_weights=sample_weights,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2,
			  visited.models=visited.models))
}


GS_RB<-function(p,gamma_start=NULL,
               T,burn_in=0,thin=1,
               hyper_par=NULL,vars_selected=c(1,2)){
  ## deterministic scan Gibbs Sampler with Rao-Blackwellization for Bayesian variable selection problems
  
  stop("not yet implemented\n")
  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    for (j in 1:p){
      fc_j<-single_full_cond(j=j,gamma=gamma,hyper_par=hyper_par)$fc_j
      # ap<-1-fc_j## not metropolized
      ap<-(1-fc_j)/fc_j## metropolized
      ## update inclusion estimators
      if(t>burn_in){
        est_inclusion_probs[j]<-est_inclusion_probs[j]+
          (gamma[j]*fc_j+(1-gamma[j])*(1-fc_j))/T
      }
      if(runif(1) < ap){
        gamma[j]<-1-gamma[j]
      }
    }
    ## store indices sequence for output analysis
    gamma_1[t]<-gamma[vars_selected[1]]
    gamma_2[t]<-gamma[vars_selected[2]]
  }
  return(list(est_inclusion_probs=est_inclusion_probs,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma))
}

HBS<-function(p,gamma_start=NULL,full_cond,
              T,burn_in=0,
              hyper_par=NULL){
  ## Hamming Ball Sampler for Bayesian variable selection problems
  stop("not yet implemented\n")

  ## initialize vector of inclusion indicators
  if(is.null(gamma_start)){
    gamma_start<-rep(0,p)
  }
  gamma<-gamma_start
  ## initialize vector to be output
  est_inclusion_probs<-rep(0,p)
  indices_sequence<-matrix(rep(NA,2*(burn_in+T)),nrow=2)
  gamma_1<-rep(NA,burn_in+T)
  gamma_2<-rep(NA,burn_in+T)
  ## MCMC iteration
  for (t in 1:(T+burn_in)){
    ## perform Hamming ball sampler move (i.e. sample U|X and X|U)
    j_aux<-sample.int(p,1)
    gamma[j_aux]<-1-gamma[j_aux]
    output_fc<-full_cond(gamma=gamma,hyper_par=hyper_par,stored=NULL)
    probs<-(1-output_fc$fc)/output_fc$fc
    if(length(which(probs==Inf))>0){
      j<-sample(which(probs==Inf),size = 1)
    }else{
      j<-which(cumsum(probs)>(runif(1)*sum(probs)))[1]
    }
    gamma[j]<-1-gamma[j]
    ## update inclusion estimators
    if(t>burn_in){
      est_inclusion_probs<-est_inclusion_probs+gamma/T
    }
    ## store indices sequence for output analysis
    indices_sequence[,t]<-c(j_aux,j)
    gamma_1[t]<-gamma[1]
    gamma_2[t]<-gamma[2]
  }
  est_inclusion_probs_1<-cumsum(gamma_1[(burn_in+1):(burn_in+T)])/c(1:T)
  est_inclusion_probs_2<-cumsum(gamma_2[(burn_in+1):(burn_in+T)])/c(1:T)
  return(list(est_inclusion_probs=est_inclusion_probs,
              indices_sequence=indices_sequence,
              gamma_1=gamma_1,
              gamma_2=gamma_2,
              final_gamma=gamma,
              est_inclusion_probs_1=est_inclusion_probs_1,
              est_inclusion_probs_2=est_inclusion_probs_2))
}


#single_full_cond contaminated version (23-11-23)...see Notes.pdf
single_full_cond_cont<-function(j, gamma, hyper_par, stored=NULL){
	gamma.1<- gamma; gamma.1[j]<- 1
	gamma.0<- gamma; gamma.0[j]<- 0
	B01<- exp(log_llh_gcont2_L(gamma.0, hyper_par)-log_llh_gcont2_L(gamma.1, hyper_par))	
	prob<- 1/(1+B01*exp(lchoose(hyper_par$p, sum(gamma)-gamma[j]+1)-lchoose(hyper_par$p, sum(gamma)-gamma[j])))
	output<-list(fc_j=prob^gamma[j]*(1-prob)^(1-gamma[j]), stored=NULL)
    return(output)
	
}

#full_cond contaminated version (23-11-23)...see Notes.pdf
full_cond_cont<- function(gamma, hyper_par, stored=NULL){
	prob<- 1/(1+(1/bf_gcont2_L(gamma, hyper_par))*exp(lchoose(hyper_par$p, sum(gamma)-gamma+1)-lchoose(hyper_par$p, sum(gamma)-gamma)))
    output<-list(fc=prob^gamma*(1-prob)^(1-gamma), stored=NULL)
    return(output)
	
}

#single_full_cond fast version using code from PARNI
single_full_cond2<-function(j,gamma,hyper_par=NULL,stored=NULL){
  gamma.1<- gamma; gamma.1[j]<- 1
  gamma.0<- gamma; gamma.0[j]<- 0
  B01<- exp(log_llh_g_L(gamma.0, hyper_par)-log_llh_g_L(gamma.1, hyper_par))	
  prob<- 1/(1+B01*exp(lchoose(hyper_par$p, sum(gamma)-gamma[j]+1)-lchoose(hyper_par$p, sum(gamma)-gamma[j])))
  output<-list(fc_j=prob^gamma[j]*(1-prob)^(1-gamma[j]), stored=NULL) 
  return(output)
}

#full_cond fast version using code from PARNI
full_cond2<-function(gamma,hyper_par=NULL,stored=NULL){
  prob<- 1/(1+(1/bf_g_L(gamma, hyper_par))*exp(lchoose(hyper_par$p, sum(gamma)-gamma+1)-lchoose(hyper_par$p, sum(gamma)-gamma)))
  output<-list(fc=prob^gamma*(1-prob)^(1-gamma), stored=NULL)
  return(output)
}

compute_flip_rates<-function(output_fc,gamma){
  # TGS version
  flip_rates<-rep(NA,length(gamma))
  flip_rates[which(gamma==0)]<-(1-output_fc$fc[which(gamma==0)])/output_fc$fc[which(gamma==0)]
  flip_rates[which(gamma==1)]<-1
  return(flip_rates)
}

