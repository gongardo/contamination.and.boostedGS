bf_gcont2_L <- function(gamma, hyper_par){
	#fast implementation of bf_cont_L based on their bf_g_L
	#this function returns a p=length(gamma) vector, where in the i-th position has the ratio 
#	#BF(gamma_1=gamma_1,..,gamma_i=1,gamma_p=gamma_p)/BF(gamma_1=gamma_1,..,gamma_i=0,gamma_p=gamma_p)
#	
#	#We can do so, with the following code
	gamma<- 1*gamma
	result<- rep(hyper_par$eta, hyper_par$p)
	
	result[hyper_par$realindex]<- bf_g_L(gamma[hyper_par$realindex], hyper_par)
	
	return(matrix(result, nr=1))

}

log_llh_gcont2_L <- function(gamma, hyper_par){
	#this function just returns the log(BF) of gamma to the null plus log(marginal(null))
	gamma<- 1*gamma
	log.m0<- -hyper_par$n*log(hyper_par$yty)/2

	#number of spureous
	nofesp<- sum(gamma[-hyper_par$realindex])
	if (sum(gamma[hyper_par$realindex])==0) return(nofesp*log(hyper_par$eta)+log.m0)	
	else	
		return(log_llh_g_L(gamma[hyper_par$realindex], hyper_par)+nofesp*log(hyper_par$eta))	
}
