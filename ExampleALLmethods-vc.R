#virtually contaminated experiment
#exact computation, PARNI, ASI, TGS, wTGS and Gibbs sampling
#here an example with n=500 so results are expected to be
#very close to the simulated contaminated experiment
#(file ExampleALLmethods-sc.R)

rm(list=ls())

mypath<- "~/Dropbox/parallelGibbs" #gonzalo
#mypath<- "~/Dropbox/Investigacion/parallelGibbs" #maria eugenia

setwd(paste(mypath,"/R", sep=""))

library(BayesVarSel)
#functions used in the parallel implementation of Gibbs plus functions to compute inclusion probabilities theoretically
#corresponding to the contaminated experiment:
source("base-functions-ParallelGibbs.R") 
#functions to create the arguments to interact with PARNI, ASI, TGS and wTGS
#(strongly based on original PARNI "make_hyper_par.R"):
source("make_hyper_par.R")
#this contains the contaminated versions of gBF:
source("more_other_supportive_functions2.R")

#almost the original functions for PARNI and ASI:
source("ASI.R")
source("PARNI.R")
source("other_supportive_functions.R")

#functions for TGS and wTGS (almost the same as original but with fast implementation and contaminated implementation)
source("functions_for_BVS.R")


########
########
#How to use all these (simulated data)
########
########

#Simulated data:
set.seed(1234)
n<- 500; kT<- 4
X<- matrix(rnorm(n*kT), nc=kT, nr=n)
y<- .1*X[,1]+.5*X[,2]+1*X[,3]+0*X[,4]+rnorm(n)
#Note: X exactly coincide with
#the first 4 columns in X used in the file ExampleALLmethods-sc.R
#and the y's are also the same
#the rest of columns in X are independent on y
#contamination factor:
eta<- 1/sqrt(dim(X)[1]+1)
#p is the number of covariates in the contaminated experiment 
p<- 1000
#

########
########
#1. Theoretical values of inclusion probabilities
########
########

Bvs1<- Bvs(y~., data=data.frame(y=y, X=X[,1:4]), prior.models="Constant", n.keep=2^kT, prior.betas="gZellner")
C<- normConst(x= Bvs1, p=p, eta=eta)
SBFbydim<- sumBF.by.dim(x= Bvs1)
contam.incl.prob(xi=1, C=C, p=p, eta=eta, matrix.sumBF.by.dim=SBFbydim) 
#[1] 0.0008142703
contam.incl.prob(xi=2, C=C, p=p, eta=eta, matrix.sumBF.by.dim=SBFbydim) 
#[1] 1
contam.incl.prob(xi=3, C=C, p=p, eta=eta, matrix.sumBF.by.dim=SBFbydim) 
#[1] 1
contam.incl.prob(xi=4, C=C, p=p, eta=eta, matrix.sumBF.by.dim=SBFbydim) 
#[1] 0.0001635408


########
########
#2. PARNI with contaminated experiment
########
########

#we 'hide' the true variables within all these virtual covariates
#and the following is a vector that states where they are
set.seed(4321)
realindex<- sort(sample(1:p, kT, replace=FALSE))

#We make hyper par, using Scott-Berger prior, and contamination with g-prior
hyper_par <- make_hyper_par(y = scale(y, scale=FALSE),  # response
                            X = scale(X),
							prior_type = 3, # 3 is for contaminated
                            eta = eta, #contamination factor
                            realindex = realindex, #positions of real variables
                            p=p #number of total variables
							)


# PARNI
alg_par <- list(N = 6000,  # number of iterations
                Nb = 2000,  # burn-in period
                full_adap = FALSE, # whether to adapt full chain
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/p,
                target_swap = 0.25, # Parallel tempering target swapping acceptance rate
                n_chain = 25,
                n_temp = 1,
                omega_adap = "kw",omega_par = c(-1, -0.5), # Kiefer-Wolfwitz
                # omega_adap = "rm",omega_par = c(-0.7, 0.65), # Robbins-Monro
                omega_init = 0.5,
                store_chains = FALSE, # whether to store the full chains
                # bal_fun = function(x) {x/(1+x)},
                bal_fun = function(x) { min(x,1) },
                # bal_fun = sqrt,
                use_rb = TRUE  # whether to use Raoo-Blackwellised PIPs
)

results <- PARNI(hyper_par = hyper_par, alg_par = alg_par)
results$rb_PIPs[realindex]
#[1] 0.0008217392 1.0000000000 1.0000000000 0.0001635173

########
########
#3. ASI with contaminated experiment
########
########

alg_par <- list(N = 6000,
                Nb = 2000,
                full_adap = FALSE,
                verbose = TRUE,
                kappa = 0.001,
                eps = 0.1/hyper_par$p,
                phi = c(-0.7, -0.7),
                target_prop = 0.2,
                target_swap = 0.25,
                n_chain = 25,
                init_zeta = 0.5,
                n_temp = 1,
                store_chains = FALSE,
				use_rb = TRUE  # whether to use Raoo-Blackwellised PIPs
				)

results <- ASI(hyper_par = hyper_par, alg_par = alg_par)
results$rb_PIPs[realindex]
#[1] 0.0008303453 1.0000000000 1.0000000000 0.0001638375


########
########
#4. TGS with contaminated experiment (note prior_type=13)
########
########

hyper_par <- make_hyper_par(y = scale(y, scale=FALSE),  # response
                            X = scale(X),
							prior_type = 13,
                            eta = eta, #contamination factor
                            realindex = realindex, #positions of real variables
                            p=p #number of total variables
)

output_TGS <- TGS(p=p, hyper_par=hyper_par, T=50000, burn_in=5000)
output_TGS$est_inclusion_probs[realindex]
#[1] 0.0011219257 1.0000000000 1.0000000000 0.0002269558

########
########
#5. wTGS with contaminated experiment (note prior_type=13)
########
########

output_wTGS <- wTGS(p=p, hyper_par=hyper_par, T=50000, burn_in=5000)
output_wTGS$est_inclusion_probs[realindex]
#[1] 0.0008166031 1.0000000000 1.0000000000 0.0001634990

#######################
#from now on the hyper_par argument is not used
#######################

########
########
#6. Gibbs sampling (component-by-component) with contaminated experiment
########
########

d.ini<- 5 #dimension of initial model (which is randomly chosen)
N<- 1000 #number of iterations
source("contamination-Gibbs.R")
incl.probRB[realindex]
#[1] 0.0008282635 1.0000000000 1.0000000000 0.0001632854

########
########
#7. Gibbs sampling (vectorized version) with contaminated experiment
########
########

#not implemented (how it works, it does not make much sense for the contaminated experiment)

########
########
#8. Gibbs sampling (parallelized version of component-by-component) with contaminated experiment
########
########

#not implemented (how it works, it does not make much sense for the contaminated experiment)



