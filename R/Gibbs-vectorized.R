#before file called GibbsPARNI-light.R
resumen<- list()
resumen$params<- c(N)
names(resumen$params)<- c("N")


incl.probRB<- rep(0,p)
names(incl.probRB)<- colnames(X)

current.model<- sample(c(rep(1, d.ini), rep(0, p-d.ini)))

visited.models<- list() #positions with the active variables in each step

pb <- txtProgressBar(min = 0, max = N, style = 3)

resumen$time.step2<- system.time(
  for (i in 1:N){
    #cat("Iteration:", i, "\n")
	setTxtProgressBar(pb, i)

	current.component<- 1
	while (current.component < p){
		current.kgamma<- sum(current.model)
		lB10<- bf_g_L_4GS_log(current.model, hyper_par, p.ini = current.component)
		lprior10<- lchoose(p, current.kgamma-current.model[current.component:p])-lchoose(p, current.kgamma+1-current.model[current.component:p])
		ratio<- 1/(1 + 1/(exp(lB10+lprior10)))
		proposal.model<- rbinom(n=p-current.component+1, prob=pmin(1,ratio), size=1)
		
		condition<- proposal.model==current.model[current.component:p]
		#sum(condition)
		
		if (sum(!condition) == 0){
			incl.probRB[current.component:p]<- incl.probRB[current.component:p]+ratio
			current.component<- p
			
		}
		if (sum(!condition) != 0){
			valid<- min(which(condition == FALSE))
			new.component<- current.component + valid -1
			current.model[new.component]<- 1 - current.model[new.component] #=proposal.model[valid]
			incl.probRB[current.component:new.component]<- incl.probRB[current.component:new.component]+ratio[1:valid]
			current.component<- new.component+1
			
		}
	}
	visited.models[[i]]<- which(current.model==1)
}
)

cat("\n")

incl.probRB<- incl.probRB/N
resumen$inclprob.step2 <- incl.probRB



	
