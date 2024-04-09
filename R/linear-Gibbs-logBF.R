#rm(incl.probRB)
resumen.gibbs<- list()
resumen.gibbs$params<- c(N)
names(resumen.gibbs$params)<- c("N")

write.output<- TRUE

incl.probRB<- rep(0,p)
current.model<- sample(c(rep(1, d.ini), rep(0, p-d.ini)))
lB.PMcurrent<- lB.PM(current.model,y,X)
proposal.model<- current.model
pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                     max = N, # Valor máximo de la barra de progreso
                     style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                     width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                     char = "=")   # Caracter usado para crear la barra


#Rao-Blackwellized inclusion probabilities:
names(incl.probRB)<- colnames(X)
incl.probRB<- rep(0,p)

incl.probRB<- 0*incl.probRB
resumen.gibbs$time<- system.time(
  for (i in 1:N){
	  for (j in 1:p){
      proposal.model<- current.model; proposal.model[j]<- 1-current.model[j]
      lB.PMproposal<- lB.PM(proposal.model,y,X)
      ratio<- 1/(1+exp(lB.PMcurrent-lB.PMproposal))
      
      if (runif(1)<ratio) {
        current.model[j]<- proposal.model[j]; lB.PMcurrent<- lB.PMproposal
      }
      #rnd<- runif(1)
      #current.model[j]<- as.numeric(rnd<ratio)*proposal.model[j]
      #lB.PMcurrent<- as.numeric(rnd<ratio)*lB.PMproposal+as.numeric(rnd>ratio)*lB.PMcurrent
      
      incl.probRB[j]<- incl.probRB[j]+proposal.model[j]*ratio+(1-proposal.model[j])*(1-ratio)	
      
    }
    setTxtProgressBar(pb, i)	  
	
	#if ((i/10-trunc(i/10))==0 & write.output){
	#	write.table(file=inclprob.file, matrix(round(incl.probRB/i, 3), nr=1), col.names=F, row.names=F, append=T)
	#	write.table(file=models.file, matrix(current.model, nr=1), col.names=F, row.names=F, append=T)
	#}	
  }
)

resumen.gibbs$inclprob<- incl.probRB/N
names(resumen.gibbs$inclprob)<-colnames(X)

