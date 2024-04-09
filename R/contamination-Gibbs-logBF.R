SSE.null<- sum((y-mean(y))^2)
current.model<- sample(c(rep(1, d.ini), rep(0, p-d.ini)))
lB.PMcurrent<- lB.PM.cont(current.model[realindex], nofesp=sum(current.model)-sum(current.model[realindex]), eta=eta, p=p,X=X,y=y)
proposal.model<- current.model

incl.probRB<- rep(0,p);

pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                     max = N, # Valor máximo de la barra de progreso
                     style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                     width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                     char = "=")   # Caracter usado para crear la barra

for (i in 1:N){
  for (j in 1:p){
    proposal.model<- current.model; proposal.model[j]<- 1-current.model[j]
    lB.PMproposal<- lB.PM.cont(proposal.model[realindex], nofesp=sum(proposal.model)-sum(proposal.model[realindex]), eta=eta, p=p,X=X,y=y)
    #ratio<- B.PMproposal/(B.PMproposal+B.PMcurrent)
    ratio<- 1/(1+exp(lB.PMcurrent-lB.PMproposal))
    if (runif(1)<ratio) {
      current.model[j]<- proposal.model[j]; lB.PMcurrent<- lB.PMproposal
    }
    
    incl.probRB[j]<- incl.probRB[j]+proposal.model[j]*ratio+(1-proposal.model[j])*(1-ratio)
  }
  setTxtProgressBar(pb, i)	  
}

incl.probRB<- incl.probRB/N
cat("\n")



