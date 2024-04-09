
resumen<- list()
cutoff<- NA
resumen$params<- c(N.ini, cutoff, nofproc, N)
names(resumen$params)<- c("N.ini", "cutoff", "nofproc", "N")

SSE.null <- sum((y-mean(y))^2)
incl.probRB<- list()
incl.probRB[[1]]<- rep(0, p)
incl.probRB[[2]]<- rep(0, p)
names(incl.probRB[[1]])<- colnames(X)
names(incl.probRB[[2]])<- colnames(X)

current.model<- list()
current.model[[1]]<- sample(c(rep(1, d.ini), rep(0, p-d.ini)))
current.model[[2]]<- sample(c(rep(1, d.ini), rep(0, p-d.ini)))

lB.PMcurrent<- list()
lB.PMcurrent[[1]]<- lB.PM(current.model[[1]], y=y, X.mat=X)
lB.PMcurrent[[2]]<- lB.PM(current.model[[2]], y=y, X.mat=X)

proposal.model<- list()
proposal.model[[1]]<- current.model[[1]]
proposal.model[[2]]<- current.model[[2]]

Xl<- list()
Xl[[1]]<- X
Xl[[2]]<- X


resumen$time.step1<- system.time(
  for (i in 1:N.ini){
    for (line in 1:2){
      #cat("Line:", line, "It:", i, "\n")
      for (j in 1:p){
        proposal.model[[line]]<- current.model[[line]]; proposal.model[[line]][j]<- 1-current.model[[line]][j]
        lB.PMproposal<- lB.PM(proposal.model[[line]], y=y, X.mat=Xl[[line]])
        #ratio<- B.PMproposal/(B.PMproposal+B.PMcurrent)
        ratio<- 1/(1+exp(lB.PMcurrent[[line]]-lB.PMproposal))
        
        if (runif(1)<ratio) {
          current.model[[line]][j]<- proposal.model[[line]][j]; lB.PMcurrent[[line]]<- lB.PMproposal
        }
        #rnd<- runif(1)
        #current.model[j]<- as.numeric(rnd<ratio)*proposal.model[j]
        #lB.PMcurrent<- as.numeric(rnd<ratio)*lB.PMproposal+as.numeric(rnd>ratio)*lB.PMcurrent
        
        incl.probRB[[line]][j]<- incl.probRB[[line]][j]+proposal.model[[line]][j]*ratio+(1-proposal.model[[line]][j])*(1-ratio)	
        
      }
      #temp.incl.probRB<- incl.probRB/i
      
    }
  }
)

#resumen$inclprob.step1 <- incl.probRB/N.ini
#Ordeno la matriz segun las inclusion probs y si tenian un uno
order.incl.probRB<- list(); esta<- list(); bloqueInicial<- list(); cl<- list()


for (line in 1:2){
  order.incl.probRB[[line]]<- order(incl.probRB[[line]]/N.ini+current.model[[line]], decreasing=T)
}

line<- 1; other.line<- 2
incl.probRB[[line]]<- incl.probRB[[line]][order.incl.probRB[[other.line]]]
Xl[[line]]<- Xl[[line]][, order.incl.probRB[[other.line]]]
current.model[[line]]<- current.model[[line]][order.incl.probRB[[other.line]]]
esta[[line]]<- esta.ini
bloqueInicial[[line]]<- 1:esta[[line]]


#line<- 2; other.line<- 2
#incl.probRB[[line]]<- incl.probRB[[line]][order.incl.probRB[[line]]]
#Xl[[line]]<- Xl[[line]][, order.incl.probRB[[other.line]]]
#current.model[[line]]<- current.model[[line]][order.incl.probRB[[other.line]]]
#esta[[line]]<- esta.ini
#bloqueInicial[[line]]<- 1:esta[[line]]

pb <- txtProgressBar(min = 0,      # Valor mínimo de la barra de progreso
                     max = N, # Valor máximo de la barra de progreso
                     style = 3,    # Estilo de la barra (también style = 1 y style = 2)
                     width = 50,   # Ancho de la barra. Por defecto: getOption("width")
                     char = "=")   # Caracter usado para crear la barra

cl <- makeCluster(nofproc)
#registerDoParallel(cl)

line<-1
resumen$time.step2<- system.time(
  for (i in 1:N){
     #Now, only for line 1:
      ####
      #Bloque I
      for (j in bloqueInicial[[line]]){
        lB.PMcurrent[[line]]<- lB.PM(current.model[[line]], y=y, X.mat=Xl[[line]])
        proposal.model[[line]]<- current.model[[line]]; proposal.model[[line]][j]<- 1-current.model[[line]][j]
        lB.PMproposal<- lB.PM(proposal.model[[line]], y=y, X.mat=Xl[[line]])
        ratio<- 1/(1+exp(lB.PMcurrent[[line]]-lB.PMproposal))		
        
        if (runif(1)<ratio) {
          current.model[[line]][j]<- proposal.model[[line]][j]; lB.PMcurrent[[line]]<- lB.PMproposal
        }
        
        #rnd<- runif(1)
        #current.model[j]<- as.numeric(rnd<ratio)*proposal.model[j]
        #lB.PMcurrent<- as.numeric(rnd<ratio)*lB.PMproposal+as.numeric(rnd>ratio)*lB.PMcurrent
        
        #
        incl.probRB[[line]][j]<- incl.probRB[[line]][j]+proposal.model[[line]][j]*ratio+(1-proposal.model[[line]][j])*(1-ratio)	
        #
      }
      
      #
      ###
      #Distribuye Bloques con ceros en el resto de procesadores
      #
      mb<- 1 
      esta[[line]]<- esta.ini
      while (mb < (p-nofproc) & distribute){
        #	
        bloques<- list()
        tam<- trunc((p-esta[[line]])/nofproc)
        cortes<- esta[[line]]+1+(1:nofproc)*tam; cortes[nofproc]<- p
        #
        bloques[[1]]<- (esta[[line]]+1):cortes[1]
        for (h in 2:nofproc){
          bloques[[h]]<- (cortes[h-1]+1):cortes[h] 
        }
        
        #Parallel with for each and dopar:
        #result<- foreach (h=1:nofproc, .combine='cbind',.packages=c("BayesVarSel")) %dopar% {
        #set.seed(h+as.numeric(Sys.time()))
        #distrib.l(bloques[[h]], current.model[[line]], incl.probRB[[line]], y=y, X.mat=Xl[[line]])
	      #}
        
        ##Parallel with parLapply:
        result<- parLapply(cl,X=bloques, fun=distrib.l, current.model[[line]], incl.probRB[[line]],y,X.mat=Xl[[line]])
        result<- do.call(cbind, result)
        
        #result<- mclapply(X=bloques, FUN=distrib.l, mc.cores=nofproc, mc.preschedule=FALSE, current.model[[line]], incl.probRB[[line]],y, X.mat=Xl[[line]])
        #result<- do.call(cbind, result)
        
        
        
        #Sin paralelizar, haciendolo con un for parece que si funciona:
        #result<- NULL
        #for(h in 1:nofproc){
        #  result<- cbind(result,distrib.l(bloques[[h]], current.model[[line]], incl.probRB[[line]], y=y, X=Xl[[line]]))
        #} 
        
        
        #
        #combine results:
        #el bloque hasta el que me sirve
        #cuidado, ¿que pasa si el 1 aparece al final, con mas procesadores para repartir que modelos?
        thisok<- min(c(trunc(min(which(result[2,]==1), p)*nofproc/p)+1, nofproc))
        #cat("real:", which(result[2,]==1), "\n")
        #cat("estimated-mean", 1/mean(incl.probRB[(esta+1):p]/(N.ini+i)), "\n")
        #cat("estimated-max", 1/max(incl.probRB[(esta+1):p]/(N.ini+i)), "\n")
        
        mb<- max(bloques[[thisok]]); 
        #cat("mb", mb, "\n")
        #
        if (mb<p){
          incl.probRB[[line]]<- c(incl.probRB[[line]][1:esta[[line]]], result[1, 1:(mb-esta[[line]])], incl.probRB[[line]][(mb+1): p])
          current.model[[line]]<- c(current.model[[line]][1:esta[[line]]], result[2, 1:(mb-esta[[line]])], current.model[[line]][(mb+1):p])
        }else{
          incl.probRB[[line]]<- c(incl.probRB[[line]][1:esta[[line]]], result[1, 1:(mb-esta[[line]])])
          current.model[[line]]<- c(current.model[[line]][1:esta[[line]]], result[2, 1:(mb-esta[[line]])])
        }
        esta[[line]]<- mb
        #
      }
      
      #Bloque final (si lo hay)
      if (mb < p & distribute){
        #result<- distrib(current.model, incl.probRB, bloques[[nofproc]])
        #incl.probRB<- c(incl.probRB[1:cortes[nofproc-1]], result[1,])
        #current.model<- c(current.model[1:cortes[nofproc-1]], result[2,])
        result<- distrib.l((mb+1):p, current.model[[line]], incl.probRB[[line]], y, Xl[[line]])
        incl.probRB[[line]]<- c(incl.probRB[[line]][1:mb], result[1,])
        current.model[[line]]<- c(current.model[[line]][1:mb], result[2,])
      }
      #
      # 
      setTxtProgressBar(pb, i)
    }#i-th iteration completed for both lines
    #
    
    
  
    #cat("esta", esta[[1]], esta[[2]], "\n")
    #cat("dimension", sum(current.model[[1]]), sum(current.model[[2]]),  "\n")
    #
    # 
    )

cat("\n")
stopCluster(cl)
incl.probRB[[1]]<- incl.probRB[[1]]/(N+N.ini)
#incl.probRB[[2]]<- incl.probRB[[2]]/(N+N.ini)

resumen$inclprob.step2 <- incl.probRB[[1]][order(order.incl.probRB[[2]])]
names(resumen$inclprob.step2)<- colnames(X)


