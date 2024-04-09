#Factor Bayes x Prior
B.PM<- function(model,y,X){
  #model is the binary vector specifying a model
  #SSE.null <- sum((y-mean(y))^2)
  SSE.model<- sum(.lm.fit(y=y, x=cbind(1, X[,model==1]))$residuals^2)
  BFi0<- .C("gBF", as.integer(dim(X)[1]), as.integer(sum(model)+1), 
            as.integer(1), as.double(SSE.model/SSE.null), 
            as.double(0))[5][[1]]		
  return(exp(log(BFi0)-lchoose(length(model), sum(model))))
  #return(runif(1,.99,1.01))
}

distrib<- function(set, cur.model, inclusion,y,X){
  B.current<- B.PM(cur.model,y,X)
  for (j in set){
    prop.m<- cur.model; prop.m[j]<- 1-cur.model[j]
    B.proposal<- B.PM(prop.m,y,X)
    ratio<- B.proposal/(B.proposal+B.current)
    
    if (runif(1)<ratio) {
      cur.model[j]<- prop.m[j]; B.current<- B.proposal
    }
    #rnd<- runif(1)
    #current.model[j]<- as.numeric(rnd<ratio)*proposal.model[j]
    #B.PMcurrent<- as.numeric(rnd<ratio)*B.PMproposal+as.numeric(rnd>ratio)*B.PMcurrent
    
    inclusion[j]<- inclusion[j]+prop.m[j]*ratio+(1-prop.m[j])*(1-ratio)
  }	
  rbind(inclusion[set], cur.model[set])
}


#Log (Factor Bayes x Prior)
lB.PM<- function(model, y, X.mat){
  #model is the binary vector specifying a model
  SSE.null <- sum((y-mean(y))^2)
  SSE.model<- sum(.lm.fit(y=y, x=cbind(1, X.mat[,model==1]))$residuals^2)
  BFi0<- .C("gBF", as.integer(dim(X.mat)[1]), as.integer(sum(model)+1), 
            as.integer(1), as.double(SSE.model/SSE.null), 
            as.double(0))[5][[1]]		
  return(log(BFi0)-lchoose(length(model), sum(model)))
}

#distrib que usa log (BF x prior)
distrib.l<- function(set, current.model, incl.probRB, y, X.mat){
  library(BayesVarSel)
  #Log (Factor Bayes x Prior)
  lB.PM<- function(model, y, X.mat){
    #model is the binary vector specifying a model
    SSE.null <- sum((y-mean(y))^2)
    SSE.model<- sum(.lm.fit(y=y, x=cbind(1, X.mat[,model==1]))$residuals^2)
    BFi0<- .C("gBF", as.integer(dim(X.mat)[1]), as.integer(sum(model)+1), 
              as.integer(1), as.double(SSE.model/SSE.null), 
              as.double(0))[5][[1]]		
    return(log(BFi0)-lchoose(length(model), sum(model)))
  }
  
  #Log (Factor Bayes x Prior)
  lB.PMcurrent<- lB.PM(current.model, y, X.mat)
  for (j in set){
    proposal.model<- current.model;  proposal.model[j]<- 1-current.model[j]
    lB.PMproposal<- lB.PM(proposal.model, y, X.mat)
    
    ratio<- 1/(1+exp(lB.PMcurrent-lB.PMproposal))
    
    if (runif(1)<ratio) {
      current.model[j]<- proposal.model[j]; lB.PMcurrent<- lB.PMproposal
    }
    
    #rnd<- runif(1)
    #current.model[j]<- as.numeric(rnd<ratio)*proposal.model[j]
    #lB.PMcurrent<- as.numeric(rnd<ratio)*lB.PMproposal+as.numeric(rnd>ratio)*lB.PMcurrent
    
    incl.probRB[j]<- incl.probRB[j]+proposal.model[j]*ratio+(1-proposal.model[j])*(1-ratio)
  }	
  rbind(incl.probRB[set], current.model[set])
}

##Calculus in the contamination experiment:


#Factor Bayes in contam. experiment x Prior
B.PM.cont<- function(modelX, nofesp, eta, p,X,y){
  if (sum(modelX==1)>0){
    SSE.model<- sum(.lm.fit(y=y, x=cbind(1, X[,modelX==1]))$residuals^2)
    BFi0<- .C("gBF", as.integer(dim(X)[1]), as.integer(sum(modelX)+1), 
              as.integer(1), as.double(SSE.model/SSE.null), 
              as.double(0))[5][[1]]		
    
    BB<- BFi0*eta^nofesp
    return(BB/choose(p, sum(modelX)+nofesp))
  }
  else return(eta^nofesp/choose(p, sum(modelX)+nofesp))
}


#LogFactor Bayes in contam. experiment x Prior
lB.PM.cont<- function(modelX, nofesp, eta, p,X,y){
  if (sum(modelX==1)>0){
    SSE.model<- sum(.lm.fit(y=y, x=cbind(1, X[,modelX==1]))$residuals^2)
    BFi0<- .C("gBF", as.integer(dim(X)[1]), as.integer(sum(modelX)+1), 
              as.integer(1), as.double(SSE.model/SSE.null), 
              as.double(0))[5][[1]]		
    lBB<- log(BFi0)+nofesp*log(eta)
    return(lBB-lchoose(p, sum(modelX)+nofesp))
  }
  else return(nofesp*log(eta)-lchoose(p, sum(modelX)+nofesp))
}




#####
#Functions
#####

normConst<- function(x, p, eta){
  kT<- length(x$inclprob)
  Sv<- x$postprobdim[2:(kT+1)]*x$C*2^kT #Suma factores bayes de cada dimension
  C<- 1
  for (d in 1:p){
    if (d < kT+1) {
      penal<- eta^(d:0)
      origSums<- c(1, Sv[1:d]) 
      combinat<- exp(lchoose(n=p-kT, k=d:0)-lchoose(n=p, k=d))
      C<- C+sum(penal*origSums*combinat)		
    }
    else
    {
      penal<- eta^(d:(d-kT))
      origSums<- c(1, Sv) 
      combinat<- exp(lchoose(n=p-kT, k=d:(d-kT))-lchoose(n=p, k=d))
      C<- C+sum(penal*origSums*combinat)						
    }
  }
  C/(p+1)
}

#another way (tested) implementing the notes in ContaminatedBF2.pdf
#This function implements the equation that appears in the draft. It is tested that gives the same results that normConst.
normConst2<- function(x, p, eta){
  kT<- length(x$inclprob)
  Sv<- x$postprobdim*x$C*2^kT
  C<- 0
  for (i in 0:kT){
    C<- C+Sv[i+1]*sum(eta^(0:(p-kT))*exp(lchoose(n=p-kT, k=0:(p-kT))-lchoose(n=p, i+(0:(p-kT)))))
  }
  C/(p+1)
}


sumBF.by.dim<- function(x){
  #x is a Bvs object with prior.models="Constant" and all models visited
  kT<- dim(x$modelsprob)[2]-1
  models<- as.matrix(x$modelsprob[,1:kT]=="*")*1
  dimens<- rowSums(models); 
  BFs<- x$modelsprob[,kT+1]*x$C*2^kT
  result<- matrix(0, ncol=kT, nrow=kT)
  for (i in 1:kT){
    estos<- which(models[,i]==1)
    aqui<- aggregate(BFs[estos]~dimens[estos], FUN=sum)
    result[i,]<- aqui[,2]		
  }
  rownames(result)<- x$variables
  colnames(result)<- paste("dim", 1:kT, sep="")
  return(result)
  
}


contam.incl.prob<- function(xi=1, C, p, eta, matrix.sumBF.by.dim){
  M<- matrix.sumBF.by.dim
  kT<- dim(M)[1]; result<- 0
  for (i in 1:kT){
    result<- result+M[xi,i]*sum(eta^(0:(p-i))/exp(lchoose(p, i:p)+log(C)-lchoose(p-kT, 0:(p-i))))
  }	
  result<- result/(p+1)
  return(result)	
}

#an alternative way (tested) implementing the notes in ContaminatedBF2.pdf
contam.incl.prob2<- function(xi=1, C, p, eta, matrix.sumBF.by.dim){
  M<- matrix.sumBF.by.dim
  kT<- dim(M)[1]; result<- 0
  for (i in 1:kT){
    result<- result+M[xi,i]*sum(eta^(0:(p-kT))/exp(lchoose(p, i+(0:(p-kT)))+log(C)-lchoose(p-kT, 0:(p-kT))))
  }	
  result<- result/(p+1)
  return(result)	
}