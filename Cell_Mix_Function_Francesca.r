
plotNet<-function(Sigma1,Sigma2){
  library(MASS)
  library(igraph)
  rho1<-ginv(Sigma1)
  rho2<-ginv(Sigma2)
  adj.tr<-rho1; adj.tr[adj.tr<1e-6]<-0; adj.tr[adj.tr>=1e-6]<-1
  adj.tr2<-rho2; adj.tr2[adj.tr2<1e-6]<-0; adj.tr2[adj.tr2>=1e-6]<-1
  graph1<-graph.adjacency(as.matrix(adj.tr), mode="undirected")
  E(graph1)$col<-"red"
  graph2<-graph.adjacency(as.matrix(adj.tr2), mode="undirected")
  E(graph2)$col<-"black"
  new.j<-graph.union(graph1 %u% graph2)
  
  t.l<-adj.tr[lower.tri(adj.tr)]; new.l<-adj.tr2[lower.tri(adj.tr2)];
  t.l[new.l==0 & t.l ==1]=2;
  t.l[new.l==1 & t.l ==1]=3;
  t.l[new.l==1 & t.l ==0]=1;
  
  t.l<-t.l[t.l!=0]; t.l<-t.l[seq(length(t.l),1)]
  
  E(new.j)$color[t.l==3]<-"green"
  E(new.j)$color[t.l==2]<-"blue"
  E(new.j)$color[t.l==1]<-"red"
  
  v.size=1
  set.seed(1)
  ll <- layout.reingold.tilford(new.j)

  # -- remove modules with less than 4 vertices
  C<-clusters(new.j, mode=c("weak"))
  bad.vs<-V(new.j)[C$csize[C$membership]<4]     # remove components with less than 4 vertices
  new.j<-delete.vertices(new.j, bad.vs) # exclude them from the graph
  
  plot(new.j, vertex.size=v.size, vertex.frame.color="white",layout=layout.fruchterman.reingold, vertex.label=NA)
  legend("topleft",c("Shared","Tumor Specific","Normal Specific"),col=c("green","red","blue"),lwd=c(.5,.5),cex=.7)
}

simulateDataN_normal<-function(dn,cluster,sn,parameterseed,dataseed, P.v=NULL,sameModuleN=0,
                        betaSize=Inf,sigmaFactor=1,sigmaRatio=1,betaModel=TRUE,type=1,n.normal){
  ## number of clusters
  ## dn: dimension (number of genes)
  ## sn: sample size
  ## c: a scale; unknown parameter to scale the percentage estimation
  ## n.normal: # of samples are sampled from the normal component
  
  library(MASS)
  sizeL<-rep(floor(dn/cluster),cluster)
  
  ####simulate parameters
  set.seed(parameterseed)
  Uy<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  Uz<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  
  Sigmay<-matrix(0,dn,dn);
  Sigmaz<-Sigmay;
  
  pos<-1;
  for(i in 1:length(sizeL)) {
    while(1)
    {
      if (type==1)  x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      if(length(x)>1) break()
    }       
    Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    pos<-pos+sizeL[i];
  }
  
  pos<-1;
  for(i in 1:length(sizeL))
  {
    while(1)
    {
      if (type==1) x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      
      if(length(x)>1) break()
    }       
    Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    if(i<=sameModuleN) { # --- normal = tumor
      Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)];
      Uz[pos:(pos+sizeL[i]-1)]<-Uy[pos:(pos+sizeL[i]-1)];
    }
    pos<-pos+sizeL[i];
  }
  
  Sigmay<-Sigmay*sigmaFactor*sqrt(sigmaRatio)
  Sigmaz<-Sigmaz*sigmaFactor/sqrt(sigmaRatio)
  
  if(is.null(P.v))
  {
    #realP<-runif(sn)*0.6+0.2
    meanP<-0.5;
    varP<-0.04;
    sizeP<-meanP*(1-meanP)/varP-1
    realP<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
    
  } else { realP=P.v}
  
  realP[sample(sn,n.normal)]<-0 # -- n.normal samples have purity equal to zero
  
  ## -- Simulate observed data
  set.seed(dataseed)
  Y<-mvrnorm(sn,Uy,Sigmay)  # -- sample Y and Z from a normal distribution
  Z<-mvrnorm(sn,Uz,Sigmaz)
  X<-sweep(Y,1,realP,"*")+sweep(Z,1,1-realP,"*");
  
  if(betaSize==Inf) P<-realP
  if(betaSize!=Inf & betaModel==FALSE){d<-rnorm(length(realP),log(realP/(1-realP)),sqrt(betaSize)); P<-exp(d)/(1+exp(d))} # -- noisy P
  if(betaSize!=Inf & betaModel==TRUE){P<-rbeta(length(realP),betaSize*realP,betaSize*(1-realP))} # -- noisy P
  
  P[P>.99]<-.99; P[P<.01]<-.01
  
  print(cor(P,realP))
  
  c<-5;
  A<-P/c;
  
  # -- P is the noisy purity
  # -- realP the real purity
  return(list(X=X,A=A,Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz,c=c,P=P,Y=Y,Z=Z,realP=realP));
}

simulateDataN_normal_log<-function(dn,cluster,sn,parameterseed,dataseed, P.v=NULL,sameModuleN=0,
                               betaSize=Inf,sigmaFactor=1,sigmaRatio=1,betaModel=TRUE,type=1,n.normal){
  ## sample from mixture of log normal
  ## number of clusters
  ## dn: dimension (number of genes)
  ## sn: sample size
  ## c: a scale; unknown parameter to scale the percentage estimation
  ## n.normal: # of samples are sampled from the normal component
  library(MASS)
  sizeL<-rep(floor(dn/cluster),cluster)
  
  ####simulate parameters
  set.seed(parameterseed)
  Uy<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  Uz<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
  
  Sigmay<-matrix(0,dn,dn);
  Sigmaz<-Sigmay;
  
  pos<-1;
  for(i in 1:length(sizeL)) {
    while(1)
    {
      if (type==1)  x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      if(length(x)>1) break()
    }       
    Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    pos<-pos+sizeL[i];
  }
  
  pos<-1;
  for(i in 1:length(sizeL))
  {
    while(1)
    {
      if (type==1) x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      
      if(length(x)>1) break()
    }       
    Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
    if(i<=sameModuleN) { # --- normal = tumor
      Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)];
      Uz[pos:(pos+sizeL[i]-1)]<-Uy[pos:(pos+sizeL[i]-1)];
    }
    pos<-pos+sizeL[i];
  }
  
  Sigmay<-Sigmay*sigmaFactor*sqrt(sigmaRatio)
  Sigmaz<-Sigmaz*sigmaFactor/sqrt(sigmaRatio)
  
  if(is.null(P.v))
  {
    #realP<-runif(sn)*0.6+0.2
    meanP<-0.5;
    varP<-0.04;
    sizeP<-meanP*(1-meanP)/varP-1
    realP<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
    
  } else { realP=P.v}
  
  realP[sample(sn,n.normal)]<-0 # -- n.normal samples have purity equal to zero
  
  ## -- Simulate observed data
  set.seed(dataseed)
  Y<-2^{mvrnorm(sn,Uy,Sigmay)}  # Y and Z are log2 Normal distribution
  Z<-2^{mvrnorm(sn,Uz,Sigmaz)}
  X<-sweep(Y,1,realP,"*")+sweep(Z,1,1-realP,"*");
  
  if(betaSize==Inf) P<-realP
  if(betaSize!=Inf & betaModel==FALSE){d<-rnorm(length(realP),log(realP/(1-realP)),sqrt(betaSize)); P<-exp(d)/(1+exp(d))} # -- noisy P
  if(betaSize!=Inf & betaModel==TRUE){P<-rbeta(length(realP),betaSize*realP,betaSize*(1-realP))} # -- noisy P
  
  P[P>.99]<-.99; P[P<.01]<-.01
  
  print(cor(P,realP))
  
  c<-5;
  A<-P/c;
  
  # -- P is the noisy purity
  # -- realP the real purity
  return(list(X=X,A=A,Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz,c=c,P=P,Y=Y,Z=Z,realP=realP));
}

simulateDataN<-function(dn,cluster,sn,parameterseed,dataseed, P.v=NULL,sameModuleN=0,
                        betaSize=Inf,sigmaFactor=1,sigmaRatio=1,betaModel=TRUE,type=1){
## number of clusters
## dn: dimension (number of genes)
## sn: sample size
## c: a scale; unknown parameter to scale the percentage estimation
library(MASS)
sizeL<-rep(floor(dn/cluster),cluster)

####simulate parameters
set.seed(parameterseed)
Uy<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal
Uz<-rnorm(dn,0,1)  # -- sample mean parameter from standard normal

Sigmay<-matrix(0,dn,dn);
Sigmaz<-Sigmay;

pos<-1;
for(i in 1:length(sizeL)) {
         while(1)
           {
      if (type==1)  x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
      if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      if(length(x)>1) break()
      }       
      Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
      pos<-pos+sizeL[i];
      }

pos<-1;
for(i in 1:length(sizeL))
       {
          while(1)
           {
            if (type==1) x<-try(PowerLaw.net(sizeL[i], umin=0.5, umax=1, power=0.5, flip=TRUE, factor=1.5, thre=0.1))
            if (type==2) x<-try(Star.net(sizeL[i], hub.no=2,hub.size=16,umin=0.5,umax=1))
      
      if(length(x)>1) break()
      }       
Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-x;
if(i<=sameModuleN) { # --- normal = tumor
  Sigmaz[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)]<-Sigmay[pos:(pos+sizeL[i]-1),pos:(pos+sizeL[i]-1)];
  Uz[pos:(pos+sizeL[i]-1)]<-Uy[pos:(pos+sizeL[i]-1)];
}
pos<-pos+sizeL[i];
      }

Sigmay<-Sigmay*sigmaFactor*sqrt(sigmaRatio)
Sigmaz<-Sigmaz*sigmaFactor/sqrt(sigmaRatio)

if(is.null(P.v))
{
  #realP<-runif(sn)*0.6+0.2
  meanP<-0.5;
  varP<-0.04;
  sizeP<-meanP*(1-meanP)/varP-1
  realP<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
   
} else { realP=P.v}


## -- Simulate observed data
set.seed(dataseed)
Y<-mvrnorm(sn,Uy,Sigmay)
Z<-mvrnorm(sn,Uz,Sigmaz)
X<-sweep(Y,1,realP,"*")+sweep(Z,1,1-realP,"*");

if(betaSize==Inf) P<-realP
if(betaSize!=Inf & betaModel==FALSE){d<-rnorm(length(realP),log(realP/(1-realP)),sqrt(betaSize)); P<-exp(d)/(1+exp(d))} # -- noisy P
if(betaSize!=Inf & betaModel==TRUE){P<-rbeta(length(realP),betaSize*realP,betaSize*(1-realP))} # -- noisy P

P[P>.99]<-.99; P[P<.01]<-.01
print(cor(P,realP))

c<-5;
A<-P/c;

# -- P is the noisy purity
# -- realP the real purity
return(list(X=X,A=A,Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz,c=c,P=P,Y=Y,Z=Z,realP=realP));
}


simulateDataN_diag<-function(dn,cluster,sn,parameterseed,dataseed, P.v=NULL,sameModuleN=0,
                        betaSize=Inf,sigmaFactor=1,sigmaRatio=1,betaModel=TRUE,type=1,Uy1,Uy2,Uz1,Uz2){
  ## number of clusters
  ## dn: dimension (number of genes)
  ## sn: sample size
  ## c: a scale; unknown parameter to scale the percentage estimation
  ## Uy1: mean of tumor cell under data 1
  ## Uz1: mean of normal/stromal cell under data1
  ## Uy2: mean of tumor cell under data 2
  ## Uz2: mean of normal/stromal cell under data2
  
  library(MASS)
  sizeL<-rep(floor(dn/cluster),cluster)
  
  ####simulate parameters
  set.seed(parameterseed)

  Sigmay<-matrix(0,dn,dn);
  Sigmaz<-Sigmay;
  diag(Sigmaz)<-diag(Sigmay)<-rep(1,dn)

  #realP<-runif(sn)*0.6+0.2
    meanP<-0.5;
    varP<-0.1;
    sizeP<-meanP*(1-meanP)/varP-1
    realP1<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
    realP2<-rbeta(sn,meanP*sizeP,(1-meanP)*sizeP)
    
  
  ## -- Simulate observed data
  set.seed(dataseed)
  Y1<-mvrnorm(sn,Uy1,Sigmay)
  Z1<-mvrnorm(sn,Uz1,Sigmaz)
  X1<-sweep(Y1,1,realP1,"*")+sweep(Z1,1,1-realP1,"*");
  
  Y2<-mvrnorm(sn,Uy2,Sigmay)
  Z2<-mvrnorm(sn,Uz2,Sigmaz)
  X2<-sweep(Y2,1,realP2,"*")+sweep(Z2,1,1-realP2,"*");
  
  P1<-rbeta(length(realP1),betaSize*realP1,betaSize*(1-realP1)) # -- noisy P
  P2<-rbeta(length(realP2),betaSize*realP2,betaSize*(1-realP2)) # -- noisy P

  
  # -- P is the noisy purity
  # -- realP the real purity
  return(list(X1=X1,X2=X2,Uy1=Uy1,Uz1=Uz1,Uy2=Uy2,Uz2=Uz2,Sigmay=Sigmay,Sigmaz=Sigmaz,
              P1=P1,P2=P2,Y1=Y1,Z1=Z1,Y2=Y2,Z2=Z2,realP1=realP1,realP2=realP2));
}

############################################################################################################################

Estep.noloop<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
  {
### Uy: vector of length p
### Uz: vector of length p
### Sigmay: covariance matrix of dimension pxp
### Sigmaz: covariance matrix of dimension pxp
### X: data matrix of dimnsion nxp
### A: abundance based percentage estimate; a vector of length n

   #Uy=r$Uy;Uz=r$Uz; Sigmay=r$Sigmay; Sigmaz=r$Sigmaz; c.scale=1; X=exprM; A=Epurity

    P<-c.scale*A  ### percentage estimate
    P.m=diag(P)
    p=ncol(X)
    n=nrow(X)


    #begin=proc.time()
 
    Ey<-matrix(NA,n, p)
    Ez<-Ey
    Esigmay<-matrix(0,p,p)
    Esigmaz<-Esigmay

    temp=eigen(Sigmay, symmetric=TRUE)
    M=temp$values
    U=t(temp$vectors)
    ### temp1=t(U)%*%diag(M)%*%U
    ### max(abs(temp1-Sigmay))
    Sy.half=t(U)%*%diag(M^0.5)
    Sy.half.inv=diag(1/M^0.5)%*%U
    ### max(abs(Sy.half%*%t(Sy.half)-Sigmay))

    W=Sy.half.inv%*%Sigmaz%*%t(Sy.half.inv)
    temp1=eigen(W, symmetric=TRUE)
    V=t(temp1$vectors)
    D=temp1$values
    ### temp2=t(V)%*%diag(D)%*%V
    ### max(abs(temp2-W))

    A=Sy.half%*%t(V)
    B=V%*%Sy.half.inv ## p by p 
    
    BX=B%*%t(X)%*%P.m ##p by n
    BUy=(B%*%Uy)%*%matrix(P^2, 1, n) ## p by n 
    BUz=(B%*%Uz)%*%matrix(P*(1-P),1,n) ## p by n
    tempB=BX-BUy-BUz ## p by n
    
    tempP2=matrix(P^2, p,n, byrow=TRUE) ##p by n
    tempC=tempP2+matrix(D, p,1)%*%matrix((1-P)^2, 1,n) ## p by n 
    tempC=1/tempC
    
    tempBC=tempB*tempC ##p by n 
    tempABC=A%*%tempBC ## p by n
    
    Ey=t(matrix(Uy, p,n)+tempABC) ##  n by p
    
    Cbar= diag(apply(tempP2*tempC, 1, mean)) ##p by p diagonal 
    Esigmay=Sigmay-A%*%Cbar%*%t(A)
    
    P.2m=diag(1/(1-P))
    Ez=P.2m%*%(X-P.m%*%Ey)
    tempP4=matrix(P^4/(1-P)^2, p,n, byrow=TRUE) ##p by n
    CbarZ= diag(apply(tempP4*tempC, 1, mean)) ##p by p 
    Esigmaz=mean(P^2/(1-P)^2)*Sigmay-A%*%CbarZ%*%t(A) ## p by p 
      
     
    Esigmay<-(Esigmay+t(Esigmay))/2  ### garantee its symmetry
    Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry

    #end=proc.time()  
    list(Ey=Ey,Ez=Ez,Esigmay=Esigmay,Esigmaz=Esigmaz)
  }



###
Estep<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
  {
### Uy: vector of length p
### Uz: vector of length p
### Sigmay: covariance matrix of dimension pxp
### Sigmaz: covariance matrix of dimension pxp
### X: data matrix of dimnsion nxp
### A: abundance based percentage estimate; a vector of length n

    P<-c.scale*A  ### percentage estimate

    p=ncol(X)
    n=nrow(X)

    Ey<-matrix(NA,n, p)
    Ez<-Ey
    Esigmay<-matrix(0,p,p)
    Esigmaz<-Esigmay

    for(i in 1:n)
      {
    Ux<-P[i]*Uy+(1-P[i])*Uz;
    Sigmax<-P[i]^2*Sigmay+(1-P[i])^2*Sigmaz
    Sigmaxinv<-ginv(Sigmax)
    Sigmaxy<-P[i]*Sigmay
    Ey[i,]<-Uy+Sigmaxy%*%Sigmaxinv%*%(X[i,]-Ux)
    Esigmayi<-Sigmay-Sigmaxy%*%Sigmaxinv%*%t(Sigmaxy)
    Esigmay<-Esigmay+Esigmayi
    Ez[i,]<-(X[i,]-P[i]*Ey[i,])/(1-P[i]);
    Esigmaz<-Esigmaz+P[i]^2/(1-P[i])^2*Esigmayi;
  }
    Esigmay<-Esigmay/nrow(X);
    Esigmaz<-Esigmaz/nrow(X);
    
    Esigmay<-(Esigmay+t(Esigmay))/2  ### garantee its symmetry
    Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry
    
    list(Ey=Ey,Ez=Ez,Esigmay=Esigmay,Esigmaz=Esigmaz)
  }

######################################

Estep.noloop<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
  {
### Uy: vector of length p
### Uz: vector of length p
### Sigmay: covariance matrix of dimension pxp
### Sigmaz: covariance matrix of dimension pxp
### X: data matrix of dimnsion nxp
### A: abundance based percentage estimate; a vector of length n

   #Uy=r$Uy;Uz=r$Uz; Sigmay=r$Sigmay; Sigmaz=r$Sigmaz; c.scale=1; X=exprM; A=Epurity

    P<-c.scale*A  ### percentage estimate
    P.m=diag(P)
    p=ncol(X)
    n=nrow(X)

    Ey<-matrix(NA,n, p)
    Ez<-Ey
    Esigmay<-matrix(0,p,p)
    Esigmaz<-Esigmay

    temp=eigen(Sigmay, symmetric=TRUE)
    M=temp$values
    U=t(temp$vectors)
    ### temp1=t(U)%*%diag(M)%*%U
    ### max(abs(temp1-Sigmay))
    Sy.half=t(U)%*%diag(M^0.5)
    Sy.half.inv=diag(1/M^0.5)%*%U
    ### max(abs(Sy.half%*%t(Sy.half)-Sigmay))

    W=Sy.half.inv%*%Sigmaz%*%t(Sy.half.inv)
    temp1=eigen(W, symmetric=TRUE)
    V=t(temp1$vectors)
    D=temp1$values
    ### temp2=t(V)%*%diag(D)%*%V
    ### max(abs(temp2-W))

    A=Sy.half%*%t(V)
    B=V%*%Sy.half.inv ## p by p 
    
    BX=B%*%t(X)%*%P.m ##p by n
    BUy=(B%*%Uy)%*%matrix(P^2, 1, n) ## p by n 
    BUz=(B%*%Uz)%*%matrix(P*(1-P),1,n) ## p by n
    tempB=BX-BUy-BUz ## p by n
    
    tempP2=matrix(P^2, p,n, byrow=TRUE) ##p by n
    tempC=tempP2+matrix(D, p,1)%*%matrix((1-P)^2, 1,n) ## p by n 
    tempC=1/tempC
    
    tempBC=tempB*tempC ##p by n 
    tempABC=A%*%tempBC ## p by n
    
    Ey=t(matrix(Uy, p,n)+tempABC) ##  n by p
    
    Cbar= diag(apply(tempP2*tempC, 1, mean)) ##p by p diagonal 
    Esigmay=Sigmay-A%*%Cbar%*%t(A)
    
    P.2m=diag(1/(1-P))
    Ez=P.2m%*%(X-P.m%*%Ey)
    
    tempP4=matrix(P^4/(1-P)^2, p,n, byrow=TRUE) ##p by n
    CbarZ= diag(apply(tempP4*tempC, 1, mean)) ##p by p 
    Esigmaz=mean(P^2/(1-P)^2)*Sigmay-A%*%CbarZ%*%t(A) ## p by p 
      
     
    Esigmay<-(Esigmay+t(Esigmay))/2  ### garantee its symmetry
    Esigmaz<-(Esigmaz+t(Esigmaz))/2  ### garantee its symmetry

    #end=proc.time()  
    list(Ey=Ey,Ez=Ez,Esigmay=Esigmay,Esigmaz=Esigmaz)
  }


#############################################################################################################################
Mstep.v2<-function(Ey,Ez,Esigmay,Esigmaz,rhoy=rhoy, rhoz=rhoz, Uy,Uz)
  {
### Ey: matrix of dimension nxp
### Ez: matrix of dimension nxp
### Esigmay: Conditional Expectation of the covariance matrix of dimension pxp
### Esigmaz: Conditional Expectation of the covariance matrix of dimension pxp

    p=ncol(Ey);
    n=nrow(Ey);

    if(missing(Uy))Uy=apply(Ey,2,mean);
    if(missing(Uz))Uz=apply(Ez,2,mean);

    #temp<-matrix(0,p, p)
    #for(i in 1:n)
    #  {
    #    temp<-temp+(Ey[i,]-Uy)%*%t(Ey[i,]-Uy)
    #  }
    # Sy<-temp/n+Esigmay;
    temp=matrix(Uy, nrow=n, ncol=p, byrow=T)
    Sy=t(Ey-temp)%*%(Ey-temp)/n+Esigmay
    Sigmay<-glasso(Sy,rhoy, penalize.diagonal=FALSE)$w;

    #temp<-matrix(0,p, p)
    #for(i in 1:n)
    #  {
    #    temp<-temp+(Ez[i,]-Uz)%*%t(Ez[i,]-Uz)
    #  }
    #Sz<-temp/n+Esigmaz;
    temp=matrix(Uz, nrow=n, ncol=p, byrow=T)
    Sz=t(Ez-temp)%*%(Ez-temp)/n+Esigmaz
    Sigmaz<-glasso(Sz,rhoz, penalize.diagonal=FALSE)$w;

    list(Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz)
  }


######
Mstep<-function(Ey,Ez,Esigmay,Esigmaz,rho=rho,Uy,Uz)
  {
### Ey: matrix of dimension nxp
### Ez: matrix of dimension nxp
### Esigmay: Conditional Expectation of the covariance matrix of dimension pxp
### Esigmaz: Conditional Expectation of the covariance matrix of dimension pxp

    p=ncol(Ey);
    n=nrow(Ey);

    if(missing(Uy))Uy=apply(Ey,2,mean);
    if(missing(Uz))Uz=apply(Ez,2,mean);

    temp<-matrix(0,p, p)
    for(i in 1:n)
      {
        temp<-temp+(Ey[i,]-Uy)%*%t(Ey[i,]-Uy)
      }
    Sy<-temp/n+Esigmay;
    Sigmay<-glasso(Sy,rho)$w;

    temp<-matrix(0,p, p)
    for(i in 1:n)
      {
        temp<-temp+(Ez[i,]-Uz)%*%t(Ez[i,]-Uz)
      }
    Sz<-temp/n+Esigmaz;
    Sigmaz<-glasso(Sz,rho)$w;
     #Sigmaz<-Sz

    list(Uy=Uy,Uz=Uz,Sigmay=Sigmay,Sigmaz=Sigmaz)
  }


calculateLL<-function(Uy,Uz,Sigmay,Sigmaz,c.scale,X,A)
  {
     P<-c.scale*A
     sum(sapply(1:nrow(X),function(i)
       {
     Ux<-P[i]*Uy+(1-P[i])*Uz;
    Sigmax<-P[i]^2*Sigmay+(1-P[i])^2*Sigmaz
    dmvnorm(X[i,], Ux, Sigmax, log = TRUE)
   }))
 }

#######################
Estep.1d<-function(uy,uz,vy,vz, c,X,A)
  {
    P<-c*A
    Ux<-P*uy+(1-P)*uz;
    Vx<-P^2*vy+(1-P)^2*vz
    Sxy<-P*vy
    Ey<-uy+(X-Ux)*Sxy/Vx
    Vy<-vy-Sxy^2/Vx
    
    Ey2<-Vy+Ey^2

    Sxz<-(1-P)*vz
    Ez<-uz+(X-Ux)*Sxz/Vx
    Vz<-vz-Sxz^2/Vx
    Ez2<-Vz+Ez^2
    #Ez<-(X-Ey*P)/(1-P);
    #Ez2<-(X^2-2*P*X*Ey+P^2*Ey2)/(1-P)^2;

    list(Ey=Ey,Ey2=Ey2, Ez=Ez, Ez2=Ez2)
  }

Mstep.1d<-function(Ey,Ey2,Ez,Ez2)
  {
    uy<-mean(Ey)
    uz<-mean(Ez)
    vy<-mean(Ey2)-uy^2
    vz<-mean(Ez2)-uz^2
    list(uy=uy,uz=uz,vy=vy,vz=vz)
  }

####
calculateLL.1d<-function(uy,uz,vy,vz,c,X,A)
  {
     P<-c*A
    Ux<-P*uy+(1-P)*uz;
    Vx<-P^2*vy+(1-P)^2*vz
    sum(log( sapply(1:length(X),function(i)dnorm(X[i],Ux[i],sqrt(Vx[i])))))
  }

#####
calMSE<-function(est, tru)
{
   result=rep(0, 4)
   result[1]=mean((est$Uy-tru$Uy)^2)
   result[2]=mean((est$Uz-tru$Uz)^2)
   result[3]=mean((est$Sigmay-tru$Sigmay)^2)
   result[4]=mean((est$Sigmaz-tru$Sigmaz)^2)
   return(result)
}

#######################################################################################################################################################################################
####################NW Decoding


deNet.purity<-function(exprM,purity,estimateP=TRUE,bootstrap=0,betaModel=TRUE)
{
  
  t<-proc.time()
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  registerDoMC(10)
  library("foreach")
  options(warn=-1)
  source("/sc/orga/projects/PGDAC/francesca/TSNet/Codes/Cell_Mix_Function_Francesca.r")
  convergecutoff=0.001 ### convergence criterion for EM
  if(bootstrap>0){
    set.seed(bootstrap);
    rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
    exprM<-exprM[,rindex];
    purity<-purity[rindex];
    geneNames<-rownames(exprM);
    exprM<-exprM[apply(exprM,1,sd)!=0,];
    curgeneNames<-rownames(exprM);
    ###standardization
  }
  exprM<-t(scale((exprM)))
  geneNames<-rownames(exprM);
 
  curgeneNames<-rownames(exprM);
  
  if(estimateP)
  {
    for(j in 1:1000)  # --- over EM iterations 
      # --- (CHECK: why do we need this loop since we have (A Loop) and (B Loop))
    {
      
      if(j==1){
        cur.P<-purity;
        cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
      }else{
         print(paste("Epurity",j,":",cor(purity,cur.P)))
      }
      
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
      
      paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
  x.v=exprM[pi,] 
  for(i in 1:10000) # -- (A Loop) over EM iterations
  {
    if(i==1)
    {
      rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
    }
    r<-rnew;
    temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
    rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
    if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
    {#print(paste("converge in step",i));
      break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}
      
      paraM<-do.call("rbind",paraM)
      colnames(paraM)<-c("uy","uz","vy","vz")
      print("EM1d")
      before.P<-cur.P;
      
      ## -- Estimate purity (P) and delta parameter in beta distribution
      for(i in 1:10000) {                          
        
        if (betaModel==FALSE){ 
          
          new.P<-sapply(1:ncol(exprM),function(ni){
          f<-function(x){
            V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
            U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
            0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dnorm(log(purity[ni]/(1-purity[ni])),log(x/(1-x)),sqrt(cur.b))) # purity
          }
          optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
          new.b<-sum((log(purity/(1-purity))-log(new.P/(1-new.P)))^2)/length(purity)
        
          } else {
          
            new.P<-sapply(1:ncol(exprM),function(ni){
              f<-function(x) {
              V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
              U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
              0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
            }
            optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
            f<-function(x) -sum(log(dbeta(purity,new.P*x,(1-new.P)*x)))
            new.b<-optim(cur.b,f,lower=0)$par        
        }
        
        if(all(abs(c(new.P-cur.P))<0.01))
        {#print(paste("converge in step",i));
          break()}
        cur.b<-new.b;
        cur.P<-new.P
        print(cur.b)
        print(mean(abs(cur.P-purity)))
      }
      print(paste("cur.b",cur.b))
      
      
      if(all((before.P-cur.P)<0.001))break();
      print(mean(abs(cur.P-before.P)))
    }
    
    if(j==10000){ stop("P not converge"); convergence="P not converge"}else{
      print("Epurity converged"); convergence="P converged"}
    
    Epurity<-cur.P
  }else{Epurity<-purity};
  
  time<-(proc.time()-t)
  return(list(Epurity,paraM,convergence,time))

}

deNet.only.mean<-function(exprM,purity,estimateP=TRUE,bootstrap=0,betaModel=TRUE)
{
  
  t<-proc.time()
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  registerDoMC(10)
  library("foreach")
  options(warn=-1)
  source("/sc/orga/projects/PGDAC/francesca/TSNet/Codes/Cell_Mix_Function_Francesca.r")
  convergecutoff=0.001 ### convergence criterion for EM
  if(bootstrap>0){
    set.seed(bootstrap);
    rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
    exprM<-exprM[,rindex];
    purity<-purity[rindex];
    geneNames<-rownames(exprM);
    exprM<-exprM[apply(exprM,1,sd)!=0,];
    curgeneNames<-rownames(exprM);
    ###standardization
  }
  exprM<-t(scale((exprM)))
  geneNames<-rownames(exprM);
  
  curgeneNames<-rownames(exprM);
  
      cur.P<-purity;
      cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
      
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
      
  
  paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
  {
    x.v=exprM[pi,] 
    for(i in 1:10000) # -- (A Loop) over EM iterations
    {
      if(i==1)
      {
        rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
      }
      r<-rnew;
      temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
      rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
      if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
      {#print(paste("converge in step",i));
        break()}
    }
    c(r$uy,r$uz,r$vy,r$vz)
  }
      
      paraM<-do.call("rbind",paraM)
      colnames(paraM)<-c("uy","uz","vy","vz")
      Epurity<-purity;
  
  return(list(Epurity,paraM))
  
}


deNet.purity.normal<-function(exprM,exprM.normal,purity,estimateP=TRUE,
                              bootstrap=0,betaModel=TRUE)
{
  # --- add normal tissue samples
  t<-proc.time()
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  registerDoMC(10)
  library("foreach")
  options(warn=-1)
  source("Codes/Cell_Mix_Function_Francesca.r")
  convergecutoff=0.001 ### convergence criterion for EM
  
  n.normal<-dim(exprM.normal)[1] # -- number sample normal
  n.tumor<-dim(exprM)[1] # -- number sample normal
  exprM<-t(scale((rbind(exprM,exprM.normal))))
  geneNames<-rownames(exprM);
  curgeneNames<-rownames(exprM);
  if(estimateP)
  {
    for(j in 1:1000)  # --- over EM iterations 
      # --- (CHECK: why do we need this loop since we have (A Loop) and (B Loop))
    {
      if(j==1){
        cur.P<-c(purity,rep(0,n.normal));
        cur.b<-mean(cur.P[seq(1,n.tumor)])*(1-mean(cur.P[seq(1,n.tumor)]))/var(cur.P[seq(1,n.tumor)])-1
      }else{
        print(paste("Epurity",j,":",cor(purity,cur.P[seq(1,n.tumor)])))
      }
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
      
      paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
      x.v=exprM[pi,] 
      for(i in 1:10000) # -- (A Loop) over EM iterations
       {
         if(i==1) {rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)}
         r<-rnew;
         temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
         rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
         if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01)){#print(paste("converge in step",i));
         break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
print("EM1d")
before.P<-cur.P;
#print(before.P)
## -- Estimate purity (P) and delta parameter in beta distribution
for(i in 1:10000) {                          
  
  if (betaModel==FALSE){ 
    new.P<-sapply(1:n.tumor,function(ni){
      f<-function(x){
        V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
        U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
        0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dnorm(log(purity[ni]/(1-purity[ni])),log(x/(1-x)),sqrt(cur.b))) # purity
      }
      optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
    new.b<-sum((log(purity[seq(1,n.tumor)]/(1-purity[seq(1,n.tumor)]))-log(new.P/(1-new.P)))^2)/n.tumor
    
  } else {
    new.P<-sapply(1:n.tumor,function(ni){
    f<-function(x) {
       V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
       U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
       0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)- log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
      }
    optim(cur.P[ni],f,lower=0,upper=1,method="Brent")$par})
    f<-function(x) -sum(log(dbeta(purity,new.P*x,(1-new.P)*x)))
    new.b<-optim(cur.b,f,lower=0)$par        
  }
  
  if(all(abs(c(new.P-cur.P[seq(1,n.tumor)]))<0.001)){#print(paste("converge in step",i));
    break()}
    cur.b<-new.b;
    cur.P<-c(new.P,rep(0,n.normal));}

print(paste("cur.b",cur.b))

if(all((before.P-cur.P)<0.001))break();
print(mean(abs(cur.P-before.P)))
    }

if(j==10000){ stop("P not converge"); convergence="P not converge"}else{
  print("Epurity converged"); convergence="P converged"}

Epurity<-cur.P
  }else{Epurity<-purity};

time<-(proc.time()-t)
return(list(Epurity,paraM,convergence,time))

}


deNet<-function(exprM,purity,paraM,rhoy,rhoz,label="test",estimateP=FALSE,bootstrap=0,
                CV=10,CVIndex=0,insilico="NO",BIC=TRUE,DIR)
  {

# -- exprM: matrix of expression (p x n) whre p is the number of genes and n the number of samples
# --- purity n-dimensional vector containing estimated purity  
setwd(DIR)
  source("Codes/Cell_Mix_Function_Francesca.r")

t<-proc.time() 
library("MASS")
library("mvtnorm")
library("clusterGeneration")
library("glasso")
library("doMC")
registerDoMC(10)
library("foreach")
options(warn=-1)
convergecutoff=0.1 ### convergence criterion for EM
if(bootstrap>0){
set.seed(bootstrap);
rindex<-sample(1:ncol(exprM),ncol(exprM),replace=TRUE)
exprM<-exprM[,rindex];
purity<-purity[rindex];
geneNames<-rownames(exprM);
exprM<-exprM[apply(exprM,1,sd)!=0,];
curgeneNames<-rownames(exprM);
###standardization
exprM<-t(scale(t(exprM)))
}
geneNames<-rownames(exprM);
curgeneNames<-rownames(exprM);

if(CVIndex>0)
  {
set.seed(0);
rIndex<-sample(1:ncol(exprM));
size<-floor(ncol(exprM)/CV);
left<-ncol(exprM)%%CV;
if(left>0)index<-c(rep(1:CV,size),1:left);
if(left==0)index<-rep(1:CV,size)
testIndex<-rIndex[which(index==CVIndex)];
trainIndex<-rIndex[which(index!=CVIndex)];

exprM<-exprM[apply(exprM[,trainIndex],1,sd)!=0,];
curgeneNames<-rownames(exprM);
###standardization
exprM<-t(scale(t(exprM)))
  }

if (CVIndex==0) exprM<-t(scale(t(exprM)))

if(estimateP){
for(j in 1:1000)  # --- over EM iterations 
   {
   
if(j==1){
  cur.P<-purity;
  cur.b<-mean(cur.P)*(1-mean(cur.P))/var(cur.P)-1
  }else{

print(paste("Epurity",j,":",cor(purity,cur.P)))
  }

## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)

paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
  {
        x.v=exprM[pi,]  
             for(i in 1:10000) # -- (A Loop) over EM iterations
               {
                 if(i==1)
                {
                  rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
                    }
                  r<-rnew;
                  temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
                  rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
               if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
               {#print(paste("converge in step",i));
                 break()}
               }
               c(r$uy,r$uz,r$vy,r$vz)
      }

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
print("EM1d")
before.P<-cur.P;

## -- Estimate purity (P) and delta parameter in beta distribution
            for(i in 1:10000)  # -- (B Loop) CHECK this loop should be the same as (A Loop)
             {                          
              new.P<-sapply(1:ncol(exprM),function(ni){
                f<-function(x)
                    {
                  V<-paraM[,"vy"]*x^2+paraM[,"vz"]*(1-x)^2;
                  U<-paraM[,"uy"]*x+paraM[,"uz"]*(1-x);
                  
                  # ** CHECK **
                  0.5*sum(log(V))+0.5*sum((exprM[,ni]-U)^2/V)-
                    log(dbeta(purity[ni],x*cur.b,(1-x)*cur.b)) # purity
              }
                 optim(cur.P[ni],f,lower=0.1,upper=0.9,method="Brent")$par})


      f<-function(x)sum(-log(dbeta(purity,new.P*x,(1-new.P)*x))) # CHECK: purity and new.P should be inverted??
      new.b<-optim(cur.b,f,lower=0)$par
              if(all(abs(c(new.P-cur.P))<0.001))
               {#print(paste("converge in step",i));
                 break()}
       cur.b<-new.b;
       cur.P<-new.P
       print(cur.b)
       print(mean(abs(cur.P-purity)))
            }
print(paste("cur.b",cur.b))


if(all((before.P-cur.P)<0.001))break();
 print(mean(abs(cur.P-before.P)))
}

if(j==10000){ stop("P not converge")}else{print("Epurity converged")}
                
Epurity<-cur.P
 

}else{ 
  # ---- ************************************************************************ ----#
  ## -- Estimate u and v using univariate model
  paraM<-foreach(pi= 1:nrow(exprM))%dopar%  # -- for each gene
{
  x.v=exprM[pi,]  
  for(i in 1:10000) # -- (A Loop) over EM iterations
  {
    if(i==1)
    {
      rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
    }
    r<-rnew;
    temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,purity)
    rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
    if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
    {break()}
  }
  c(r$uy,r$uz,r$vy,r$vz)
}

paraM<-do.call("rbind",paraM)
colnames(paraM)<-c("uy","uz","vy","vz")
print("EM1d")
Epurity<-purity

# ---- ************************************************************************ ----#
};

exprM<-t(exprM)
## ---  Analyze through Multiv-EM
## --- Estimate mean and covariance matrix for y and z
  LL.temp=c()
  i<-1
  while(i<5000)
  {
   print(paste("Multiv-EM",i))
   if(i==1)
     {
       # -- initialize values
       rnew<-list(Uy=paraM[,"uy"],Uz=paraM[,"uz"],Sigmay=diag(paraM[,"vy"]),Sigmaz=diag(paraM[,"vz"]))
       LL.temp<-calculateLL(rnew$Uy,rnew$Uz,rnew$Sigmay,rnew$Sigmaz,c=1,exprM,Epurity)
     }
   r<-rnew;
   cur.para<-Estep.noloop(r$Uy,r$Uz,r$Sigmay,r$Sigmaz,c.scale=1,exprM,Epurity)
  
   rnew<-Mstep.v2(cur.para$Ey,cur.para$Ez,cur.para$Esigmay,cur.para$Esigmaz,rhoy=rhoy, rhoz=rhoz)

    # --- remove bad nodes
    v<-cbind(diag(rnew$Sigmay),diag(rnew$Sigmaz))
    badpoints<-which(v[,1]>20|v[,2]>20|is.na(v[,1])|is.na(v[,2]))

    if(length(badpoints)>0)
      {
        print(paste("remove bad points",badpoints));
        exprM<-exprM[,-c(badpoints)];
        curgeneNames<-curgeneNames[-c(badpoints)];
        rnew$Uy<-rnew$Uy[-c(badpoints)];
        rnew$Uz<-rnew$Uz[-c(badpoints)];
        rnew$Sigmay<-rnew$Sigmay[-c(badpoints),-c(badpoints)];
        rnew$Sigmaz<-rnew$Sigmaz[-c(badpoints),-c(badpoints)];

        r$Uy<-r$Uy[-c(badpoints)];
        r$Uz<-r$Uz[-c(badpoints)];
        r$Sigmay<-r$Sigmay[-c(badpoints),-c(badpoints)];
        r$Sigmaz<-r$Sigmaz[-c(badpoints),-c(badpoints)];
#        next;
         }
      
     LL.temp<-c(LL.temp, calculateLL(rnew$Uy,rnew$Uz,rnew$Sigmay,rnew$Sigmaz,c=1,exprM,Epurity))
     print(paste("edgeN",sum(abs(ginv(rnew$Sigmay)[upper.tri(rnew$Sigmay)])>0.001)))
     print(paste("edgeN",sum(abs(ginv(rnew$Sigmaz)[upper.tri(rnew$Sigmaz)])>0.001)))
#     print(paste("LL",LL.temp[length(LL.temp)]))
    
    # if(all(abs(c(rnew$Uy-r$Uy,rnew$Uz-r$Uz))<convergecutoff)&all(abs((rnew$Sigmay-r$Sigmay))<convergecutoff)&all(abs((rnew$Sigmaz-r$Sigmaz))<convergecutoff))
     
# if(all(abs(c(rnew$Uy-r$Uy,rnew$Uz-r$Uz))<convergecutoff))
     if(abs(LL.temp[i+1]-LL.temp[i])<convergecutoff) # -- break when log likelihood converges
     { break()}else{
      print(abs(LL.temp[i+1]-LL.temp[i]))
#      print(max(abs(rnew$Sigmay-r$Sigmay)))
#      print(max(abs(rnew$Sigmaz-r$Sigmaz)))
 #     print(sum(abs(rnew$Sigmay)))
 #     print(sum(abs(rnew$Sigmaz)))      
# if(i==50){save(cur.para,r,exprM,Epurity,rho,file=paste(label,"debug.rda",sep="_"));stop()}
    }
    i<-i+1;
    
  }
   if(i==5000) {converge<-FALSE}else{converge<-TRUE}

#if(CVIndex>0)
#   {
#     exprMorig<-exprMorig[curgeneNames,]
#     LLtest<-calculateLL(r$Uy,r$Uz,r$Sigmay,r$Sigmaz,c=1,t(exprMorig[curgeneNames,testIndex]),Epurityorig[testIndex]);
#     LLtrain<-calculateLL(r$Uy,r$Uz,r$Sigmay,r$Sigmaz,c=1,t(exprMorig[curgeneNames,trainIndex]),Epurity);

#     Sigmax<-glasso(cov(exprM),rhoy)$w;
#     Ux<-apply(exprMorig[curgeneNames,],1,mean);
#     LLtrain0<-sum(sapply(1:length(trainIndex),function(i){dmvnorm(exprMorig[curgeneNames,trainIndex[i]], Ux, Sigmax, log = TRUE)}))
#     LLtest0<-sum(sapply(1:length(testIndex),function(i){dmvnorm(exprMorig[curgeneNames,testIndex[i]], Ux, Sigmax, log = TRUE)}))  
#   }

if(length(curgeneNames)!=length(geneNames)){
  t<-rep(NA,length(geneNames));
  names(t)<-geneNames;
  t[curgeneNames]<-r$Uy;
  r$Uy<-t;
  t<-rep(NA,length(geneNames));
  names(t)<-geneNames;
  t[curgeneNames]<-r$Uz;
  r$Uz<-t;
  t<-matrix(NA,length(geneNames),length(geneNames),dimnames=list(geneNames,geneNames));
  t[curgeneNames,curgeneNames]<-r$Sigmay;
  r$Sigmay<-t;
  t<-matrix(NA,length(geneNames),length(geneNames),dimnames=list(geneNames,geneNames));
  t[curgeneNames,curgeneNames]<-r$Sigmaz;
  r$Sigmaz<-t;
}
 

time=(proc.time()-t)
# if(CVIndex>0)   {
#     save(r,converge,Epurityorig,purity,exprMorig,LLtest,LLtrain,LLtrain0,LLtest0,trainIndex,testIndex,time,file=paste(paste("~/work/mixCell5/deNetCV",label,rhoy,rhoz,CVIndex,sep="_"),".rda",sep=""))
#   }else{

#save(r,edge.y,edge.z,converge,Epurity,exprM,purity,time,LL.All.Data,file=paste(paste("InSilico.highdim/",label,rhoy,rhoz,bootstrap,sep="_"),".rda",sep=""))
#   }
edge.y<-sum(abs(ginv(rnew$Sigmay)[upper.tri(rnew$Sigmay)])>0.001)
edge.z<-sum(abs(ginv(rnew$Sigmaz)[upper.tri(rnew$Sigmaz)])>0.001)
BIC<-(-2*LL.temp[i]+(edge.y+edge.z)*log(dim(exprM)[1]))

return(list(r=r,LL.temp,list(BIC=BIC,LL=LL.temp[i],EY=edge.y,EZ=edge.z),
            converge=converge,purity=Epurity,data=exprM,time=time))
 }



deNet.diagonal<-function(exprM,purity,label="test",estimateP=FALSE,bootstrap=0,
                CV=10,CVIndex=0,insilico="NO",BIC=TRUE,DIR)
{
  
  # -- exprM: matrix of expression (p x n) whre p is the number of genes and n the number of samples
  # --- purity n-dimensional vector containing estimated purity  
  setwd(DIR)
  source("Codes/Cell_Mix_Function_Francesca.r")
  
  t<-proc.time() 
  library("MASS")
  library("mvtnorm")
  library("clusterGeneration")
  library("glasso")
  library("doMC")
  registerDoMC(10)
  library("foreach")
  options(warn=-1)
  convergecutoff=0.1 ### convergence criterion for EM
  geneNames<-rownames(exprM);
  curgeneNames<-rownames(exprM);
  


        cur.P<-purity;
      
      ## -- Estimate u and v (mean and variance of latent variables, i.e., uy uz vy vz)
     converge=FALSE
       paraM<-converge.step<-NULL
       
       for (pi in 1:nrow(exprM))  # -- for each gene
      {
        x.v=exprM[pi,]  
        for(i in 1:10000) # -- (A Loop) over EM iterations
        {
          if(i==1)
          {
            rnew<-list(uy=mean(x.v),uz=mean(x.v),vy=var(x.v),vz=var(x.v),c=1)
          }
          r<-rnew;
          temp<-Estep.1d(r$uy,r$uz,r$vy,r$vz,1,x.v,cur.P)
          rnew<-Mstep.1d(temp$Ey,temp$Ey2,temp$Ez, temp$Ez2)
          if(all(abs(c(rnew$uy-r$uy,rnew$uz-r$uz,rnew$vy-r$vy,rnew$vz-r$vz))<0.01))
          {#print(paste("converge in step",i));
            converge.step=c(converge.step,i)
            break()}
          
        }
        paraM<-rbind(paraM,c(r$uy,r$uz,r$vy,r$vz))
      }
      
      colnames(paraM)<-c("uy","uz","vy","vz")

  return(list(r=paraM,purity=purity,data=exprM,converge=converge.step))
}



