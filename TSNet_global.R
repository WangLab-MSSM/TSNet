
  setwd("/sc/orga/projects/PGDAC/francesca/TSNet/Codes")
  source("Cell_Mix_Function_Francesca.r")
  source("network_function.txt")

  setwd("/sc/orga/projects/PGDAC/francesca/TSNet/ccRCC/")
  rna<-read.table("data/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB_clean_1105.tsv",header=TRUE,sep="\t",stringsAsFactors=FALSE,row.names=1)
  tissue<-as.character(rna[2,])
  rna<-rna[-c(1,2),]
  
  # -- load purity
  purity<-read.table("data/Global_signature_nonimputed.tsv",skip=2,header=TRUE,sep="\t",row.names=1)
  mg<-match(colnames(purity),colnames(rna))
  
  rna<-rna[,mg[!is.na(mg)]]
  tissue<-tissue[mg[!is.na(mg)]]
  purity<-purity["TumorPurity",!is.na(mg)]

  # -- delete missing values
  rna<-rna[rowSums(is.na(rna))==0,]
  
  X<-t(rna)
  index<-(apply(X[tissue=="Tumor",],2,sd)>1E-1 & apply(X[tissue=="Normal",],2,sd)>1E-1)
  X<-X[,index]
  geneID<-rownames(rna)[index]
  mean.x<-apply(X,2,mean)
  sd.x<-apply(X,2,sd)
  # --- Find purity considering the whole data
  # --- Check code Purity_fit.R 
   
     # --- Cross Validation to find best penalty parameters
     class<-rep(0,dim(X)[1])
     class[tissue=="Tumor"]<-1
     class[tissue=="Normal"]<-2
     
     X<-apply(X,2,as.numeric)
     X<-apply(X,2,function(x)(x-mean(x))/sd(x))
     P<-as.numeric(purity); P[P<.01]<-.01; P[P>.99]<-.99
     
     X1<-X[class==1,]; P1<-P[class==1]
     X2<-X[class==2,]; P2<-P[class==2]
     
     out1<-deNet.only.mean(X1,P1,betaModel=TRUE)
     out2<-deNet.only.mean(X2,P2,betaModel=TRUE)
     
     bulk1<-colSums(X1)/dim(X1)[1]
     bulk2<-colSums(X2)/dim(X2)[1]
     save(mean.x,sd.x,geneID,out1,out2,bulk1,bulk2,file=paste0("/sc/orga/projects/PGDAC/francesca/TSNet/Results/ccRCC/global/TSNet_global.rda"))
