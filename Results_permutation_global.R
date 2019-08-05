# --- load results
library(ggplot2)
Sim=2

DIR<-"/Users/francescapetralia/Box Sync/Kidney/Revision/DEG"
DIR<-"/sc/orga/projects/PGDAC/francesca/TSNet/"

# -- function FDR
FDR<-function(perm,score,TH,FC){
  score_or<-sort(score,decreasing=TRUE)
  
    for (j in 5:dim(perm)[1]){
    fp<-sum(perm>score_or[j])/dim(perm)[2]
    fdr<-fp/j
    if (fdr > TH) break;
    print(j)
  }
  index<-(score>score_or[j] & FC>1)
  diff.estimated<-seq(1,dim(perm)[1])[index]
  return(diff.estimated)
}
  
# --- load results
all<-list.files(path=paste0(DIR,"Results/ccRCC/global/"),pattern="Permutation*",full.names = TRUE)

perm.bulk1<-perm.bulk2<-Uy1.tsnet<-Uz1.tsnet<-Uy2.tsnet<-Uz2.tsnet<-NULL
for (j in 1:length(all)){
  load(all[j])
  Uy1.tsnet<-cbind(Uy1.tsnet,out1[[2]][,"uy"])
  Uz1.tsnet<-cbind(Uz1.tsnet,out1[[2]][,"uz"])
  
  Uy2.tsnet<-cbind(Uy2.tsnet,out2[[2]][,"uy"])
  Uz2.tsnet<-cbind(Uz2.tsnet,out2[[2]][,"uz"])
  
  perm.bulk1<-cbind(perm.bulk1,bulk1)
  perm.bulk2<-cbind(perm.bulk2,bulk2)
}
  
# --- compute mean difference based on permutation
perm.tsnet<-(Uy1.tsnet-Uy2.tsnet)
perm.bulk<-(perm.bulk1-perm.bulk2)

# -- load tsnet original
load(paste0("/sc/orga/projects/PGDAC/francesca/TSNet/Results/ccRCC/global/TSNet_global.rda"))
diff.tsnet<-(out1[[2]][,"uy"]-out2[[2]][,"uy"])
diff.bulk<-(bulk1-bulk2)

# --- compute differential edges
TH=.1 # -- FDR cut-off
tsnet.diff<-FDR(perm.tsnet,diff.tsnet,TH=TH,FC=(out1[[2]][,"uy"]*sd.x+mean.x)-(out2[[2]][,"uy"]*sd.x+mean.x))
bulk.diff<-FDR(perm.bulk,diff.bulk,TH=TH,FC=(bulk1*sd.x+mean.x)-(bulk2*sd.x+mean.x))

length(geneID[tsnet.diff])
length(geneID[bulk.diff])

# --- list of genes upregulated 
write.table(geneID[tsnet.diff],file=paste0(DIR,"Results/ccRCC/Global_TSnet_upregulated_tumor_th_",TH,".tsv"),sep="\t",append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)
write.table(geneID[bulk.diff],file=paste0(DIR,"Results/ccRCC/Global_Original_upregulated_tumor_th_",TH,".tsv"),sep="\t",append=FALSE,quote=FALSE,col.names=FALSE,row.names=FALSE)
