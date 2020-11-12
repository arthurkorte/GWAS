### Function to simulate traits with causative markers for GWAS 

# need genomic data . A sample of 22k random SNPs from the Arabidopsis 1001G Projekt is available as example for Arabidopsis traits 
# a file with the Information of all 2029 accessions for which genotype data exist is also provided (lat_long_2029.csv). This information has been downloaded from the AraPheno database 

#A<-read.csv('~/git/GWAS/data/lat_long_2029.csv')
#If you use own data ensure that the id should be in a column called 'accession_id'
#load('~/git/GWAS/data/X_Mix.rda')
#rownames(X) needs to be the ids of the accessions
# for models with multiple causative markers, it is recommended to have them orderd by the amount of variance explained (ve)

# the simulation create random, but reproducible (seed) phenotypes.
# to simulate only in specific accession provide an integer of accessions ids in sp_acc (length need to be >1). 
# to siumulate distinct markers, the script can be modified manually.

sim_Y<-function(n=100,acc=A,sp_acc=0,no_acc=200,fix_acc=TRUE,SNPs=X,no_snps=1,ve=.1,mac=5,h2=0.7,seed=42,bk=1000) {
  stopifnot(no_snps>0)
  stopifnot(length(ve)==no_snps)
  set.seed(seed)
  
  Sim<-list()
  Caus<-list()
if (length(sp_acc)>1){
  a<-sp_acc
  no_acc=length(a)
}else {
a<-sample(acc$accession_id,no_acc) }
  
X_<-SNPs[rownames(SNPs)%in%a,]

af<-apply(X_,2,sum)
X_ok<-X_[,which(af>mac&af<(no_acc-mac))]

u<-1
# set the seed 
set.seed(seed+u)
  
caus<-X_ok[,sample(1:ncol(X_ok),(no_snps+1))]

#generating polygenic background

X3<-X_ok[,!colnames(X_ok)%in%colnames(caus)]

back<-X3[,sample(1:ncol(X3),bk)]

betas<-rnorm(bk,mean=0,sd=0.1)
first<- back %*% betas
### adding genetic background to data
sim<-data.frame(ecot_id=as.integer(rownames(back)),value=first)
### set heritability 
dat<-var(sim[,2])

h_2<-dat/h2-dat
fix1<-rnorm(nrow(back),0,sqrt(h_2))
sim_<-data.frame(ecot_id=as.integer(rownames(back)),value=first+fix1)


for ( t in 1:length(ve)) {
beta<-sqrt((ve[t]/(1-ve[t]))*(var(sim_[,2])/var(caus[,t])))

cand<-beta*caus[,t]
sim_$value<-sim_$value+cand

}
Sim[[u]]<-sim_
Caus[[u]]<-colnames(caus)[1:t]

if (fix_acc==FALSE) {
  for ( u in 2:n ) {
    set.seed(seed+u)
  a<-sample(A$accession_id,no_acc)
  X_<-subset(X,rownames(X)%in%a)
  
  af<-apply(X_,2,sum)
  X_ok<-X_[,which(af>mac&af<(no_acc-mac))]
  

  caus<-X_ok[,sample(1:ncol(X_ok),(no_snps+1))]
  
  #generating polygenic background
  
  X3<-X_ok[,!colnames(X_ok)%in%colnames(caus)]
  
  back<-X3[,sample(1:ncol(X3),bk)]
  
  betas<-rnorm(bk,mean=0,sd=0.1)
  first<- back %*% betas
  ### adding genetic background to data
  sim<-data.frame(ecot_id=as.integer(rownames(back)),value=first)
  ### set heritability 
  dat<-var(sim[,2])
  
  h_2<-dat/h2-dat
  fix1<-rnorm(nrow(back),0,sqrt(h_2))
  sim_<-data.frame(ecot_id=as.integer(rownames(back)),value=first+fix1)
  
  
  for ( t in 1:length(ve)) {
    beta<-sqrt((ve[t]/(1-ve[t]))*(var(sim_[,2])/var(caus[,t])))
    
    cand<-beta*caus[,t]
    sim_$value<-sim_$value+cand
    
  }
  Sim[[u]]<-sim_
  Caus[[u]]<-colnames(caus)[1:t]
}
  }else {
  for ( u in 2:n ) {
    set.seed(seed+u)
 
  caus<-X_ok[,sample(1:ncol(X_ok),(no_snps+1))]
  
  #generating polygenic background
  
  X3<-X_ok[,!colnames(X_ok)%in%colnames(caus)]
  
  back<-X3[,sample(1:ncol(X3),bk)]
  
  betas<-rnorm(bk,mean=0,sd=0.1)
  first<- back %*% betas
  ### adding genetic background to data
  sim<-data.frame(ecot_id=as.integer(rownames(back)),value=first)
  ### set heritability 
  dat<-var(sim[,2])
  
  h_2<-dat/h2-dat
  fix1<-rnorm(nrow(back),0,sqrt(h_2))
  sim_<-data.frame(ecot_id=as.integer(rownames(back)),value=first+fix1)
  
  
  for ( t in 1:length(ve)) {
    beta<-sqrt((ve[t]/(1-ve[t]))*(var(sim_[,2])/var(caus[,t])))
    
    cand<-beta*caus[,t]
    sim_$value<-sim_$value+cand
    
  }
  Sim[[u]]<-sim_
  Caus[[u]]<-colnames(caus)[1:t]
}
}

return(list(Y=Sim,Caus=Caus))
}
