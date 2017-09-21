t2029_par<-function(Y,n=2,incl.lm=FALSE,incl.beta=FALSE,save.output=F,generate.plot=F,cores=11) {
    
    require("doMC")
    registerDoMC(cores=cores)
    source('emma.r')
    source('gwas.r')
    source('plots_gwas.r')
    load('K_2029.rda')

    K<-K_2029
    
    results <- foreach (u = 1:22) %dopar% {
        if(u != 2){
            filename<-paste('X_2029_',u,'.rda',sep='')
            load(filename)
             amm_gwas(Y,X,K,m=n,include.lm=incl.lm,calculate.effect.size=incl.beta)
        }else{
            filename<-paste('X_2029_',u,'.rda',sep='')
            load(filename)
            amm_gwas(Y,X,K,m=n,include.lm=incl.lm,calculate.effect.size=incl.beta)
        }
    }
    
    results <- do.call("rbind",results)
    return(results)
    
    
    if (save.output==T) {

name1<-paste(colnames(Y)[n],'_gwas2029.rda',sep='')
save(results,file=name1)
}
if (generate.plot==T) {
 name2<-paste(colnames(Y)[n],'_gwas2029.pdf',sep='')
pdf(file=name2)
plot_gwas(results)
dev.off()
}

}







