##############################################################################################################################################
###   AMM   -- R script for GWAS corecting for population structure (similar to EMMAX and P3D)
###
#######
#
##
## 
##		
#


##REQUIRED DATA & FORMAT


#library(doMC) # only required for multi-core
#library(gaston) # only required if AI algorithm is used


#requires functions from the original emma function (Kang et al. 2008, Genetics) to calculate the exact p-value if update.top_snps=TRUE or if AI=FALSE 
# source ('emma.r')  needs to be sourced before.
#PHENOTYPE - Y: a n by m matrix, where n=number of individuals and the rownames(Y) contains the individual names

#GENOTYPE - X: a n by m matrix, where n=number of individuals, m=number of SNPs, with rownames(X)=individual names, and colnames(X)=SNP names
#KINSHIP - K: a n by n matrix, with rownames(K)=colnames(K)=individual names, you can calculate K as IBS matrix using the emma package K<-emma.kinship(t(X)) 
# or use the vanRaden Kinship matrix from the function G.matrix in the R-package snpReady
#each of these data being sorted in the same way, according to the individual name
#
#
#SNP INFORMATION - SNP_INFO: a data frame having at least 3 columns:
# - 1 named 'SNP', with SNP names (same as colnames(X)),
# - 1 named 'Chr', with the chromosome number to which belong each SNP
# - 1 named 'Pos', with the position of the SNP onto the chromosome it belongs to.
#######
###WARNINGS####

# mc enabled to use on linux 

amm_gwas<-function(Y,X,K,p=0.001,m=2,run=TRUE,AI=TRUE,calculate.effect.size=FALSE,include.lm=FALSE,include.kw=FALSE,use.SNP_INFO=FALSE,update.top_snps=FALSE,report=TRUE,plot_h=FALSE, mc = FALSE, cores ='all') {

    stopifnot(is.numeric(Y[,1]))
    Y_ <- Y[order(Y[,1]),]
    Y <- as.matrix(Y_[,m])
    rownames(Y) <- as.integer(Y_[,1])
    Y <- na.omit(Y)
    XX <- X[rownames(X) %in% rownames(Y),]
    XX<-XX[order(as.numeric(rownames(XX))),]
    if(report == TRUE) {
        cat ('GWAS performed on', length(which(rownames(Y)%in%rownames(X))),'ecotypes, ', nrow(Y)-length(which(rownames(Y)%in%rownames(X))),'values excluded','\n')
    }

    if(use.SNP_INFO == FALSE) { 
        options(stringsAsFactors = FALSE)
        if(report == TRUE) {
            cat('SNP_INFO file created','\n')
        }
        SNP_INFO<-data.frame(cbind(colnames(X),matrix(nrow=ncol(X),ncol=2,data=unlist(strsplit(colnames(X),split='- ')),byrow=TRUE)))
        colnames(SNP_INFO) <-c('SNP','Chr','Pos')
        SNP_INFO[,2]<-as.numeric(SNP_INFO[,2])
        SNP_INFO[,3]<-as.numeric(SNP_INFO[,3])
    } else {
        if(report == TRUE) { cat('User definied SNP_INFO file is used','\n')}}


    Y1 <- as.matrix(Y[rownames(Y) %in% rownames(XX),])
    colnames(Y1) <- colnames(Y)
    rownames(Y1) <- rownames(Y)[rownames(Y)%in%rownames(XX)]
    ecot_id <- as.integer(rownames(Y1))

    K1 <- K[rownames(K) %in% ecot_id,]
    K2 <- K1[,colnames(K1) %in% ecot_id]
    K_ok <- as.matrix(K2)
    K_ok<-K_ok[order(as.numeric(rownames(K_ok))),order(as.numeric(colnames(K_ok)))]

    a <- rownames(K_ok)

    n <- length(a)
    K_stand <- (n-1)/sum((diag(n)-matrix(1,n,n)/n)*K_ok)*K_ok

    Y <- Y1[which(rownames(Y1)%in%a),]
    X_ <- XX[which(rownames(XX)%in%a),]

    rm(X,XX)
    gc()

    
    AC_1 <- data.frame(colnames(X_),apply(X_,2,sum))
    colnames(AC_1) <- c('SNP','AC_1')
    MAF_1 <- data.frame(AC_1,AC_0=nrow(X_)-AC_1$AC_1)
    MAF_2 <- data.frame(MAF_1,MAC=apply(MAF_1[,2:3],1,min))
    MAF_3 <- data.frame(MAF_2,MAF=(MAF_2$MAC/nrow(X_)))
    MAF_ok <- merge(SNP_INFO,MAF_3,by='SNP')
    rm(AC_1,MAF_1,MAF_2,MAF_3)

#Filter for MAF

    MAF <- subset(MAF_ok,MAF<p)[,1]
    X_ok <- X_[,!colnames(X_) %in% MAF]
    rm(MAF,X_)
    
# REML

if (AI==TRUE) {
  require(gaston)
    null <-lmm.aireml(Y,K=K_stand)
    herit<-null$tau/(null$tau+null$sigma2)
    M <- solve(chol(null$tau*K_stand+null$sigma2*diag(dim(K_stand)[1])))
    H2 <- seq(0,1,length=101)
    lik <- lmm.diago.likelihood(h2 = H2, Y = Y, eigenK = eigen(K_stand))
   minh<-min(H2[which(exp(lik$likelihood-max(lik$likelihood))>0.9)])
   maxh<-max(H2[which(exp(lik$likelihood-max(lik$likelihood))>0.9)])
   # covH1<-matrix(nrow=2,ncol=2,data=c(null$vg,null$ve,null$ve,null$vg))
  #  se_h1<-deltamethod( ~ x1/(x1+x2),c(null$vg,null$ve),covH1)
      
      if (report==TRUE) {
        cat('pseudo-heritability estimate is',round(herit,digits=3),'with a 90% CI between',minh,'and',maxh,'.\n')
      }
   if(plot_h==TRUE) {
     pdf(file='heritability_likelihood.pdf')
         plot(H2, exp(lik$likelihood-max(lik$likelihood)), type="l", ylab="relative likelihood")
        points(herit,exp(null$logL-max(lik$likelihood)),pch=16,col='red')
        abline(h=0.90,lty=2,col='red')
        cat('A heritability likelihood plot is generated','\n')
      dev.off()
   }
} else {
  Xo<-rep(1,nrow(X_ok))
  ex<-as.matrix(Xo)
  null<-emma.REMLE(Y,ex,K_stand)
  
  herit<-null$vg/(null$vg+null$ve) 
  M <- solve(chol(null$vg*K_stand+null$ve*diag(dim(K_stand)[1])))
  
    if (report==TRUE) {
      cat('pseudo-heritability estimate is',round(herit,digits=3),'\n')
    }
   
}
      if (run==FALSE&AI==TRUE) {
        
        cat('no GWAS performed','\n')
  
        return(c(herit,minh,maxh))   
      }
if (run==FALSE&AI==FALSE) {
       return(herit)}
   if (run==TRUE) {
      
        Y_t <- crossprod(M,Y)
        int_t <- crossprod(M,(rep(1,length(Y))))
        
        if(calculate.effect.size == TRUE) {
            models1 <- apply(X_ok,2,function(x){summary(lm(Y_t~0+int_t+crossprod(M,x)))$coeff[2,]})
            bet <- models1[1,]
            se <- models1[2,]
### variance explained from betas veb  
                                        # bet^2*var(X_ok[,t])/var(Y[,n]) =  veb/(1-veb)
            veb <- bet[1]^2*var(X_ok[,1])/var(Y)/(1+bet[1]^2*var(X_ok[,1])/var(Y))
            for( t in 2:ncol(X_ok)) {
                veb <- c(veb,bet[t]^2*var(X_ok[,t])/var(Y)/(1+bet[t]^2*var(X_ok[,t])/var(Y)))
            }
          
# similar (but slightly different) to RSS method below or to adjusted R2 from lm call 
            out_models<-data.frame(SNP=colnames(models1),Pval=models1[4,],variance_explained=veb,beta=bet,se_beta=se)

        } else {
             #EMMAX SCAN
            RSS_env <- rep(sum(lsfit(int_t,Y_t,intercept = FALSE)$residuals^2),ncol(X_ok))
            
            if(mc == FALSE) {
              
                R1_full <- apply(X_ok,2,function(x){sum(lsfit(cbind(int_t,crossprod(M,x)),Y_t,intercept = FALSE)$residuals^2)})
            
            } else {
                require(doMC)
                if(cores == "all") { 
                    registerDoMC(system("cat /proc/cpuinfo  | grep processor  | wc -l", intern = TRUE))
                } else {
                    registerDoMC(cores)
                }
                    R1_full <- foreach(i = 1:ncol(X_ok), .combine = "c") %dopar%
                    sum(lsfit(cbind(int_t,crossprod(M,X_ok[,i])),Y_t,intercept = FALSE)$residuals^2)
            
            }
            
            pa<-nrow(Y1)
            
            F_1 <- ((RSS_env-R1_full)/1)/(R1_full/(pa-3))
            pval_Y1 <- pf(F_1,1,(pa-2),lower.tail=FALSE)
            
            snp <- colnames(X_ok)
            out_models <- data.frame(SNP=snp,Pval=pval_Y1,variance_explained=1-R1_full/RSS_env)
        }
        
        output <- merge(MAF_ok,out_models,by='SNP')
        
        if(include.lm == TRUE) {
            RSS_env_ <- rep(sum(lsfit(rep(1,length(Y)),Y,intercept = FALSE)$residuals^2),ncol(X_ok))
            R1_full_ <- apply(X_ok,2,function(x){sum(lsfit(x,Y,intercept = TRUE)$residuals^2)})
            pa<-nrow(Y1)
            
            F_1_ <- ((RSS_env_-R1_full_)/1)/(R1_full_/(pa-3))
            pval_Y1_lm <- pf(F_1_,1,(pa-3),lower.tail=FALSE)
            
            snp <- colnames(X_ok)
            out_models_lm <- data.frame(SNP=snp,Pval_lm=pval_Y1_lm)
            output <- merge(output,out_models_lm,by='SNP')
        }

        if(include.kw == TRUE) {
            kw <- apply(X_ok,2,function(x){kruskal.test(Y,x)$p.value})
            KW <- cbind(names(kw),kw)
            colnames(KW) <- c('SNP','Pval_kw')
            output <- merge(output,KW,by='SNP')
            
        }
        
        ## update tp SNPs with correct model 
        if(update.top_snps == FALSE) {
            return(output)
        } else {
            oi <- output[order(output[,8]),][1:update.top_snps,1]
            xs <- t(X_ok[,colnames(X_ok)%in%oi])
            auto <- emma.ML.LRT (Y, xs, K_ok)
            up <- data.frame(SNP=rownames(xs),update=auto$ps[,1])
            
            for(i in 1:update.top_snps) {
                output[which(output[,1]==up[i,1]),]$Pval <- up[i,2]
            }
            return(output)
        }
    }

}





