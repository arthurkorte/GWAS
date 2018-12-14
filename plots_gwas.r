## function for plotting GWAS results, default filter is MAF 5%, to change for MAC filter change maf_or_mac to a different value then 1

plot_gwas<-function(output,h=8,maf_or_mac=1,maf=0.05,mac=NA,black=T,thres=NA,max.y=NA,show=NA,lower.limit=0.01,title=TRUE,name=NA) {
 colnames(output)[h]<-'Plot'
if(maf_or_mac==1) {
mafi=paste('maf',maf)
new_<-subset(output,output$MAF>=maf&output$Plot<lower.limit)
if (is.na(thres)==T) {thres<--log10(0.05/nrow(subset(output,output$MAF>=maf)))}

if (is.na(mac)==FALSE) {cat ('MAF Filter of',maf,' have been used, for MAC set maf_or_mac to 2','\n')}}
 else { 
mafi<-paste('mac',mac)
new_<-subset(output,output$MAC>=mac&output$Plot<lower.limit)
if (is.na(thres)==T) {thres<--log10(0.05/nrow(subset(output,output$MAC>=mac)))}}

if(is.na(max.y)==T) {
max.y<-ceiling(-log10(min(new_[,h]))+1)} 

output_<-new_[order(new_$Pos),]
output_ok<-output_[order(output_$Chr),]

maxpos<-c(0,cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x+max(cumsum(aggregate(output_ok$Pos,list(output_ok$Chr),max)$x))/100))
if(black==T) {
plot_col<-rep(c('gray10','gray60'),ceiling(max(unique(output_ok$Chr))/2))
} else { 
require(RColorBrewer)
colorCount = length(unique(output_ok$Chr))
getPalette = colorRampPalette(brewer.pal(colorCount, "Set1"))

plot_col<-getPalette(colorCount)}
size<-aggregate(output_ok$Pos,list(output_ok$Chr),length)$x

a<-rep(maxpos[1],size[1])
b<-rep(plot_col[1],size[1])
for (i in 2:max(unique(output_ok$Chr))){
a<-c(a,rep(maxpos[i],size[i]))
b<-c(b,rep(plot_col[i],size[i]))}

output_ok$xpos<-output_ok$Pos+a
output_ok$col<-b

d<-(aggregate(output_ok$xpos,list(output_ok$Chr),min)$x+aggregate(output_ok$xpos,list(output_ok$Chr),max)$x)/2

plot(output_ok$xpos,-log10(output_ok$Plot),col=output_ok$col,pch=16,ylab='-log10(pval)',xaxt='n',xlab='chromosome',axes=FALSE,cex=1.2,ylim=c(-log10(lower.limit),max.y))
axis(1,tick=FALSE,at=d,labels=c(1:max(unique(output_ok$Chr))))
axis(2,lwd=2)
abline(h=thres,lty=9,col='black',lwd=2)
if(title==TRUE) {
title(main=name)}
if(is.character(show)==T) {
poso<-list()
for (z in 1:length(show)) {
poso[[z]]<-tair9[which(tair9[,3]==show[[z]]),1]+maxpos[tair9[which(tair9[,3]==show[[z]]),4]]}
abline(v=poso,lty=2,col=2)}

}


qq_plot<-function(output,h=8,farbe='red',max.y=NA,maf=0.05,sex=2) {
colnames(output)[h]<-'Plot'
out_<-subset(output,output$MAF>maf)

out<-subset(output,output$MAF>maf&output$Plot<0.05)

e<--log10(ppoints(nrow(out_)))[1:nrow(out)]
o<--log10(sort(out$Plot))
if(is.na(max.y)==T) {
max.y<-ceiling(max(unlist(o)))+1} 


plot(e,o,type='l',col=farbe,xlab=expression(Expected~~-log[10](italic(p))), ylab=expression(Observed~~-log[10](italic(p))),xlim=c(0,max(e)+1),ylim=c(0,max.y),lwd=sex,axes=FALSE)
abline(0,1,col="dark grey")
axis(1,lwd=2)
axis(2,lwd=2)
}




## set parameters
#high=T
local_plot<-
function(output,h=8,chr=0,pos=0,rank=1,high=FALSE,maf=0.05,window_size=50,cex=1,cex.c=2,thres=8){
colnames(output)[h]<-'Plot'
new_<-subset(output,output$MAF>maf)
colnames(new_)[2]<-'chromosomes'
jet.colors = colorRampPalette(c("black", "black", "black", "blue",
"#7FFF7F", "yellow", "#FF7F00", "red")) 

### or provide a specific SNP
if (chr!=0) {
z<-which(new_[order(new_[,h]),2]==chr& new_[order(new_[,h]),3]==pos)
} else {z<-rank}


wind<-window_size*1000
a<-new_[order(new_[,h]),][z,2]
b<-new_[order(new_[,h]),][z,3]-wind
c<-b+2*wind
stopifnot(new_[order(new_[,h]),][z,h]<0.01)
localo <- subset(new_,chromosomes==a&Pos>b&Pos<c&Plot<0.01)

vu<--log10(min(localo[,h]))+2

xx<-X_ok[,colnames(X_ok)%in%localo$SNP]
snp_mean<-apply(xx,2,mean)
snp_sd<-apply(xx,2,sd)
snpm<-matrix(nrow=nrow(xx),ncol=ncol(xx),data=snp_mean,byrow=T)
snpd<-matrix(nrow=nrow(xx),ncol=ncol(xx),data=snp_sd,byrow=T)
xx_stand<-(xx-snpm)/snpd
r2<-(crossprod(xx_stand,xx_stand)/nrow(xx))^2

r2<<-r2/r2[1,1]


if(high==T) {
lead_snp<-as.character(localo[which(localo[,h]==min(localo[,h])),1][1])

} else { lead_snp<-new_[order(new_[,h]),][z,1]}

ld<-r2[rownames(r2)%in%lead_snp,]


ld_col<-jet.colors(21)[round(ld*20)+1]


plot(localo[,3],-log10(localo[,h]),type='p',pch=19,xlab='',ylab='-log(p)',col=ld_col,ylim=c(2,vu),xlim=c(b,c),axes=FALSE,cex=cex)
par(new=T)
plot(localo[which(localo[,1]%in%lead_snp),3],-log10(localo[which(localo[,1]%in%lead_snp),h]),type='p',pch=18,xlab='',ylab='-log(p)',ylim=c(2,vu),cex=cex.c,col='red',xlim=c(b,c),axes=FALSE,main=localo[which(localo[,1]%in%lead_snp),1])
abline(h=thres,lty=2,lwd=2)
#abline(h=-log10(0.05/ncol(X_ok)),lty=2,lwd=2)
axis(1,lwd=2)
axis(2,lwd=2)

}

### plot color scale 
color.scale<-function(x){
jet.colors = colorRampPalette(c("black", "black", "black", "blue",
"#7FFF7F", "yellow", "#FF7F00", "red")) 

#hallo<-c(0,0,0.1,0.1,0.2,0.2,0.3,0.3,0.4,0.4,0.5,0.5,0.6,0.6,0.7,0.7,0.8,0.8,0.9,0.9,1,1)
plot(1,1,col='white',xlim=c(0,2),ylim=c(0,1),xlab='',ylab='',axes=FALSE)
axis(4,lwd=4)

for(z in 1:x) {
 rect(0.9,(0+((z-1)/x)),1.1,(0+((z)/x)),col=jet.colors(x)[z],border='gray')
 par(new=T)}
# if( z %% 2 !=0) {text(c,(+((z-1)/5)),hallo[z])}}
 #text(c,(2+((z)/5)),hallo[21])
}

		
		
