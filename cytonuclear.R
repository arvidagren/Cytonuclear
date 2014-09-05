#Chromosomal distribution of cyto-nuclear genes in a dioecious plant with sex chromosomes
#Genome Biology and Evolution; DOI:  10.1093/gbe/evu197

#Script for determining the representation of cyto-nuclear genes  on X and autosomal chromosomes in R. hastatulus

#import data
#data should be a .csv annotation file with a column called "GO" containing GO ID's for each gene
#for R. hastatulus these files are in the directory called "cytonuclear". 
#change 'yourdirectory' below to the directory containing "cytonuclear" 
goX=read.csv("/cytonuclear/goX.csv")
goA=read.csv("/cytonuclear/goA.csv")
gohemi=read.csv("/cytonuclear/gohemi.csv")

#numbers, proportions, and expected proportions of cyto-nuclear (and non-cyto-nuclear) genes
XN=length(unique(goX$gene));XN
Xmito=length(grep("0005739",goX$GO));Xmito
X_nomito=XN-Xmito;X_nomito
Xchl=length(grep("0009507",goX$GO));Xchl
X_nochl=XN-Xchl;X_nochl

autN=length(unique(goA$gen));autN
autmito=length(grep("0005739",goA$GO));autmito
aut_nomito=autN-autmito;aut_nomito
autchl=length(grep("0009507",goA$GO));autchl
aut_nochl=autN-autchl;aut_nochl

hemiN=length(unique(gohemi$gene));hemiN
hemimito=length(grep("0005739",gohemi$GO));hemimito
hemi_nomito=hemiN-hemimito;hemi_nomito
hemichl=length(grep("0009507",gohemi$GO));hemichl
hemi_nochl=hemiN-hemichl;hemi_nochl

totN=XN+autN+hemiN;totN
mito_tot=Xmito+autmito+hemimito;mito_tot
exp_mito=mito_tot/totN;exp_mito
xexp_mito=exp_mito*XN;xexp_mito
autexp_mito=exp_mito*autN;autexp_mito
hemiexp_mito=exp_mito*hemiN;hemiexp_mito

chl_tot=Xchl+autchl+hemichl;chl_tot
exp_chl=chl_tot/totN;exp_chl
xexp_chl=exp_chl*XN;xexp_chl
autexp_chl=exp_chl*autN;autexp_chl
hemiexp_chl=exp_chl*hemiN;hemiexp_chl

xratio_mito=Xmito/xexp_mito;xratio_mito
xratio_chl=Xchl/xexp_chl;xratio_chl
autratio_mito=autmito/autexp_mito;autratio_mito
autratio_chl=autchl/autexp_chl;autratio_chl
hemiratio_mito=hemimito/hemiexp_mito;hemiratio_mito
hemiratio_chl=hemichl/hemiexp_chl;hemiratio_chl

#Fisher's Exact test for count data: Mito-nuclear
m_mito=matrix(c(Xmito,X_nomito,autmito,aut_nomito,hemimito,hemi_nomito), nrow=2,ncol=3, 
              dimnames=list(c("mito-nuclear","Not mito-nuclear"),
                            c("X-genes","Autosomes","HemiX")));m_mito
fisher.test(m_mito,simulate.p.value=T,B=10000)

#Fisher's Exact test for count data: Chloro-nuclear
m_chl=matrix(c(Xchl,X_nochl,autchl,aut_nochl,hemichl,hemi_nochl),nrow=2,ncol=3, 
             dimnames=list(c("Chloro-nuclear","Not Chloro-nuclear"),
                           c("X-genes","Autosomes","HemiX")));m_chl
fisher.test(m_chl,simulate.p.value=T,B=10000)


#BOOTSTRAPPING & CONFIDENCE INTERVALS

#MITONUCLEAR
#Bootstrapped CI's for mito-X-genes (numbers)
library(boot)
set.seed(123)
b=c(rep(0, X_nomito),rep(1, Xmito));b
fun<-function(b,i){
  sum(b[i]==1)
}
bootprop <- boot(b,fun,10000); bootprop
plot(bootprop)
boot.ci(bootprop, conf=0.95, type="basic")
#X-linked CI's from bootstrapping, divided by expected values (mito-nuclear)
CIlowXm=boot.ci(bootprop,conf=0.95,type="basic")$basic[4]
CIhighXm=boot.ci(bootprop,conf=0.95,type="basic")$basic[5]
CIlowX_ratiom=CIlowXm/xexp_mito;CIlowX_ratiom
CIhighX_ratiom=CIhighXm/xexp_mito;CIhighX_ratiom

#Bootstrapped CI's for Autosomes (mito-nuclear)
set.seed(123)
b2=c(rep(0, aut_nomito),rep(1, autmito));b2
fun2<-function(b2,i){
  sum(b2[i]==1)
}
bootprop2 <- boot(b2,fun2,10000); bootprop2
plot(bootprop2)
boot.ci(bootprop2, conf=0.95, type="basic")
#Autosomal CI's from bootstrapping divided by expected values (chloro-nuclear)
CIlowAm=boot.ci(bootprop2,conf=0.95,type="basic")$basic[4]
CIhighAm=boot.ci(bootprop2,conf=0.95,type="basic")$basic[5]
CIlowA_ratiom=CIlowAm/autexp_mito;CIlowA_ratiom
CIhighA_ratiom=CIhighAm/autexp_mito;CIhighA_ratiom

#Bootstrapped CI's for hemizygous genes (mito-nuclear)
set.seed(123)
b3=c(rep(0, hemi_nomito),rep(1, hemimito));b3
fun3<-function(b3,i){
  sum(b3[i]==1)
}
bootprop3 <- boot(b3,fun3,10000); bootprop3
plot(bootprop3)
boot.ci(bootprop3, conf=0.95, type="basic")
#hemizygous CI's from bootstrapping divided by expected values (chloro-nuclear)
CIlowhemim=boot.ci(bootprop3,conf=0.95,type="basic")$basic[4]
CIhighhemim=boot.ci(bootprop3,conf=0.95,type="basic")$basic[5]
CIlowhemi_ratiom=CIlowhemim/hemiexp_mito;CIlowhemi_ratiom
CIhighhemi_ratiom=CIhighhemim/hemiexp_mito;CIhighhemi_ratiom


#CHLORONUCLEAR
#Bootstrapped CI's for chloro-X-genes (numbers)
library(boot)
set.seed(123)
b4=c(rep(0, X_nochl),rep(1, Xchl));b4
fun4<-function(b4,i){
  sum(b4[i]==1)
}
bootprop4 <- boot(b4,fun4,10000); bootprop4
plot(bootprop4)
boot.ci(bootprop4, conf=0.95, type="basic")
#X-linked CI's from bootstrapping, divided by expected values (chloro-nuclear)
CIlowXC=boot.ci(bootprop4,conf=0.95,type="basic")$basic[4]
CIhighXC=boot.ci(bootprop4,conf=0.95,type="basic")$basic[5]
CIlowX_ratioC=CIlowXC/xexp_chl;CIlowX_ratioC
CIhighX_ratioC=CIhighXC/xexp_chl;CIhighX_ratioC

#Bootstrapped CI's for Autosomes (chloro-nuclear)
set.seed(123)
b5=c(rep(0, aut_nochl),rep(1, autchl));b5
fun5<-function(b5,i){
  sum(b5[i]==1)
}
bootprop5 <- boot(b5,fun5,10000); bootprop5
plot(bootprop5)
boot.ci(bootprop5, conf=0.95, type="basic")
#Autosomal CI's from bootstrapping divided by expected values (chloro-nuclear)
CIlowAC=boot.ci(bootprop5,conf=0.95,type="basic")$basic[4]
CIhighAC=boot.ci(bootprop5,conf=0.95,type="basic")$basic[5]
CIlowA_ratioC=CIlowAC/autexp_chl;CIlowA_ratioC
CIhighA_ratioC=CIhighAC/autexp_chl;CIhighA_ratioC

#Bootstrapped CI's for hemizygous genes (chloro-nuclear)
set.seed(123)
b6=c(rep(0, hemi_nochl),rep(1, hemichl));b6
fun6<-function(b6,i){
  sum(b6[i]==1)
}
bootprop6 <- boot(b6,fun6,10000); bootprop6
plot(bootprop6)
boot.ci(bootprop6, conf=0.95, type="basic")
#hemizygous CI's from bootstrapping divided by expected values (chloro-nuclear)
CIlowhemiC=boot.ci(bootprop6,conf=0.95,type="basic")$basic[4]
CIhighhemiC=boot.ci(bootprop6,conf=0.95,type="basic")$basic[5]
CIlowhemi_ratioC=CIlowhemiC/hemiexp_chl;CIlowhemi_ratioC
CIhighhemi_ratioC=CIhighhemiC/hemiexp_chl;CIhighhemi_ratioC

#plot
par(mfrow=c(1,2))
dotchart(c(xratio_mito,autratio_mito,hemiratio_mito),cex=1,pch=19, col="black",
         labels=c("X chromosome","Autosomes","Hemizygous X"),xlim=c(0,2),
         xlab="observed/expected", main="Mito-nuclear genes")
abline(v=1,col="red",lwd=2,lty=2)
segments(CIlowX_ratiom,1,CIhighX_ratiom,1)
segments(CIlowA_ratiom,2,CIhighA_ratiom,2)
segments(CIlowhemi_ratiom,3,CIhighhemi_ratiom,3)
dotchart(c(xratio_chl,autratio_chl,hemiratio_chl),cex=1,pch=19, col="black",
         labels=c("X chromosome","Autosomes","Hemizygous X"),xlim=c(0,2),
         xlab="observed/expected",main="Chloro-nuclear genes")
abline(v=1,col="red",lwd=2,lty=2)
segments(CIlowX_ratioC,1,CIhighX_ratioC,1)
segments(CIlowA_ratioC,2,CIhighA_ratioC,2)
segments(CIlowhemi_ratioC,3,CIhighhemi_ratioC,3)
