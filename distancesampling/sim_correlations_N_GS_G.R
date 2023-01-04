# functions and plotting to demosntrate relationships between group size and density under several simulated populaiton sizes



#Compare populations of 1 10 and 100 individuals
#Group sizes of mean of 1, 2 and 10 individuals
# compare model 1:constant gs, G must increase with N, 2: constat G, gs must increas with N, 3: inverse effects of site covaraiate on N (positive) and gs (negative)
rm(list=ls())
set.seed(1234567)
# set mean abundance 

# function for extraxting glm coefficent for effect of N on group size
GS_N_correlation<-function(mean.lam=1, # mean site abundance
                           beta1=0.5, # effects of spatial covariate ('habitat') on log mean abundance
                           nsites=100, # number of spatial replicates
                           sizeN=5, # dispersion parameter in negative binomial abundance
                           gamma0=log(1), # log mean group size (-1)
                           gamma1=0, # effect of spatial covariate ('habitat') on log mean group size 
                           mu_G=NA) # fixed mean group abundance (ignores group size model) 
  {
  habitat<-rnorm(nsites)
  beta0<-log(mean.lam)
  lambda<-exp(beta0 + beta1*habitat)
  N<-rnbinom(nsites,size=sizeN, mu=lambda)
  
  if(is.na(mu_G)){
    mu_GS_tr<-exp(gamma0 + gamma1*habitat)
    G<-apply(cbind(N,rpois(nsites, N/(mu_GS_tr))+1),1,min) 
  }else{
    G<-apply(cbind(N,rpois(nsites, mu_G)+1),1,min)
  }
  real_mu_site_GS<-N[N!=0]/G[N!=0]
  mean_real_mu_site_GS<-mean(real_mu_site_GS,na.rm=T)
  #coef(glm(real_mu_site_GS~N[N!=0],family="quasipoisson"))[2]
  GS_N_cor<-(cor(real_mu_site_GS,N[N!=0]))
  data<-list(GS_N_cor=GS_N_cor,G=G,N=N, sitemean_GS=real_mu_site_GS,meanGS=mean_real_mu_site_GS,meanG=mean(G[N!=0],na.rm=T))
}

mydata<-GS_N_correlation()

nsims<-1000
GS_N_cor_arr<-array(NA,dim=c(nsims,3,3,3))
GS_N_meanG_arr<-array(NA,dim=c(nsims,3,3,3))
testlam<-c(1,10,100) # mean individual abundances to evaluate
testmuGS<-log(c(1,2,10)) # mean group sizes -1  to evaluate
testmodel<-c(1,2,3) # number of model formulations to test
testG<-c(NA,5,NA) # evaluate fixed mean G (mean Group Abundance) or G allowed to vary (NULL)?
testgamma1<-c(0,0,-0.5) # evaluate an inverse effect of spatial covarate in group size

for (i in 1:nsims){
  for(n in 1:length(testlam)){
    for(gs in 1:length(testmuGS)){
      for(m in 1:length(testmodel)){
        GS_N_cor_arr[i,n,gs,m]<-GS_N_correlation(mean.lam=testlam[n],nsites=100,gamma0=testmuGS[gs],mu_G=testG[m],gamma1=testgamma1[m])[[1]]
        GS_N_meanG_arr[i,n,gs,m]<-GS_N_correlation(mean.lam=testlam[n],nsites=100,gamma0=testmuGS[gs],mu_G=testG[m],gamma1=testgamma1[m])[[6]]
      }
    }
  }
}
# result is an array where the rows are sims, the columns are abundances, the matrices are first mean group sizes, then models


# now plot
# negative correlation is betwewn GS and N is highly unlikely
dev.off()
par(mar=c(10,10,1,1))
plot(c(4,20),c(-1,1),type="n",xlab="model",xaxt="n",yaxt="n",ylab="")
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
for(m in 1:length(testmodel)){
  for(n in 1:length(testlam)){
    curcol<-c("red","black","blue")[n]
    for(gs in 1:length(testmuGS)){
      curpt<-c(15,16,17)[gs]
      x<-c(-.4,0,.4)[n]+c(-1.5,0,1.5)[gs]+(m*6)
      nna<-sum(is.na(GS_N_cor_arr[1:nsims,n,gs,m]))
      segments(x,sort(GS_N_cor_arr[1:nsims,n,gs,m])[floor(0.025*(nsims-nna))],x,sort(GS_N_cor_arr[1:nsims,n,gs,m])[floor(0.975*(nsims-nna))],col=curcol,lwd=1)
      points(x,mean(GS_N_cor_arr[1:nsims,n,gs,m],na.rm=T),col="black",cex=1.5,pch=curpt,bg='white')
    }
  }
}
legend(5,-0.5,legend=c("1","10","100"),lwd=1,lty=1,col=c("red","green","blue"),title="mean\nsite\nabundance",bty="n")
legend(10,-.7,legend=c("1","2","10"),pt.cex=2,pch=c(15,16,17),title="mean\ngroup\nsize",bty="n",ncol=3)
axis(1,at=c(6,12,18),labels=c("fixed GS","fixed G","inverse"))
axis(2,at=c(-1,-0.5,0,0.5,1),las=1)
mtext("correlation\nbetween mean\ngroup size\nand N across\n100 sites",2,las=1,5,adj=0.5)
abline(a=0,b=0)
abline(v=c(9,15))



