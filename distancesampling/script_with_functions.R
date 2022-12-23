
# Create function to simlate a population meeting the parameters desired



library(extraDistr)
library(dplyr)
#'make_group_pis() function estimates negative hypergeometric probabilities of breakpoints between sorted groups occurring at each ranked individual in the population.  The negative hypergometric will maintain its mass in within 1 to n+m range with most mass nearest the smallest possible group size (r) but a very long tail (a few very large groups).  Here n should be relatively large to accommodate the largest possible group size *but can be larger to accomodate a different shape to the distribution*. r should be the smallest possible group size (usually 1), and m should be manipuated to achieve the skew desired.  With r = 1 the distribution goes from nearly uniform to increasingly right skewed, as m increases from r to n.  
# This function also produces a histogram of the group size distribution assumed by the settings for n,m, and r. THis allows for a visual assessment of the fit. To create the histogram, you must provide 
# 'gamma' is the expected groupsize. 
# 'lambda' is the expected population size (negatively binomially distributed) at a site. 
# 'size' is the dispersion parameter for the negative binomial . 
# this function fails multiple times and you may need to keep trying it or reduce nn to get it to work.
rm(list=ls())

sample.int_v<-Vectorize(sample.int)

library(extraDistr)
library(dplyr)
#


make_group_pis<-function (plothist=T, # inpecition of the groupsize distribution resulting from selection of parameters below. this likes to fail so may want to make F.
                          r=1, 
                          nn=100000, # number of random draws from the negative hypergeometric from which to calculate probabilities.
                          NN=10000, # how many simulated populations?
                          tgamma=6,  # idealized mean group size in an unrestricted populatoin - 1. This is used to constrain the expected number of groups at a site based on population size. however number of groups is a poisson random variable at rate N/(gamma+1). When gammma is near N, the realized mean group size will be much smaller because sites with much < gamma individuals will still count as at least 1 group (decreasing mean) and because the right tail is much longer than left tail in Poissons centered near 1. 
                          lambda=25, # for fitting purposes, mean abundance at a site
                          size=5) # for fitting purposes, ind. abundance dispersion parameter in neg binomial, the equivalent of r in JAGS parameterization.
  {
  N<-rnbinom(NN,size=size,mu=lambda) # realized abundance at NN sites
  break_pis0<-table(factor(rnhyper(nn=nn,m=max(N),n=max(N)*tgamma,r=r), levels=seq(1,(max(N)+r))))/nn
  
  if(plothist){

  G<-apply(cbind(rpois(n=NN,lambda=N/(tgamma))+1,N),1,min) # number of groups <= N and no groups when N = 0.
  break_pis_list<-lapply(N[N>0], function(x) head(break_pis0, n=x)) # ragged list of pi's
  breaks_list<-try(sample.int_v(n=N[N>0],size=G[N>0]-1,replace=F,prob=break_pis_list),silent=T) # there are G-1 breaks in a population separated into G groups
  breaks_list_sort<-lapply(breaks_list,sort)
  gs_list<-lapply(mapply(function(x,y) c(0,x,y),breaks_list_sort,N[N>0]),diff)
  par(fig = c(0,1,0,1))  
  hist(unlist(gs_list),breaks=seq(0,max(unlist(gs_list),na.rm=T),by=1),main="Distribution of groups sizes with given N and idealized gamma",xlab="Group size")
  figdims<-par("usr")
  text(0.25*figdims[2],0.9*figdims[4],paste0("mean group size: ",round(mean(unlist(gs_list)[unlist(gs_list)!=0],na.rm=T),1)))
  text(0.25*figdims[2],0.85*figdims[4],paste0("max group size: ",round(max(unlist(gs_list)[unlist(gs_list)!=0],na.rm=T),1)))
  text(0.25*figdims[2],0.8*figdims[4],paste0("var group size: ",round(var(unlist(gs_list)[unlist(gs_list)!=0],na.rm=T),1)))
  text(0.25*figdims[2],0.75*figdims[4],paste0("mean pop size: ",round(mean(unlist(N),na.rm=T),0)))
  text(0.25*figdims[2],0.7*figdims[4],paste0("max pop size: ",round(max(unlist(N),na.rm=T),0)))
  text(0.25*figdims[2],0.65*figdims[4],paste0("var pop size: ",round(var(unlist(N),na.rm=T),0)))
  text(0.25*figdims[2],0.6*figdims[4],paste0("mean # of groups: ",round(mean(unlist(G),na.rm=T),1)))
  text(0.25*figdims[2],0.55*figdims[4],paste0("mean # of groups (occupied sites): ",round(mean(unlist(G)[N>0],na.rm=T),1)))
  text(0.25*figdims[2],0.5*figdims[4],paste0("var # of groups",round(var(unlist(G),na.rm=T),1)))
  par(fig = c(0.5,1, 0.5, 0.9), new = T)  
  hist(unlist(N),main=paste0("Distr. of ", NN," pop. sizes"),xlab='N') 
  par(fig = c(0.5,1, 0.1, 0.5), new = T)  
  hist(unlist(G),main=paste0("Distr. of ", NN," group abund."),breaks=seq(0,max(unlist(G))),xlab='G') 
  }
  return(break_pis0)
}



# likes to fail when visualizing group size distribution. try 100 times, reducing NN by 5% each time, before giving up or set plothist=F (which always works).
dev.off()
my_pis<-NULL
attempt<-0
NN<-10526
while( is.null(my_pis) && attempt <= 100 ) {
  attempt <- attempt + 1
  NN<-ceiling(0.95*NN)
  print(paste0("attempt ",attempt,": NN = ", NN))
  try(
    my_pis<- make_group_pis(NN=NN,lambda=25,tgamma=10),silent=T
  )
} 
hist(rnbinom(100000,size=1000,mu=30))
hist(rnhyper(100000,1000,10,1))

