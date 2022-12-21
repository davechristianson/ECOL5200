
# Create function to simlate a population meeting the parameters desired



library(extraDistr)
library(dplyr)
#'make_group_pis() function estimates negative hypergeometric probabilities of breakpoints between sorted groups occurring at each ranked individual in the population.  The negative hypergometric will maintain its mass in within 1 to n+m range with most mass nearest the smallest possible group size (r) but a very long tail (a few very large groups).  Here n should be relatively large to accommodate the largest possible group size *but can be larger to accomodate a different shape to the distribution*. r should be the smallest possible group size (usually 1), and m should be manipuated to achieve the skew desired.  With r = 1 the distribution goes from nearly uniform to increasingly right skewed, as m increases from r to n.  
# This function also produces a histogram of the group size distribution assumed by the settings for n,m, and r. THis allows for a visual assessment of the fit. To create the histogram, you must provide 
# 'gamma' is the expected groupsize. 
# 'lambda' is the expected population size (negatively binomially distributed) at a site. 
# 'size' is the dispersion parameter for the negative binomial . 



sample.int_v<-Vectorize(sample.int)

make_group_pis<-function (plothist=T,
                          n=500, 
                          r=1, 
                          m= 25,
                          nn=10000,
                          gamma=10,
                          lambda=50,
                          size=1) {
  break_pis<-table(factor(rnhyper(nn=nn,n=n,m=m,r=r), levels=seq(1,(n+r))))/nn
  break_pis<-break_pis+(0.5/nn)

  if(plothist){
    N<-rnbinom(nn,size=size,mu=lambda)
    G<-apply(cbind(rpois(n=nn,lambda=N/gamma)+1,N),1,min) # truncation when N is not 0 (at least 1 group when N>0)
    
    break_pis_list<-lapply(N[N>0], function(x) head(break_pis, n=x))
    breaks_list<-sample.int_v(n=N[N>0],size=G[N>0]-1,replace=F,prob=break_pis_list) # there are G-1 breaks in a population separated into G groups
    breaks_list_sort<-lapply(breaks_list,sort)
    gs_list<-lapply(mapply(function(x,y) c(0,x,y),breaks_list_sort,N[N>0]),diff)
    hist(unlist(gs_list))
  }
  else{
  }
  return(break_pis)

}

test<-make_group_pis()


