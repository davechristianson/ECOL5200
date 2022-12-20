
# Create function to simlate a population meeting the parameters desired



library(extraDistr)

#estimate probabilities of breakpoints between groups occuring at ranked individuals in the populatoin.  The negative hypergometric will maintain its mass in within 1 to n+m range with most mass nearest the smallest possible group size (r) but a very long tail (a few very large groups).  Here n should be relatively large to accomodate the largest possible group size *but can be larger*. r should be the smallest possible group size (usually 1), and m should be manipuated to achieve the skew desired.  With r = 1 the disribution goes from nearly uniform to increasingly right skewed, mased on r as m increases from r to n.  
# For visual assessment of the fit, 
# 'gamma' is the expected groupsize. 
# 'lambda' is the expected population size (negatively binomially distrubted) at a site. 
# 'size' is the dispersion parameter for the negative bimiall . 
# this function will display a histogram of the resulting group size distributions you 
# nn needs to be very high to avoid non-zero probs
library(dplyr)

sample.int_v<-Vectorize(sample.int)

make_group_pis<-function (plothist=T,
                          n=200, 
                          r=1 , 
                          m= 50,
                          nn=10000,
                          gamma=5,
                          lambda=50,
                          size=1) {
  break_pis<-table(factor(rnhyper(nn=nn,n=n,m=m,r=r), levels=seq(1,(n+r))))/nn
  break_pis<-break_pis+(0.5/nn)

  if(plothist){
    N<-apply(cbind(rnbinom(nn,size=size,mu=lambda),n+r),1,min)
    G<-apply(cbind(rpois(n=nn,lambda=N/gamma)+1,N),1,min) # truncation when N is not 0 (at least 1 group when N>0)
    
    break_pis_list<-lapply(N[N>0], function(x) head(my_pis, n=x))
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


