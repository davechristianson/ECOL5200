rm(list=ls())

library(extraDistr)
library(dplyr)
library(optimx)

# this script simulates group sizes in a way that will allow a mean group size, mean abundance, and overdispersion in abundance to be specificied


# when their is only one instance of each case, sample.int is Wallenius' multinomial noncentral negative hypergeometric
# but we need it vectorized
sample.int_v <- Vectorize(sample.int)


make_group_pis <- function(r=1, nn=10000, tgamma=6, mabundance=30, sizeN=1, igamma=1.2) {
  N <- rnbinom(nn, size=sizeN, mu=mabundance) # realized abundance at nn sites
  G <- apply(cbind(rpois(n=nn, lambda=N/(tgamma*igamma))+1, N), 1, min) # number of groups is poisson draw at rate N/GS but must me no larger than N (number of groups cannot exceed population size)
  break_pis0 <- (table(factor(rnhyper(nn=nn*10, n=max(N), m=mabundance, r=r), levels=seq(1,((max(N))+r))))/(nn*10))+((1/(max(N)+r))/(nn*10))  # bin counts from mega sample of each integer from negative hypergeo.  Also 1/nn*10 individuals added to each bin (to avoid 0 mass in some bins). # this histogram used to estimate probabilites 
  break_pis_list <-list()
  break_pis_list[N>1] <- lapply(N[N>1], function(x) head(break_pis0, n=x-1)) # ragged list of pi's: unique length N for each OCCUPIED site
  breaks_list <-list()
  breaks_list[N>1] <- try(sample.int_v(n=N[N>1]-1, size=G[N>1]-1, replace=F, prob=break_pis_list[N>1]), silent=T) # sample G-1 RANKS (or cumulative sums) from Wallenius' Multinomial noncentral negative hybergeometric using pis from above.  
  breaks_list_sort<-list()
  breaks_list_sort[N>1] <- lapply(breaks_list[N>1], sort) # sort the G-1 RANKS (aka cumulative sums)
  breaks_list_sort_comp<-list()
  breaks_list_sort_comp<-mapply(function(x,y) c(x,y), breaks_list_sort,N) # Add N as the rank of the final group (total cumulative sum)
  breaks_list_sort_comp<-lapply(breaks_list_sort_comp,function(x)replace(x,sum(x)==0,NA))
  gs_list <- lapply(mapply(function(x,y) c(0,x), breaks_list_sort_comp), diff) # extract each groups size as the difference between RANKED members
  gs_list <- lapply(gs_list,function(x)replace(x,sum(x)==0,NA)) # group sizes of 0 can throw off estimates
  
  mgs <- mean(unlist(gs_list),na.rm=T) # means group size
  return(list(Rank_Pis_list=break_pis0, MGS=mgs, muGS=tgamma,GS_list=gs_list,Rank_list=breaks_list_sort_comp,N=N))
}

make_group_pis()


# Change tgamma and mabundance to be your sample gs and your expected mean abundance.
# optimize the dispersion of the population abundances and the adjustment to tgamma

# define mean group size and mean site abundance
tgamma=3
mabundance=30

optim_fun <- function(x) {
  # Extract the values of igamma and size from the input vector
  igamma <- x[1]
  sizeN <- x[2]
  tgamma <-tgamma
  mabundance<-mabundance

  # Calculate the mean group size using the make_group_pis function
  result <- make_group_pis(igamma=igamma, sizeN=sizeN,tgamma=tgamma,mabundance=mabundance)
  mgs <- result[[2]]
  tgamma <- result[[3]]
  
  # Return the difference between tgamma and mgs as the objective function to be minimized
  return(abs(tgamma - mgs))
}

# Set the lower and upper bounds for the optimization
lower_bounds <- c(0.5, 0.1)
upper_bounds <- c(3, mabundance)

# Set the initial values for the optimization
initial_values <- c(1,1)

# Run the optimization
optim_result <- optimr(par=initial_values, fn=optim_fun,  lower=lower_bounds, upper=upper_bounds, method="L-BFGS-B")
optim_result

#run again with new initial values
initial_values <- c(optim_result$par[[1]],optim_result$par[[2]])
optim_result <- optimr(par=initial_values, fn=optim_fun,  lower=lower_bounds, upper=upper_bounds, method="L-BFGS-B")
optim_result


# Extract the optimized values of m, igamma, and size
optimized_igamma <- optim_result$par[[1]]
optimized_sizeN <- optim_result$par[[2]]

# Calculate the mean group size using the make_group_pis function with the optimized values

my_pis<-make_group_pis(igamma=optimized_igamma, sizeN=optimized_sizeN,tgamma=tgamma,mabundance=mabundance)

# is the ration of mean of the group size distribution to the mean of the rank distribution the same as the ration of the meangs 
sum(my_pis$N)
mean(my_pis$N)
my_pis$muGS

my_pis$MGS

my_pis$MGS/mean(my_pis$N)

mean_typ_GS<-(sum(unlist(my_pis$GS_list)*unlist(my_pis$GS_list),na.rm=T)/sum(unlist(my_pis$GS_list),na.rm=T))


mean(unlist(my_pis$GS_list),na.rm=T)
var(unlist(my_pis$GS_list),na.rm=T)

mean(unlist(my_pis$Rank_list),na.rm=T)
var(unlist(my_pis$Rank_list),na.rm=T)

plot(lapply(my_pis$Rank_list,length),
lapply(my_pis$Rank_list,mean,na.rm=T))
abline(lm(unlist(lapply(my_pis$Rank_list,mean,na.rm=T))~unlist(lapply(my_pis$Rank_list,length))))
lm(unlist(lapply(my_pis$Rank_list,mean,na.rm=T))~unlist(lapply(my_pis$Rank_list,length)))

mean(unlist(my_pis$GS_list),na.rm=T)/mean(unlist(my_pis$Rank_list),na.rm=T)

mean(unlist(lapply(my_pis$GS_list,length))[unlist(lapply(my_pis$GS_list,length))!=0])




dev.off()
par(mfrow=c(2,1))
hist(unlist(my_pis$GS_list),xlab="Group sizes",breaks=seq(0,max(unlist(my_pis$Rank_list),na.rm=T)+2,by=2))
hist(unlist(my_pis$Rank_list),xlab="Ranks",breaks=seq(0,max(unlist(my_pis$Rank_list),na.rm=T)+2,by=2))



plot(unlist(my_pis$GS_list)~unlist(my_pis$Rank_list))
lm(log(unlist(my_pis$GS_list))~log(unlist(my_pis$Rank_list)))


table(unlist(my_pis$GS_list))
table(unlist(my_pis$Rank_list))


# weighting the rank score of the last individual in the group by the group size
wtd_ranks<-Map(function(x, y) { x * y }, my_pis[[5]], my_pis[[4]])

# weighting the group size of every group by the size of the group (typical group)
wtd_groupsizes<-Map(function(x, y) { (x * y) /sum(x) }, my_pis$GS_list, my_pis$GS_list)

#typical group size
typ_groupsizes<-Map(function(x) { sum(x * x) /sum(x) }, my_pis$GS_list)

dev.off()
par(mfrow=c(2,1))
hist(unlist(lapply(my_pis$GS_list,mean,na.rm=T)),breaks=seq(0,max(unlist(my_pis$GS_list),na.rm=T)+2,by=2))
hist(unlist(typ_groupsizes),breaks=seq(0,max(unlist(my_pis$GS_list),na.rm=T)+2,by=2))

hist(unlist(typ_groupsizes))




sum(unlist(typ_groupsizes),na.rm=T)

# relationship between the rank of individuals and their group size

allranks<-unlist(my_pis[[5]])
allgroupsizes<-unlist(my_pis[[4]])
allgroupsizes2<-allgroupsizes^2
plot(allranks~allgroupsizes)
abline(lm(allranks~allgroupsizes))

# the ranks of ramked individuals weighted by the group size
cumul_wtd_ranks<-tapply(unlist(wtd_ranks),factor(unlist(my_pis[[5]]),levels=1:max(unlist(my_pis[[5]]))),sum,na.rm=T)

# the group size of individuals weighted by their gorup size (crowding)
cumul_wtd_gs<-tapply(unlist(wtd_groupsizes),factor(unlist(my_pis[[4]]),levels=1:max(unlist(my_pis[[4]]))),sum,na.rm=T)

# distribution of ranks, wtd by gorup sizes
plot(cumul_wtd_ranks~c(1:89))

#distribution of group sizeswetd by indivdiuals (crowding)
plot(c(cumul_wtd_gs/sum(cumul_wtd_gs,na.rm=T))~c(1:75))

# distibution of group sizes across population vs distribution of indivdiuals across group sizes
hist(allgroupsizes,freq=F,xlab="group size")
points(1:75,c(cumul_wtd_gs/sum(cumul_wtd_gs,na.rm=T)))

# distribution of ranks across population of groups vs distribution of individuals across ranks
hist(allranks,freq=F,xlab="rank")
points(1:89,c(cumul_wtd_ranks/sum(cumul_wtd_ranks,na.rm=T)))

# take-home
# when populations are negatively binomially distributed, group sizes follow a negative hypergeometric
# distribution.  N, gamma, and G are all interdependent.  THere is information about N in the distrubiotn of group sizes. Large rare groups infer large large rare populaiton sizes.  THe existence of rare large groups is evidence for the existence of rare large populatoin sizes. missing a single large gorup can greatly effect abundanc models.  There is little penalty or value in identifyig small gorups that hold a small fraction of the population.  The distribution of indivdiuals is more hypergeometric or binomially distrubted in the population over the range of possible group sizes.  rank of individuals is correalted with group size.  Group size observations hold information about the distribution of the population across groups.  This could be useful to infer detection probability of individuals in the population rather than gorups.  next steps. 




