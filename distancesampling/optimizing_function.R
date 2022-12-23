rm(list=ls())

library(extraDistr)
library(dplyr)

sample.int_v <- Vectorize(sample.int)


make_group_pis <- function(r=1, nn=10000, tgamma=6, mabundance=30, size=3, igamma=1.3) {
  N <- rnbinom(nn, size=size, mu=mabundance) # realized abundance at nn sites
  G <- apply(cbind(rpois(n=nn, lambda=N/(tgamma*igamma))+1, N), 1, min) 
  break_pis0 <- (table(factor(rnhyper(nn=nn*10, n=max(N), m=mabundance, r=r), levels=seq(1,((max(N))+r))))/(nn*10))+((1/(max(N)+r))/(nn*10))
  break_pis_list <- lapply(N[N>0], function(x) head(break_pis0, n=x)) # ragged list of pi's
  breaks_list <- try(sample.int_v(n=N[N>0], size=G[N>0]-1, replace=F, prob=break_pis_list), silent=T) 
  breaks_list_sort <- lapply(breaks_list, sort)
  breaks_list_sort_comp<-lapply(seq_along(breaks_list_sort), function(i) c(breaks_list_sort[[i]], N[i]))
  gs_list <- lapply(mapply(function(x,y) c(0,x,y), breaks_list_sort, N[N>0]), diff)
  mgs <- mean(unlist(gs_list))
  return(list(break_pis0, mgs, tgamma,gs_list,breaks_list_sort_comp,N))
}

make_group_pis()



# Change tgamma and mabundance to be your sample gs and your expected mean abundance.
optim_fun <- function(x) {
  # Extract the values of m, igamma, and size from the input vector
  #m <- x[1]
  igamma <- x[1]
  size <- x[2]

  # Calculate the mean group size using the make_group_pis function
  result <- make_group_pis(igamma=igamma, size=size,tgamma=6,mabundance=30)
  mgs <- result[[2]]
  tgamma <- result[[3]]
  
  # Return the difference between tgamma and mgs as the objective function to be minimized
  return(abs(tgamma - mgs))
}

# Set the lower and upper bounds for the optimization
lower_bounds <- c(0.3, 1)
upper_bounds <- c(3, 100)

# Set the initial values for the optimization
initial_values <- c(1,1)

# Run the optimization
optim_result <- optim(par=initial_values, fn=optim_fun,  lower=lower_bounds, upper=upper_bounds, method="L-BFGS-B")

# Extract the optimized values of m, igamma, and size
optimized_igamma <- optim_result$par[1]
optimized_size <- optim_result$par[2]

# Calculate the mean group size using the make_group_pis function with the optimized values
# specify the mean group size observed in the sample
mean_gs <- 6
# specify the expexted mean abundance at a site, this can be a guess
mean_abundance <- 30

# now fit with optimized parameters

my_pis<-make_group_pis(igamma=optim_result$par[1], size=optim_result$par[2],tgamma=6,mabundance=30)

# weighting the rank score of the last individual in the group by the group size
wtd_ranks<-Map(function(x, y) { x * y }, my_pis[[5]], my_pis[[4]])

# weighting the group size of every group by the size of the group (crowding)
wtd_groupsizes<-Map(function(x, y) { x * y }, my_pis[[4]], my_pis[[4]])

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




