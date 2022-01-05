## IICR with K classes along the genome, each following a symmetrical n-island model with stationary (i.e. constant over time) deme size.
## Used in Figure 4
## IICR is computed analytically from Equations (5-6) in Boitard et al (2022).

library(tidyverse)
library(latex2exp)

# generic function to compute R (numerator of the IICR) in a given class
Rstruct.nislandc <- function (t,c,M,n) {
	# t: time where the R is computed (scalar)
	# c: coalescence rate (inverse of local deme size) in each class
	# M: scaled migration rate
	# n: number of islands (or demes)
	gamma <- M/(n-1)
	delta <- (c+n*gamma)^2-4*c*gamma
	alpha <- (c+n*gamma+sqrt(delta))/2
	beta <- (c+n*gamma-sqrt(delta))/2
	R <- c*(1/alpha)*(gamma-alpha)/(beta-alpha)*exp(-alpha*t)+c*(1-(gamma-alpha)/(beta-alpha))*(1/beta)*exp(-beta*t)
	return(R)
}

# generic function to compute f (denominator of the IICR) in a given class
# parameters as in Rstruct.nislandc
fstruct.nislandc <- function(t,c,M,n) {
	gamma <- M/(n-1)
	delta <- (c+n*gamma)^2-4*c*gamma
	alpha <- (c+n*gamma+sqrt(delta))/2
	beta <- (c+n*gamma-sqrt(delta))/2
	f <- c*(gamma-alpha)/(beta-alpha)*exp(-alpha*t)+c*(1-(gamma-alpha)/(beta-alpha))*exp(-beta*t)
	return(f)
}

## Figure 4 - K=2 classes
n=10 # number of islands (or demes)
N=1000 # number of haploids in one deme
v_t=seq(from=0,to=10*n,by=0.01) # time vector (when the IICR will be computed)
lambda=c(0.1,1) # relative Ne (deme size) in each class
c=1/lambda
v_a2=c(0,0.01,0.1,0.5,0.9,0.99,1) # possible values for the proportion of class 2
v_M=c(0.5,1,5) # possible values for M
p=length(v_a2)*length(v_t)*length(v_M)
# builds the table of IICR values
res=data.frame(matrix(nrow=p,ncol=4))
colnames(res)=c('a2','M','t','IICR')
i=1
for (i_a2 in 1:length(v_a2)){
	for (i_M in 1:length(v_M)){
		for (i_t in 1:length(v_t)){
			res$M[i]=v_M[i_M]
			res$a2[i]=v_a2[i_a2]
			res$t[i]=v_t[i_t]
			a=c(1-res$a2[i],res$a2[i]) # proportion of each class
			num <- a[1]*Rstruct.nislandc(res$t[i],c[1],res$M[i],n) + a[2]*Rstruct.nislandc(res$t[i],c[2],res$M[i],n)
			denom <- a[1]*fstruct.nislandc(res$t[i],c[1],res$M[i],n) + a[2]*fstruct.nislandc(res$t[i],c[2],res$M[i],n)
			res$IICR[i] <- num/denom
			i=i+1
		}
	}
}

# add the IICR under panmixia
IICR.pan.varN <- function (time, lambda.freq, lambda) {
	R <- sum(lambda.freq*exp(-time/lambda))
   	f <- sum( (lambda.freq/lambda) * exp(-time/lambda) )
   	return(R/f)
}

lambda=n*lambda # meta-population size
p=length(v_a2)*length(v_t)
res2=data.frame(matrix(nrow=p,ncol=4))
colnames(res2)=c('a2','M','t','IICR')
i=1
for (i_a2 in 1:length(v_a2)){
	for (i_t in 1:length(v_t)){
		res2$a2[i]=v_a2[i_a2]
		res2$t[i]=v_t[i_t]
		res2$M[i]=1000
		a=c(1-res2$a2[i],res2$a2[i]) # proportion of each class
		res2$IICR[i]=IICR.pan.varN(res2$t[i],a,lambda)
		i=i+1
	}
}
res=rbind(res,res2)

# rescale t with respect to total meta-population size
res$t=res$t/n
res=res%>% filter(t==0 | t>=0.01)

# change panel names
M.labs=c(paste("M=",v_M,sep=''),'panmixia')
names(M.labs)=c(v_M,1000)
# plots in log scale
p=ggplot(res,aes(x=t,y=IICR,color=as.factor(a2)))+geom_line()+theme_bw()+theme(legend.position='bottom')+ylim(0,max(res$IICR)) +scale_color_discrete(name = TeX("$a_2$")) +facet_grid(~M,labeller=labeller(M =M.labs))+xlim(0,2)+scale_x_log10()
ggsave('Figure4.pdf',plot = p, width = 10, height = 4)



