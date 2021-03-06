## IICR with K classes along the genome, each following a panmictic model whose population size can vary according to a piecewise constant process. 
## Used in Figures 3 & 6
## IICR is computed analytically from Equation (5) in Boitard et al (2022), combined with standard derivations of the IICR under a panmictic model with piecewise constant population size.  

library(tidyverse)
library(latex2exp)
library(gridExtra)

# generic function to compute R (numerator of the IICR) in a given class
coal.cum=function(t,tmodif,lambda){
	# lambda: population size in epoch 0,1,..,E (vector of size E+1)
	# tmodif: times where population size changes (vector of size E)
	# t: time where R is evaluated (scalar)
	ind=which(t>tmodif)
	if (length(ind)==0){# t included in the first period
		R=exp(-t/lambda[1])
	} else {
		iMax=max(ind)
		K=length(tmodif)
		d=tmodif-c(0,tmodif[1:(K-1)])
		R=prod(exp(-d[1:iMax]/lambda[1:iMax]))*exp(-(t-tmodif[iMax])/lambda[iMax+1])
	}
	return(R)
}

# generic function to compute d (numerator of the IICR) in a given class
# parameters as in coal.cum
coal.dens=function(t,tmodif,lambda){
	ind=which(t>tmodif)
	if (length(ind)==0){
		R=exp(-t/lambda[1])/lambda[1]
	} else {
		iMax=max(ind)
		K=length(tmodif)
		d=tmodif-c(0,tmodif[1:(K-1)])
		R=prod(exp(-d[1:iMax]/lambda[1:iMax]))*exp(-(t-tmodif[iMax])/lambda[iMax+1])
		R=R/lambda[iMax+1]
	}
	return(R)
}

## Test: IICR in a panmicitc model with stepwise population size and K=1
lambda=c(3,1,6,10,2)
tmodif=c(1,3,4,5)
t=seq(from=0,to=12,by=0.01)
IICR=rep(0,length(t))
for (i in 1:length(t)){
	IICR[i]=coal.cum(t[i],tmodif,lambda)/coal.dens(t[i],tmodif,lambda)
}
plot(t,IICR)

## Figure 3: IICR in a panmictic model with population size change and K=2
lambda=c(0.1,1) # relative Ne in each class
v_tmodif=c(0.1,0.5,1) # time where population size changes
v_t=seq(from=0,to=10,by=0.01) # time vector (when the IICR will be computed)
v_a2=c(0,0.01,0.1,0.5,0.9,0.99,1) # possible values for the proportion of class 2
 g=5 
# g=100
n=length(v_a2)*length(v_t)*length(v_tmodif)
# builds the table of IICR values
res=data.frame(matrix(nrow=n,ncol=4))
colnames(res)=c('a2','tmodif','t','IICR')
i=1
for (i_a2 in 1:length(v_a2)){
	for (i_g in 1:length(v_tmodif)){
		for (i_t in 1:length(v_t)){
			res$a2[i]=v_a2[i_a2]
			res$t[i]=v_t[i_t]
			res$tmodif[i]=v_tmodif[i_g]
			a=c(1-res$a2[i],res$a2[i]) # proportion of each class
			R=rep(0,2)
			dens=rep(0,2)
			for (j in 1:2){
				R[j]=coal.cum(v_t[i_t],v_tmodif[i_g],lambda[j]*c(g,1))
				dens[j]=coal.dens(v_t[i_t],v_tmodif[i_g],lambda[j]*c(g,1))
			}
			res$IICR[i]=sum(a*R)/sum(a*dens)
			i=i+1
		}
	}
}
# change panel names
#tmodif.labs=paste(TeX("$t_1$"),"=",v_tmodif,sep='')
tmodif.labs=paste0("T=",v_tmodif)
names(tmodif.labs)=v_tmodif
# plots in log scale
p=ggplot(res,aes(x=t,y=IICR,color=as.factor(a2)))+geom_line()+theme_bw()+ylim(0,lambda[2]*g)+xlim(0,10)+scale_x_log10()+scale_color_discrete(name = TeX("$a_2$"))+facet_grid(~tmodif,labeller=labeller(tmodif=tmodif.labs))
ggsave('Figure3.pdf',plot = p, width = 10, height = 4)
#ggsave('FigureS5.pdf',plot = p, width = 10, height = 4) 

## Figure 6: selective sweep in a panmictic population - transient reduction of Ne around the selected site

# function to define a model with several classes of non stationary Ne along the genome, allowing to approximate a given selective sweep scenario. Based on derivations in section 3.2 of Boitard et al (under review).
approx.sweeps=function(N,alpha,r,nL){
	# N diploid population size
	# alpha=2Ns, s selective advantage of homozygote mutants
	# r per site recombination rate
	# nL number of classes within a sweep region
	# returns:
	#	-a: size of each class (vector of size nL)
	# 	-lambda: expected decrease of Ne during the sweep for each class (vector of size nL)
	#	-times=rep(tfix,nL): duration of the sweep (vector of size nL)

	# distance such that only 5% heterozygosity is lost at the end of the sweep
	Lmax=-log(0.05)/log(alpha)*(alpha/(8*N))/r # p0=1/alpha
	dL=Lmax/nL
	vL=((1:nL)-0.5)*dL
	a=rep(dL,nL) # only accounts for one side of the selected site
	# sweep phase - assumes that initial frequency is p0=1/alpha
	tfix=2*N*4*log(alpha)/alpha # fixation time (duration of the sweep)
	q=1-exp(-2*vL*r*log(alpha)*2*N/alpha) # probability to escape the sweep
	times=rep(tfix,nL)
	rcoal=(1-q)*(1-q)/tfix+q*q/(2*N)
	lambda=1/(2*N*rcoal)
	return(rbind(a,times,lambda))
}

# application to the sweep scenario of Schrider et al (2016)
N=10000 # diploid size
v_t0=c(0.2) # end of the sweep (fixation time), time before present in 2N units
v_alpha=c(200,1000,10000) # possible values for the alpha=2Ns
r=10**(-8) # per site recombination rate
nL=10 # number of classes within the sweep region
# builds the table of IICR values
v_t=c(seq(from=0,to=1000/(2*N),by=0.0001),seq(from=1000/(2*N)+0.001,to=100000/(2*N),by=0.001)) # time vector (when the IICR will be computed)
n=length(v_t0)*length(v_alpha)*length(v_t)
res=data.frame(matrix(nrow=n,ncol=4))
colnames(res)=c('t0','alpha','t','IICR')
i=1
for (i_t0 in 1:length(v_t0)){
	for (i_alpha in 1:length(v_alpha)){
		params=approx.sweeps(N,v_alpha[i_alpha],r,nL)
		a=params[1,]*2/(15*10**6) # proportion of each sweep class
		times=params[2,1]/(2*N)+v_t0[i_t0]
		lambda=params[3,]
		# define points where Ne changes
		tmodif=c(v_t0[i_t0],times)
		lambda=rbind(rep(1,nL),lambda,rep(1,nL))
		# add the neutral class
		a=c(a,1-sum(a))
		lambda=cbind(lambda,rep(1,3))
		# compute IICR	
		for (i_t in 1:length(v_t)){
			res$t0[i]=v_t0[i_t0]
			res$alpha[i]=v_alpha[i_alpha]			
			res$t[i]=v_t[i_t]
			R=rep(0,nL+1)
			dens=rep(0,nL+1)
			for (j in 1:(nL+1)){
				R[j]=coal.cum(v_t[i_t],tmodif,lambda[,j])
				dens[j]=coal.dens(v_t[i_t],tmodif,lambda[,j])
			}
			res$IICR[i]=sum(a*R)/sum(a*dens)
			i=i+1
		}
	}
}
# rescale t in generations
res$t=res$t*2*N
# shift curves to improve visibilty
ind=which(res$alpha==1000)
res$IICR[ind]=res$IICR[ind]-0.01
ind=which(res$alpha==10000)
res$IICR[ind]=res$IICR[ind]-0.02
# plots in log scale
p0=ggplot(res,aes(x=t,y=IICR,color=as.factor(alpha)))+geom_line()+theme_bw()+scale_x_continuous(limits=c(1000,100000),trans='log10')+scale_color_discrete(name = TeX("$\\alpha$"))+ylim(0,2)+xlab('generations')
# combines with IICRs obtained from simulated coalescence times (under selection)
source('iicr_selection_loop.R')
gs=list(p0,p1,p2)
lay=rbind(c(NA,1,1,NA),c(2,2,3,3))
select_grobs <- function(lay) {
	id <- unique(c(t(lay))) 
  	id[!is.na(id)]
} 
p=grid.arrange(grobs=gs[select_grobs(lay)], layout_matrix=lay)
ggsave('Figure6.pdf',plot = p, width = 10, height = 8)



