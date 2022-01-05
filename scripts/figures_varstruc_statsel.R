## IICR with K classes along the genome, each following a symmetrical n-island model whose deme size or migration rate can very over time. 
## Used in Figure 5 (or to replicate Figure 4)
## IICR is computed numerically from Equation (5) in Boitard et al (2022) combined with the Q matrix approach described in Rodriguez et al (2018).

library(tidyverse)
library(matlib)
library(latex2exp)
library(gridExtra)

# function to compute the Q matrix (or rate matrix) of a symmetrical n-island model
Qmat=function(M,n,c){
	# M: scaled migration rate
	# n: number of islands (or demes)
	# c: coalescence rate (inverse of the local deme size)
    	Q=matrix(0,nrow=3,ncol=3)
    	Q[1,2]=M
    	Q[1,3]=c
    	Q[1,1]=-M-c
    	Q[2,1]=M/(n-1)
    	Q[2,2]=-M/(n-1)
    	return(Q)
}

# generic function to compute the numerator and denominator of the IICR in a given class
coal.distrib=function(t,tmodif,lambda,M,n){
	# lambda: deme size in epoch 0,1,..,E (vector of size E+1)
	# M: scaled migration rate in epoch 0,1,..,E for a deme size of 1 (vector of size E+1)
	# tmodif: times where parameters change (vector of size E)
	# t: time where the IICR terms are evaluated (scalar)
	# n: number of demes
	ind=which(t>tmodif)
	if (length(ind)==0){ # t included in the first period
		Q=Qmat(M[1],n,1/lambda[1])
		temp=eigen(Q)
		Pt=temp$vectors%*%diag(exp(t*temp$values))%*%inv(temp$vectors)
		D=Pt%*%Q
		res=c(1-Pt[1,3],D[1,3])
	} else {
		iMax=max(ind)
		K=length(tmodif)
		d=tmodif-c(0,tmodif[1:(K-1)])
		Pt=diag(rep(1,3))
		for (i in 1:iMax){
			Q=Qmat(M[i],n,1/lambda[i])
			temp=eigen(Q)
			Pt=Pt%*%temp$vectors%*%diag(exp(d[i]*temp$values))%*%inv(temp$vectors)
		}
		Q=Qmat(M[iMax+1],n,1/lambda[iMax+1])
		temp=eigen(Q)
		Pt=Pt%*%temp$vectors%*%diag(exp((t-tmodif[iMax])*temp$values))%*%inv(temp$vectors)
		D=Pt%*%Q
		res=c(1-Pt[1,3],D[1,3])
	}
	return(res)
}

## Test 1: stationary n-island model (K=1)
n=10
lambda=rep(1,2)
M=rep(10,2)
tmodif=c(5)
t=seq(from=0,to=12,by=0.01)
IICR=rep(0,length(t))
for (i in 1:length(t)){
	temp=coal.distrib(t[i],tmodif,lambda,M,n)
	IICR[i]=temp[1]/temp[2]
}
plot(t,IICR)

## Test 2: non stationary n-island model (K=1), similar to that in Mazet et al (2016) for humans.
n=10
tmodif=c(6000,12000,30000)/1000
lambda=rep(1,4)
M=c(1,5,0.8,5)
t=seq(from=0,to=100,by=0.1)
IICR=rep(0,length(t))
for (i in 1:length(t)){
	temp=coal.distrib(t[i],tmodif,lambda,M,n)
	IICR[i]=temp[1]/temp[2]
}
plot(t*1000,IICR)

## Test 3: replicates Figure 4 (stationary n-island model with K=2) obtained analytically in figures_statstruc_statsel.R

# generic function to compute the numerator and denominator of the IICR for a stationary model
# simplified version of coal.distrib()
coal.distrib.stat=function(t,lambda,M,n){
	# lambda: deme size (scalar)
	# M: scaled migration rate for a deme size of 1 (scalar)
	# t: time where the IICR terms are evaluated (scalar)
	# n: number of demes
	Q=Qmat(M,n,1/lambda)
	temp=eigen(Q)
	Pt=temp$vectors%*%diag(exp(t*temp$values))%*%inv(temp$vectors)
	D=Pt%*%Q
	res=c(1-Pt[1,3],D[1,3])
	return(res)
}

n=10
N=1000
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
			res1=coal.distrib.stat(res$t[i],lambda[1],res$M[i],n)
			res2=coal.distrib.stat(res$t[i],lambda[2],res$M[i],n)
			num=a[1]*res1[1]+a[2]*res2[1]
			denom=a[1]*res1[2]+a[2]*res2[2]
			res$IICR[i]=num/denom
			i=i+1
		}
	}
}
# rescale t with respect to total meta-population size
res$t=res$t/n
# change panel names
M.labs=paste("M=",v_M,sep='')
names(M.labs)=v_M
# plots in natural scale
p=ggplot(res,aes(x=t,y=IICR,color=as.factor(a2)))+geom_line()+theme_bw()+theme(legend.position='bottom')+ylim(0,max(res$IICR))+scale_color_discrete(name = TeX("$a_2$")) +facet_grid(~M,labeller=labeller(M =M.labs))+xlim(0,10)+scale_x_log10()
ggsave('Figure4_verif.pdf',plot = p, width = 8, height = 4)

## Figure 5a-b - non stationary n-island model of Arredondo et al (2021), K=3 classes

# neutral parameters
N=1380 # diploid deme size
n=11
tmodif=c(24437,82969,107338,179666)/(2*N)
lambda_n=rep(1,5)
M=c(0.905519,17.7122,2.49889,0.721108,1.09619)
# selection parameters
v_a1=c(0,0.01,0.1,0.2,0.3,0.5,0.8) # possible values for the proportion of class 1
lambda=cbind(lambda_n/10,lambda_n,3*lambda_n)
# builds the table of IICR values
v_t=seq(from=0,to=10*n,by=0.1) # time vector (when the IICR will be computed)
p=length(v_t)*length(v_a1)
res=data.frame(matrix(nrow=p,ncol=3))
colnames(res)=c('a1','t','IICR')
i=1
for (i_t in 1:length(v_t)){
	for (i_a1 in 1:length(v_a1)){
		res$a1[i]=v_a1[i_a1]
		a=c(res$a1[i],1-res$a1[i]-0.01,0.01)
		temp=matrix(nrow=2,ncol=length(a))
		res$t[i]=v_t[i_t]
		for (j in 1:length(a)){
			temp[,j]=coal.distrib(res$t[i],tmodif,lambda[,j],M,n)
		}
		res$IICR[i]=sum(a*temp[1,])/sum(a*temp[2,])
		i=i+1
	}
}
# add a neutral IICR for comparison
res2=data.frame(matrix(nrow=length(v_t),ncol=2))
colnames(res2)=c('t','neutral')
for (i in 1:length(v_t)){		
	res2$t[i]=v_t[i]
	temp=coal.distrib(res2$t[i],tmodif,lambda[,2],M,n)
	res2$neutral[i]=temp[1]/temp[2]
}
res=left_join(res,res2)
# rescale t with respect to local deme size
res$t=res$t*2*N
# plots in log scale
pa=ggplot(res,aes(x=t,y=IICR,color=as.factor(a1)))+geom_line()+theme_bw()+scale_x_continuous(trans='log10',limits=c(100,300000))+scale_color_discrete(name = TeX("$a_1$"))+ylim(0,40)+geom_line(aes(x=t,y=neutral),color='black')+xlab('generations')

# similar plot varying a3 instead of a1
N=1380
n=11
tmodif=c(24437,82969,107338,179666)/(2*N)
lambda_n=rep(1,5)
M=c(0.905519,17.7122,2.49889,0.721108,1.09619)
v_a3=c(0,0.0001,0.001,0.01,0.1) # possible values for a1
lambda=cbind(lambda_n/10,lambda_n,3*lambda_n)
v_t=seq(from=0,to=10*n,by=0.1)
p=length(v_t)*length(v_a3)
res=data.frame(matrix(nrow=p,ncol=3))
colnames(res)=c('a3','t','IICR')
i=1
for (i_t in 1:length(v_t)){
	for (i_a3 in 1:length(v_a3)){
		res$a3[i]=v_a3[i_a3]
		a=c(0.5,0.5-res$a3[i],res$a3[i])
		temp=matrix(nrow=2,ncol=length(a))
		res$t[i]=v_t[i_t]
		for (j in 1:length(a)){
			temp[,j]=coal.distrib(res$t[i],tmodif,lambda[,j],M,n)
		}
		res$IICR[i]=sum(a*temp[1,])/sum(a*temp[2,])
		i=i+1
	}
}
res2=data.frame(matrix(nrow=length(v_t),ncol=2))
colnames(res2)=c('t','neutral')
for (i in 1:length(v_t)){		
	res2$t[i]=v_t[i]
	temp=coal.distrib(res2$t[i],tmodif,lambda[,2],M,n)
	res2$neutral[i]=temp[1]/temp[2]
}
res=left_join(res,res2)
res$t=res$t*2*N
pb=ggplot(res,aes(x=t,y=IICR,color=as.factor(a3)))+geom_line()+theme_bw()+scale_x_continuous(trans='log10',limits=c(100,300000))+scale_color_discrete(name = TeX("$a_3$"))+ylim(0,40)+geom_line(aes(x=t,y=neutral),color='black')+xlab('generations')

## Figure 5c - non stationary n-island model of Arredondo et al (2021), K=25 classes from Gossmann et al

# neutral parameters
N=1380 # diploid deme size
n=11
tmodif=c(24437,82969,107338,179666)/(2*N)
lambda_n=rep(1,5)
M=c(0.905519,17.7122,2.49889,0.721108,1.09619)
# selection parameters
obs=rlnorm(100000,-(0.682**2)/2,0.682)
h=hist(obs[which(obs<=5)],nclass=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
lambda=array(dim=c(5,length(hres$a)),data=0)
for (j in 1:length(hres$a)){
	lambda[,j]=lambda_n*hres$lambda[j]/hres$lambda[6]
}
# builds the table of IICR values
v_t=seq(from=0,to=10*n,by=0.1) # time vector (when the IICR will be computed)
p=length(v_t)
res=data.frame(matrix(nrow=p,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){	
	temp=matrix(nrow=2,ncol=length(hres$a))
	res$t[i]=v_t[i]
	for (j in 1:length(hres$a)){
		temp[,j]=coal.distrib(res$t[i],tmodif,lambda[,j],M,n)
	}
	res$IICR[i]=sum(hres$a*temp[1,])/sum(hres$a*temp[2,])
}
# add a neutral IICR for comparison
res2=data.frame(matrix(nrow=length(v_t),ncol=2))
colnames(res2)=c('t','neutral')
for (i in 1:length(v_t)){		
	res2$t[i]=v_t[i]
	temp=coal.distrib(res2$t[i],tmodif,lambda_n,M,n)
	res2$neutral[i]=temp[1]/temp[2]
}
res=left_join(res,res2)
# rescale t with respect to local deme size
res$t=res$t*2*N
# plots in log scale
pc=ggplot(res,aes(x=t,y=IICR))+geom_line(color='red')+theme_bw()+scale_x_continuous(limits=c(100,300000),trans='log10')+ylim(0,40)+geom_line(aes(x=t,y=neutral),color='black')+xlab('generations')

# combines the 3 panels
gs=list(pa,pb,pc)
lay=rbind(c(1,1,2,2),c(NA,3,3,NA))
select_grobs <- function(lay) {
	id <- unique(c(t(lay))) 
  	id[!is.na(id)]
} 
p=grid.arrange(grobs=gs[select_grobs(lay)], layout_matrix=lay)
ggsave('Figure5.pdf',plot = p, width = 10, height = 8)

## Figure S7: similar to Figure 5c, removing large values of lambda (as if regions under balancing selection were filtered out)

# figure S7a
N=1380
n=11
tmodif=c(24437,82969,107338,179666)/(2*N)
lambda_n=rep(1,5)
M=c(0.905519,17.7122,2.49889,0.721108,1.09619)
obs=rlnorm(100000,-(0.682**2)/2,0.682)
h=hist(obs[which(obs<=5)],nclass=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
# remove high lambdas
hres=hres%>%filter(lambda<=2)
hres$a=hres$a/sum(hres$a)
# again similar as above
lambda=array(dim=c(5,length(hres$a)),data=0)
for (j in 1:length(hres$a)){
	lambda[,j]=lambda_n*hres$lambda[j]/hres$lambda[6]
}
v_t=seq(from=0,to=10*n,by=0.1)
p=length(v_t)
res=data.frame(matrix(nrow=p,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){	
	temp=matrix(nrow=2,ncol=length(hres$a))
	res$t[i]=v_t[i]
	for (j in 1:length(hres$a)){
		temp[,j]=coal.distrib(res$t[i],tmodif,lambda[,j],M,n)
	}
	res$IICR[i]=sum(hres$a*temp[1,])/sum(hres$a*temp[2,])
}
res2=data.frame(matrix(nrow=length(v_t),ncol=2))
colnames(res2)=c('t','neutral')
for (i in 1:length(v_t)){		
	res2$t[i]=v_t[i]
	temp=coal.distrib(res2$t[i],tmodif,lambda_n,M,n)
	res2$neutral[i]=temp[1]/temp[2]
}
res=left_join(res,res2)
res$t=res$t*2*N
pa=ggplot(res,aes(x=t,y=IICR))+geom_line(color='red')+theme_bw()+scale_x_continuous(limits=c(100,300000),trans='log10')+ylim(0,40)+geom_line(aes(x=t,y=neutral),color='black')+xlab('generations')

# figure S7b
N=1380
n=11
tmodif=c(24437,82969,107338,179666)/(2*N)
lambda_n=rep(1,5)
M=c(0.905519,17.7122,2.49889,0.721108,1.09619)
obs=rlnorm(100000,-(0.682**2)/2,0.682)
h=hist(obs[which(obs<=5)],nclass=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
# remove high lambdas
hres=hres%>%filter(lambda<=3)
hres$a=hres$a/sum(hres$a)
# again similar as above
lambda=array(dim=c(5,length(hres$a)),data=0)
for (j in 1:length(hres$a)){
	lambda[,j]=lambda_n*hres$lambda[j]/hres$lambda[6]
}
v_t=seq(from=0,to=10*n,by=0.1)
p=length(v_t)
res=data.frame(matrix(nrow=p,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){	
	temp=matrix(nrow=2,ncol=length(hres$a))
	res$t[i]=v_t[i]
	for (j in 1:length(hres$a)){
		temp[,j]=coal.distrib(res$t[i],tmodif,lambda[,j],M,n)
	}
	res$IICR[i]=sum(hres$a*temp[1,])/sum(hres$a*temp[2,])
}
res2=data.frame(matrix(nrow=length(v_t),ncol=2))
colnames(res2)=c('t','neutral')
for (i in 1:length(v_t)){		
	res2$t[i]=v_t[i]
	temp=coal.distrib(res2$t[i],tmodif,lambda_n,M,n)
	res2$neutral[i]=temp[1]/temp[2]
}
res=left_join(res,res2)
res$t=res$t*2*N
pb=ggplot(res,aes(x=t,y=IICR))+geom_line(color='red')+theme_bw()+scale_x_continuous(limits=c(100,300000),trans='log10')+ylim(0,40)+geom_line(aes(x=t,y=neutral),color='black')+xlab('generations')

# combines plots
p=grid.arrange(pa,pb,ncol=2)
ggsave('FigureS7.pdf',plot = p, width = 10, height = 4)

