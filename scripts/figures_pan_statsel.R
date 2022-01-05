## IICR with K classes of Ne along the genome, all panmictic and stationary (i.e. constant over time)
## Used for Figures 1-2
## IICR is computed analytically from Equation (2) of Boitard et al (2022)

library(tidyverse)
library(latex2exp)
library(gridExtra)

# generic function to compute the IICR with K classes
IICR.pan.varN <- function (time, lambda.freq, lambda) {
	# time: vector of times points when the IICR is computed
	# lambda.freq: proportion of each class
	# lambda: Ne in each class
	R <- sum(lambda.freq*exp(-time/lambda))
   	f <- sum( (lambda.freq/lambda) * exp(-time/lambda) )
   	return(R/f)
}

## Figure 1 - K=2 classes
lambda=c(0.1,1) # relative Ne in each class
v_t=seq(from=0,to=10,by=0.01) # time vector (when the IICR will be computed)
v_a2=c(0,0.01,0.1,0.5,0.9,0.99,1) # possible values for the proportion of class 2
n=length(v_a2)*length(v_t)
# builds the table of IICR values
res=data.frame(matrix(nrow=n,ncol=3))
colnames(res)=c('a2','t','IICR')
i=1
for (i_a2 in 1:length(v_a2)){
	for (i_t in 1:length(v_t)){
		res$a2[i]=v_a2[i_a2]
		res$t[i]=v_t[i_t]
		a=c(1-res$a2[i],res$a2[i]) # proportion of each class
		res$IICR[i]=IICR.pan.varN(res$t[i],a,lambda)
		i=i+1
	}
}
# plots in log scale
p=ggplot(res,aes(x=t,y=IICR,color=as.factor(a2)))+geom_line()+theme_bw()+ylim(0,lambda[2])+scale_color_discrete(name = TeX("$a_2$"))+xlim(0,10)+scale_x_log10()
ggsave('Figure1_log.pdf',plot = p, width = 6, height = 4)
# similar plot with rescaled values of lambda
lambda=c(1,10)
v_t=seq(from=0,to=100,by=0.1)
v_a2=c(0,0.01,0.1,0.5,0.9,0.99,1)
n=length(v_a2)*length(v_t)
res=data.frame(matrix(nrow=n,ncol=3))
colnames(res)=c('a2','t','IICR')
i=1
for (i_a2 in 1:length(v_a2)){
	for (i_t in 1:length(v_t)){
		res$a2[i]=v_a2[i_a2]
		res$t[i]=v_t[i_t]
		a=c(1-res$a2[i],res$a2[i])
		res$IICR[i]=IICR.pan.varN(res$t[i],a,lambda)
		i=i+1
	}
}
p=ggplot(res,aes(x=t,y=IICR,color=as.factor(a2)))+geom_line()+theme_bw()+ylim(0,lambda[2])+scale_color_discrete(name = TeX("$a_2$"))+xlim(0,100)+scale_x_log10()
ggsave('Figure1_log_rescaled.pdf',plot = p, width = 6, height = 4)

## Figure 2 - K=3 classes
lambda=c(0.1,1,3) # relative Ne in each class
v_t=seq(from=0,to=10,by=0.01) # time vector (when the IICR will be computed)
v_a1=c(0,0.01,0.1,0.2,0.3,0.5,0.8) # possible values for the proportion of class 1
n=length(v_a1)*length(v_t)
# builds the table of IICR values
res=data.frame(matrix(nrow=n,ncol=3))
colnames(res)=c('a1','t','IICR')
i=1
for (i_a1 in 1:length(v_a1)){
	for (i_t in 1:length(v_t)){
		res$a1[i]=v_a1[i_a1]
		res$t[i]=v_t[i_t]
		a=c(res$a1[i],1-res$a1[i]-0.01,0.01) # proportion of each class
		res$IICR[i]=IICR.pan.varN(res$t[i],a,lambda)
		i=i+1
	}
}
# plots in natural scale
p1=ggplot(res,aes(x=t,y=IICR,color=as.factor(a1)))+geom_line()+theme_bw()+xlim(0,10)+ylim(0,lambda[3]) +scale_color_discrete(name = TeX("$a_1$"))
# plots in log scale
p1_log=ggplot(res,aes(x=t,y=IICR,color=as.factor(a1)))+geom_line()+theme_bw()+xlim(0,10)+ylim(0,lambda[3]) +scale_color_discrete(name = TeX("$a_1$"))+scale_x_log10()
# similars plot varying a3 instead of a1
lambda=c(0.1,1,3)
v_t=seq(from=0,to=10,by=0.01)
v_a3=c(0,0.0001,0.001,0.01,0.1) # possible values for the proportion of class 1
n=length(v_a3)*length(v_t)
res=data.frame(matrix(nrow=n,ncol=3))
colnames(res)=c('a3','t','IICR')
i=1
for (i_a3 in 1:length(v_a3)){
	for (i_t in 1:length(v_t)){
		res$a3[i]=v_a3[i_a3]
		res$t[i]=v_t[i_t]
		a=c(0.5,0.5-res$a3[i],res$a3[i])
		res$IICR[i]=IICR.pan.varN(res$t[i],a,lambda)
		i=i+1
	}
}
p2=ggplot(res,aes(x=t,y=IICR,color=as.factor(a3)))+geom_line()+theme_bw()+xlim(0,10)+ylim(0,lambda[3]) +scale_color_discrete(name = TeX("$a_3$"))
p2_log=ggplot(res,aes(x=t,y=IICR,color=as.factor(a3)))+geom_line()+theme_bw()+xlim(0,10)+ylim(0,lambda[3]) +scale_color_discrete(name = TeX("$a_3$"))+scale_x_log10()
# combine plots
p=grid.arrange(p1,p2,ncol=2)
ggsave('Figure2.pdf',plot = p, width = 10, height = 4)
p_log=grid.arrange(p1_log,p2_log,ncol=2)
ggsave('Figure2_log.pdf',plot = p_log, width = 10, height = 4)

## Figure 3 - K=25 classes defined Elyashiv et al results

# load local Ne data and compute empirical distribution
# these files can be obatine at https://github.com/sellalab/LinkedSelectionMaps/
u1=read.table('melanogaster_maps/LSmap_BS1234_SW123_a5_2L.LS',comment.char="",skip=2,head=F)
u2=read.table('melanogaster_maps/LSmap_BS1234_SW123_a5_2R.LS',comment.char="",skip=2,head=F)
u3=read.table('melanogaster_maps/LSmap_BS1234_SW123_a5_3L.LS',comment.char="",skip=2,head=F)
u4=read.table('melanogaster_maps/LSmap_BS1234_SW123_a5_3R.LS',comment.char="",skip=2,head=F)
u=rbind(u1,u2,u3,u4)
# scales to get a mean value of 1
u[,2]=u[,2]/mean(u[,2])
u = u %>% mutate(curve='Ne distribution') #%>% mutate(scenario='Elyashiv et al (2016), Drosophila')
# plots the densitystrip.text.y = element_text(size=12)
p1=ggplot(u,aes(x=V2))+geom_line(stat="density")+theme_bw()+xlab(TeX("$\\lambda$"))+ylim(0,1)+xlim(0,5)+facet_grid(~curve)+theme(strip.text.x = element_text(size=12))
# breaks into 25 classes of lambda
h=hist(u[,2],breaks=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
# builds the table of IICR values
v_t=seq(from=0,to=500,by=0.01) # time vector (when the IICR will be computed)
n=length(v_t)
res=data.frame(matrix(nrow=n,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){
	res$t[i]=v_t[i]		
	res$IICR[i]=IICR.pan.varN(res$t[i],hres$a,hres$lambda)
}

# plots in log scale scale for t <= 10
res = res %>% mutate(curve='Recent IICR')
p2=ggplot(res%>%filter(t<=10),aes(x=t,y=IICR))+geom_line()+theme_bw()+scale_x_log10()+ylim(0,5)+facet_grid(~curve)+theme(strip.text.x = element_text(size=12))
# plots in log scale for all t
res = res %>% mutate(curve='Full IICR') %>% mutate(scenario='Elyashiv et al (2016), Drosophila')
p3=ggplot(res,aes(x=t,y=IICR))+geom_line()+theme_bw()+scale_x_log10()+ylim(0,5)+facet_grid(scenario~curve)+theme(strip.text.y = element_text(size=12),strip.text.x = element_text(size=12))

# similar plots after removing the 'mode' in 0
ind=which(u[,2]>=0.25)
u[ind,2]=u[ind,2]/mean(u[ind,2])
p10=ggplot(u[ind,],aes(x=V2))+geom_line(stat="density")+theme_bw()+xlab(TeX("$\\lambda$"))+xlim(0,5)+ylim(0,1)+facet_grid(~curve)+theme(strip.text.x = element_text(size=12))
h=hist(u[ind,2],breaks=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
v_t=seq(from=0,to=500,by=0.01)
n=length(v_t)
res=data.frame(matrix(nrow=n,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){
	res$t[i]=v_t[i]		
	res$IICR[i]=IICR.pan.varN(res$t[i],hres$a,hres$lambda)
}
res = res %>% mutate(curve='Recent IICR')
p11=ggplot(res%>%filter(t<=10),aes(x=t,y=IICR))+geom_line()+theme_bw()+scale_x_log10()+ylim(0,5)+facet_grid(~curve)+theme(strip.text.x = element_text(size=12))
res = res %>% mutate(curve='Full IICR')
p12=ggplot(res,aes(x=t,y=IICR))+geom_line()+theme_bw()+scale_x_log10()+ylim(0,5)+facet_grid(~curve)+theme(strip.text.x = element_text(size=12))
# plot a column of figures
p=grid.arrange(p10,p11,p12,ncol=3)
ggsave('Figure3_elyashiv_nomode.pdf',plot = p, width = 12, height = 4)

## Figure 3 - K=25 classes defined Gossmann et al results

# compute density from their Table 1 (drosophila estimations)
x=seq(from=0,to=5,by=0.1)
d_th=dlnorm(x,-(0.743**2)/2,0.743)
res=data.frame(cbind(x,d_th))
p4=ggplot(res,aes(x=x,y=d_th))+geom_line()+theme_bw()+xlab(TeX("$\\lambda$"))+ylab('density')+ylim(0,1)+xlim(0,5)
# simulate lambda values and breaks into 25 classes of lambda
obs=rlnorm(100000,-(0.743**2)/2,0.743)
h=hist(obs[which(obs<=5)],nclass=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
# builds the table of IICR values
v_t=seq(from=0,to=500,by=0.01) # time vector (when the IICR will be computed)
n=length(v_t)
res=data.frame(matrix(nrow=n,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){
	res$t[i]=v_t[i]		
	res$IICR[i]=IICR.pan.varN(res$t[i],hres$a,hres$lambda)
}

# plots in log scale for t <= 10
res = res %>% mutate(scenario='Gossmann et al (2011), Drosophila')
p5=ggplot(res%>%filter(t<=10),aes(x=t,y=IICR))+geom_line()+theme_bw()+ scale_x_log10()+ylim(0,5)
# plots in log scale for all t
p6=ggplot(res,aes(x=t,y=IICR))+geom_line()+theme_bw() + scale_x_log10()+ylim(0,5)+facet_grid(scenario~.)+theme(strip.text.y = element_text(size=12))

# similar plots but for humans
x=seq(from=0,to=5,by=0.1)
d_th=dlnorm(x,-(0.682**2)/2,0.682)
res=data.frame(cbind(x,d_th))
p7=ggplot(res,aes(x=x,y=d_th))+geom_line()+theme_bw()+xlab(TeX("$\\lambda$"))+ylab('density')+ylim(0,1)+xlim(0,5)
obs=rlnorm(100000,-(0.682**2)/2,0.682)
h=hist(obs[which(obs<=5)],nclass=25)
hres=cbind(h$mid,h$counts/sum(h$counts))
colnames(hres)=c('lambda','a')
hres=data.frame(hres)
v_t=seq(from=0,to=500,by=0.01)
n=length(v_t)
res=data.frame(matrix(nrow=n,ncol=2))
colnames(res)=c('t','IICR')
for (i in 1:length(v_t)){
	res$t[i]=v_t[i]		
	res$IICR[i]=IICR.pan.varN(res$t[i],hres$a,hres$lambda)
}
res = res %>% mutate(scenario='Gossmann et al (2011), humans')
p8=ggplot(res%>%filter(t<=10),aes(x=t,y=IICR))+geom_line()+theme_bw()+ scale_x_log10()+ylim(0,5)
p9=ggplot(res,aes(x=t,y=IICR))+geom_line()+theme_bw() + scale_x_log10()+ylim(0,5)+facet_grid(scenario~.)+theme(strip.text.y = element_text(size=12))

# plot an array of figures
p=grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3,widths=c(9,9,10),heights=c(10,9,9))
ggsave('FigureS3.pdf',plot = p, width = 12, height = 12)

## Variance of T2
Vexp <- function (a,lambda) {
	# a: proportion of each class
	# lambda: Ne in each class
	S1=sum(a*lambda)
	S2=sum(a*lambda**2)
	V=2*S2-S1**2
	return(sqrt(V)/sum(a*lambda))
}
lambda=c(0.1,1)
a=c(0.1,0.9)
Vexp(a,lambda)


