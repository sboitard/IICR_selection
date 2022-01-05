# reads a table providing the pairwise coalescence times (T2) for all non recombining segments of several independent loci, and plots the empirical IICR (per locus and / or averaged over all loci) as in Chikhi et al (2018).

library(tidyverse)

estim.iicr <- function (u,N,E) {
	# u: table of simulated pairwise coalescence times for non recombining segments of several independent samples
	# N: diploid population size
	# E: number of time windows used to discretize the IICR 
	nb_s=max(u$sample)

	# computes the IICR for each locus
	T=5
	v_t=0.1*exp((1:E)/E*log(1+10*T))-0.1 # time vector (when the IICR will be computed)
	nb_t=length(v_t)
	res=data.frame(matrix(nrow=nb_s*(nb_t),ncol=3))
	colnames(res)=c('sample','t','IICR')
	i=1
	for (i_s in 1:nb_s){
		utemp=u %>% filter(sample==i_s)
		for (i_t in 1:(nb_t)){
			res$sample[i]=i_s
			res$t[i]=v_t[i_t]
			ind_P=which(utemp$time>=v_t[i_t])
			ind_f=which((utemp$time>=v_t[i_t-1])&(utemp$time<v_t[i_t+1]))
			if (length(ind_f)>0){
				Num=sum(utemp$length[ind_P])
				Den=sum(utemp$length[ind_f])/(v_t[i_t+1]-v_t[i_t-1])
				res$IICR[i]=Num/Den
			} else {
				res$IICR[i]=NA
			}
			i=i+1
		}
	}

	# computes the average IICR over all loci
	res_av=data.frame(matrix(nrow=nb_s*(nb_t),ncol=2))
	colnames(res_av)=c('t','IICR_av')
	for (i_t in 1:(nb_t)){
		res_av$t[i_t]=v_t[i_t]
		ind_P=which(u$time>=v_t[i_t])
		ind_f=which((u$time>=v_t[i_t-1])&(u$time<v_t[i_t+1]))
		if (length(ind_f)>0){
			Num=sum(u$length[ind_P])
			Den=sum(u$length[ind_f])/(v_t[i_t+1]-v_t[i_t-1])
			res_av$IICR_av[i_t]=Num/Den
		} else {
			res_av$IICR_av[i_t]=NA
		}
	}
	res=left_join(res,res_av)

	# rescale time in generations
	res$t=res$t*2*N

	# selects only a few samples to improve vizualization
	res = res %>% filter(sample<=5)
	
	return(res)
}

# Figure 6, bottom panels
u=read.table('schrider_1000_02_200rep.times',head=T)
u$time=2*u$time # rescale time in standard 2N units (msms unit is 4N)
res=estim.iicr(u,10000,25)
p1=ggplot(res,aes(x=t,y=IICR,color=as.factor(sample)))+geom_line()+theme_bw()+scale_x_continuous(limits=c(1000,100000),trans='log10')+scale_color_discrete(name = "region")+ylim(0,2)+geom_line(aes(x=t,y=IICR_av),color='black')+xlab('generations')
res=estim.iicr(u,10000,200)
p2=ggplot(res,aes(x=t,y=IICR,color=as.factor(sample)))+geom_line()+theme_bw()+scale_x_continuous(limits=c(1000,100000),trans='log10')+scale_color_discrete(name = "region")+ylim(0,2)+geom_line(aes(x=t,y=IICR_av),color='black')+xlab('generations')



