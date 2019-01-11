##Test model with hospitalization data


#for(ag.select in c('2-59m','2-11m','12-23m','24-59m','<2m')){
#for(ag.select in c('2-11m','12-23m','24-59m','<2m', '2-23m')){
ag.select<-'2-59m'
subnat='F'
    print(ag.select)
    print(subnat)
rm(list=ls()[-which(ls() %in% c('ag.select', 'subnat'))]) #for instance 

input.dir <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/nat'   #Directory where results will be saved.
input.dir <- paste(input.dir,  '/',ag.select,'/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.

sim.counts<-readRDS( file=paste0(input.dir,"sim.counts.", ag.select,".rds"))

library(RCurl)
library(reshape2)
library(lubridate)
library(abind)
library(dummies)

max.time.points=48

#Filter to obtain relevant age groups
#prelog_data<-prelog_data[substr(prelog_data$age_group,4,8)==ag.select,]  #Only <12m
if(subnat){
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/subnat'   #Directory where results will be saved.
  
}else{
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/nat'   #Directory where results will be saved.
}
if(ag.select=='<2m'){ ag.select<-'u2m'}
output_directory <- paste(output_directory,  '/',ag.select,'.SIM','/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)


##########################################################################################################################################
##########################################################################################################################################
##########################################################################################################################################
######################################################################
###############NEXT PERFORM THE SMOOTHING META_MODEL
######################################################################
##########################################################################################################################################
##########################################################################################################################################

#for(ag.select in c('2-59m','2-11m','12-23m','24-59m','<2m', '2-24m')){
#for(ag.select in c('12-23m','24-59m','<2m')){
  
  
  log.rr<- readRDS( file=paste0(output_directory,"log.rr.0.05.q_", ag.select,".rds"))
  log.rr.sd<- readRDS( file=paste0(output_directory,"log.rr.0.05.sd_", ag.select,".rds"))
  
  aux.output1<- readRDS( file=paste0(output_directory,"in.data.", ag.select,".rds"))
  
  log.rr.prec<- 1/log.rr.sd^2
  log.rr.med<-log.rr['50%',,,]
  log.rr.med<-array(log.rr.med, dim=c(dim(log.rr.med),1))

  post.index.array<-aux.output1$post.index.array
  q<-4
  I_Omega<-diag(q)
  N.countries<-aux.output1$N.countries
  N.states.country<-aux.output1$N.states.country
  ts.length_mat<-aux.output1$n.times
  pre.vax.time<-aux.output1$n.times.pre  
    
  N.states<-aux.output1$N.states.country
  time.index<-aux.output1$post.index.array
  
  #########################
  #CALL JAGS POOLING MODEL
  #########################
  max.time.points=48
  source('model.pool.log.rr.R')
  
  countries<-aux.output1$countries
  
  
  ##########################################################################################
  ##### test for Convergence - trace plots (specify which to test in "variable.names") #####
  ##########################################################################################
  #par(ask=TRUE)
  #dev.off()
  #pdf(file = "Trace Plots.pdf")
  #par(mfrow = c(3,2))
  #plot(posterior_samples)
  #dev.off()
  
  ##################################################################
  ##### model fit (DIC) - must change n.chains=2 in model_jags #####
  ##################################################################
  # dic.samples(model_jags, 500) #with pooling josh's fix (2 runs w 500 each): deviance: -717, -715; penalty=269,266, DIC=-447, -449
  #dic.samples(model_jags, 500)
  
  # par(mfrow=c(1,1))
  # caterplot(posterior_samples)
  
  #############################################################################
  ##### create posterior estimates into data frames (w.true and reg.mean) #####
  #############################################################################
  #############################################################################
  ##### create posterior estimates into data frames (w.true and reg.mean) #####
  #############################################################################
  x.func <- function(x){ sub(".*\\[(.*)\\].*", "\\1", x, perl=TRUE)} 
  
  beta1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
  beta1.lab<-x.func(dimnames(beta1)[[2]]) 
  beta1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
  names(beta1.lab.spl)<-c('country','state')
  
  #Extract first change point time:
  beta3<-posterior_samples[[1]][,grep("^beta.*,3]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
  cp1<-max.time.points*exp(beta3) #Intercept
  beta3.lab<-x.func(dimnames(beta3)[[2]]) 
  beta3.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta3.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
  names(beta3.lab.spl)<-c('country','state')
  #for(i in c(1:ncol(beta3))){hist(beta3[,i])}
  quant.cp1<-t(max.time.points*exp(apply(beta3,2,quantile, probs=c(0.025,0.5,0.975))))
  var.cp1<-apply(beta3,2,var)
  cp1.lab<-x.func(dimnames(quant.cp1)[[1]]) 
  cp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(cp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
  names(cp1.lab.spl)<-c('country','state')
  par(mfrow=c(2,2))
  plot(y=1:nrow(quant.cp1), x=quant.cp1[,'50%'], bty='l')
  library(ggplot2)
  plot.cp<-cbind.data.frame('strata'=1:nrow(quant.cp1),'median.cp'=quant.cp1[,'50%'],'lcl.cp'=quant.cp1[,'2.5%'],'ucl.cp'=quant.cp1[,'97.5%'], 'inv.var.cp'=1/var.cp1,beta3.lab.spl)
  plot.cp<-plot.cp[order(plot.cp$country),]
  plot.cp$order2<-1:nrow(plot.cp)
  plot.cp$country2<-NA
  for(i in 1:length(countries)){plot.cp$country2[plot.cp$country==i] <-countries[i] }
  tiff(paste0(output_directory,'cp by country.log.rr.tiff'), width = 7, height = 8, units = "in",res=200)
  cp.plot<-ggplot(data=plot.cp, aes(x=median.cp, y=order2, color=country2)) +
    geom_point(aes(size=inv.var.cp)) +
    scale_size_continuous(range=c(1,15)) +
    theme_bw()+
    guides( size = FALSE)+
    # theme(legend.position = "none")+
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
  print(cp.plot)
  dev.off()
  saveRDS(plot.cp, file=paste0(output_directory,'CP1 pool model.log.rr.rds'))
  
  
  ############################
  #Extract slope
  beta2<-posterior_samples[[1]][,grep("^beta.*,2]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
  slp1<-beta2 #Intercept
  beta2.lab<-x.func(dimnames(beta2)[[2]]) 
  beta2.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta2.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
  names(beta2.lab.spl)<-c('country','state')
  #for(i in c(1:9)){hist(beta2[,i])}
  quant.slp1<-t(apply(beta2,2,quantile, probs=c(0.025,0.5,0.975)))
  var.slp1<-apply(beta2,2,var)
  slp1.lab<-x.func(dimnames(quant.slp1)[[1]]) 
  slp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(slp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
  names(slp1.lab.spl)<-c('country','state')
  par(mfrow=c(2,2))
  plot(y=1:nrow(quant.slp1), x=quant.slp1[,'50%'], bty='l')
  library(ggplot2)
  plot.slp<-cbind.data.frame('strata'=1:nrow(quant.slp1),'median.slp'=quant.slp1[,'50%'], 'inv.var.slp'=1/var.slp1,beta2.lab.spl)
  plot.slp<-plot.slp[order(plot.slp$country),]
  plot.slp$order2<-1:nrow(plot.slp)
  plot.slp$country2<-NA
  for(i in 1:length(countries)){plot.slp$country2[plot.slp$country==i] <-countries[i] }
  tiff(paste0(output_directory,'country.slope.log.rr.tiff'), width = 7, height = 8, units = "in",res=200)
  slp.plot<-ggplot(data=plot.slp, aes(x=median.slp, y=order2, color=country2)) +
    geom_point(aes(size=inv.var.slp)) +
    scale_size_continuous(range=c(1,15)) +
    theme_bw()+
    guides( size = FALSE)+
    # theme(legend.position = "none")+
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
  print(slp.plot)
  dev.off()
  saveRDS(plot.slp, file=paste0(output_directory,"slope pool model.log.rr.rds"))
  
  ##########################################
  
  #Relationship between change point location and slope
  tiff(paste0(output_directory,'cp vs slope.log.rr.tiff'), width = 7, height = 8, units = "in",res=200)
  par(mfrow=c(1,1), mar=c(4,4,1,1))
  plot(quant.slp1[,'50%'], quant.cp1[,'50%'], col=beta2.lab.spl$country,bty='l', ylab="Change point", xlab="Slope")
  dev.off()
  
  ##melt and cast predicted values into 4D array N,t,i,j array
  reg_mean<-posterior_samples[[1]][,grep("reg_mean",dimnames(posterior_samples[[1]])[[2]])]
  pred1<-reg_mean
  pred.indices<- x.func(dimnames(pred1)[[2]]) 
  pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
  pred.indices.spl<-as.data.frame(pred.indices.spl)
  names(pred.indices.spl)<- c('country','state','time')
  reg_mean2<-cbind.data.frame(pred.indices.spl,t(reg_mean))
  reg_mean_m<-melt(reg_mean2, id=c('time','country','state'))
  reg_mean_c<-acast(reg_mean_m, variable~time~country~state)
  preds.q<-apply(reg_mean_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
  dimnames(preds.q)[[2]]<- as.numeric(as.character(dimnames(preds.q)[[2]]))
  
  #Unbias
  reg_unbias<-posterior_samples[[1]][,grep("reg_unbias",dimnames(posterior_samples[[1]])[[2]])]
  pred1<-reg_unbias
  pred.indices<- x.func(dimnames(pred1)[[2]]) 
  pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
  pred.indices.spl<-as.data.frame(pred.indices.spl)
  names(pred.indices.spl)<- c('country','state','time')
  reg_unbias2<-cbind.data.frame(pred.indices.spl,t(reg_unbias))
  reg_unbias_m<-melt(reg_unbias2, id=c('time','country','state'))
  reg_unbias_c<-acast(reg_unbias_m, variable~time~country~state)
  reg_unbias_c<-reg_unbias_c[,order(as.numeric(dimnames(reg_unbias_c)[[2]])),order(as.numeric(dimnames(reg_unbias_c)[[3]])),order(as.numeric(dimnames(reg_unbias_c)[[4]]))]
  preds.unbias.q<-apply(reg_unbias_c,c(2,3),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
  dimnames(preds.unbias.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))
  
  tiff(paste0(output_directory,'country.trend.log.rr.0.05.tiff'), width = 7, height = 8, units = "in",res=200)
  par(mfrow=c(5,2), mar=c(4,2,1,1))
  for(i in 1:length(countries)){
    # for(j in 1:N.states[i]){
    plot.data<-t(preds.unbias.q[,,i])
    if( abs(sum(plot.data, na.rm=T))>0){
      plot.data<-plot.data[complete.cases(plot.data),]
      final.rr<-paste0(round(exp(plot.data[nrow(plot.data),'50%', drop=F]),2),
                       ' (' ,round(exp(plot.data[nrow(plot.data),'2.5%', drop=F]),2),',',
                       round(exp(plot.data[nrow(plot.data),'97.5%', drop=F]),2),")")
    }
    tot_time<-nrow(plot.data)
    max.x<-max(((1:tot_time)-pre.vax.time[i]), na.rm=T)
    matplot(post.index.array[i,1,][1:nrow(plot.data)]*max.time.points, plot.data+1,type='l',yaxt='n', xlim=c(0, 48), xlab='months post-PCV introduction', ylim=c(-0.2,1.5), col='black', lty=c(2,1,2), bty='l')
    abline(h=1, col='gray')
    axis(side=2, at=c(0,0.5,1, 1.5, 2), las=1,labels=T)
     abline(v=0)
    text(44, 1.4, final.rr)
    title(countries[i])
  }
  dev.off()
  
  saveRDS(preds.unbias.q, file=paste0(output_directory,"reg_mean_with_pooling CP nobias.log.rr.0.05.rds"))
  saveRDS(reg_unbias_c, file=paste0(output_directory, "reg_mean_unbias_with_pooling CP.log.rr.0.05.rds"))
  #saveRDS(state.labels, file=paste0(output_directory,"state labels.rds"))
  saveRDS(posterior_samples, file=paste0(output_directory, "posterior samples pooling with CP.log.rr.0.05.rds"))



