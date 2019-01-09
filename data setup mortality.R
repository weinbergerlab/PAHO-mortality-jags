##Test model with hospitalization data

for(ag.select in c('2-59m','2-11m','12-23m','24-59m','<2m')){
#for(ag.select in c('12-23m','24-59m','<2m', '2-23m')){
  for(subnat in c(FALSE)){
    print(ag.select)
    print(subnat)
rm(list=ls()[-which(ls() %in% c('ag.select', 'subnat'))]) #for instance 

library(RCurl)
library(reshape2)
library(lubridate)
library(abind)
library(dummies)

max.time.points=48


input_directory  <- "C:/Users/dmw63/Your team Dropbox/PAHO mortality/Data/" #Directory or URL page containing input data file
file_name="PAHO all age cuts_SubChapters.csv"
data_file <- paste0(input_directory, file_name)
#prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM LOCAL
prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM URL

#Filter to obtain relevant age groups
#prelog_data<-prelog_data[substr(prelog_data$age_group,4,8)==ag.select,]  #Only <12m
prelog_data<-prelog_data[grep(ag.select,prelog_data$age_group ),]
if(subnat){
  prelog_data<-prelog_data[!grepl(' A',prelog_data$age_group ),]  #filter out summary categories
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/subnat'   #Directory where results will be saved.
  
}else{
  prelog_data<-prelog_data[grepl(' A',prelog_data$age_group ),]  #filter out summary categories
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/nat'   #Directory where results will be saved.
}
if(ag.select=='<2m'){ ag.select<-'u2m'}
output_directory <- paste(output_directory,  '/',ag.select,'/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
prelog_data<-prelog_data[,c('age_group', 'monthdate','J12_J18_prim','acm_noj_prim','A00_B99_prim','P00_P96_prim',"E00_E90_prim" )]
prelog_data$monthdate<-as.Date(prelog_data$monthdate)
prelog_data$country<-substr(prelog_data$age_group,1,2)

##Set Vax intro date for each country
prelog_data$intro.date<-as.Date('1900-01-01')
prelog_data$intro.date[prelog_data$country=='ar'] <-as.Date('2010-03-01')
prelog_data$intro.date[prelog_data$country=='br']<-as.Date('2010-03-01')
prelog_data$intro.date[prelog_data$country=='co'] <-as.Date('2011-11-01')
prelog_data$intro.date[prelog_data$country=='dr'] <-as.Date('2013-09-01')
prelog_data$intro.date[prelog_data$country=='ec'] <-as.Date('2010-08-01')
prelog_data$intro.date[prelog_data$country=='gy'] <-as.Date('2011-01-01')
prelog_data$intro.date[prelog_data$country=='hr'] <-as.Date('2011-01-01')
prelog_data$intro.date[prelog_data$country=='mx'] <-as.Date('2008-01-01')
prelog_data$intro.date[prelog_data$country=='nc'] <-as.Date('2011-01-01')
prelog_data$intro.date[prelog_data$country=='pr']  <-as.Date('2009-08-01') 
prelog_data$interval.date<- prelog_data$intro.date %--% prelog_data$monthdate
prelog_data$post.index<-round(as.numeric(as.duration(prelog_data$interval.date),'months'))
prelog_data<-prelog_data[prelog_data$post.index<=max.time.points,]
countries<-unique(prelog_data$country)
prelog_data$J12_J18_prim.pre<-prelog_data$J12_J18_prim
prelog_data$J12_J18_prim.pre[prelog_data$post.index>0]<-NA
##

prelog_data$interval.date<-NULL
prelog_data$hdi<-substr(prelog_data$age_group,10,11)

#Assign a time index, starting at 1:N time points
prelog_data$age_group<-factor(prelog_data$age_group)
prelog_data_spl<-split(prelog_data, prelog_data$age_group)
prelog_data_spl<-lapply(prelog_data_spl, function(x) cbind.data.frame(x,index=1:nrow(x) )  )
prelog_data_spl<-do.call("rbind", prelog_data_spl)
prelog_data_spl<-split(prelog_data_spl, prelog_data_spl$country)
prelog_data_spl<-lapply(prelog_data_spl, function(x) cbind.data.frame(x,grp.index=as.numeric(factor(x$hdi)) ))  
prelog_data_spl<-do.call("rbind", prelog_data_spl)
prelog_data_spl
hdi.index<-unique(cbind.data.frame(prelog_data_spl$country, prelog_data_spl$grp.index, prelog_data_spl$hdi))
names(hdi.index)<-c('country','grp.index','hdi')
prelog_data_spl$country_index<-as.numeric(as.factor(prelog_data_spl$country))
prelog_data_spl<-prelog_data_spl[,c( 'index','country_index','grp.index',"J12_J18_prim","J12_J18_prim.pre",'acm_noj_prim','A00_B99_prim','P00_P96_prim',"E00_E90_prim","post.index","monthdate")]

#Reshape data into 3D arrays--separate arrays for outcome, controls, time index
ds.m<-melt(prelog_data_spl,id.vars=c( 'index','country_index','grp.index'))
outcome.array<- acast(ds.m[ds.m$variable=="J12_J18_prim",], country_index~grp.index~index )   #Outcome array
outcome.array.pre<- acast(ds.m[ds.m$variable=="J12_J18_prim.pre",], country_index~grp.index~index )   #Outcome array
control1.array<- acast(ds.m[ds.m$variable=="acm_noj_prim",],  country_index~grp.index~index )  #Control array
control1.array<-log(control1.array+0.5) #Log transform
control2.array<- acast(ds.m[ds.m$variable=="A00_B99_prim",],  country_index~grp.index~index )  #Control array
control2.array<-log(control2.array+0.5) #Log transform
control3.array<- acast(ds.m[ds.m$variable=="P00_P96_prim",],  country_index~grp.index~index )  #Control array
control4.array<-log(control3.array+0.5) #Log transform
control4.array<- acast(ds.m[ds.m$variable=="E00_E90_prim",],  country_index~grp.index~index )  #Control array
control4.array<-log(control4.array+0.5) #Log transform
post.index.array <- acast(ds.m[ds.m$variable=="post.index",],  country_index~grp.index~index ) #time index array
post.index.array[post.index.array<0]<-0 #Pre-vax indices are 0
post.index.array<-post.index.array/max.time.points


date.array<- as.Date(acast(ds.m[ds.m$variable=="monthdate",], country_index~grp.index~index ),origin="1970-01-01")   #Outcome array
month.array<-array(month(date.array)   , dim=dim(date.array))
month.dummy<-array(NA, dim=c(dim(month.array),11))
for(i in c(1:11)){
  month.dummy[,,,i]<- as.numeric(month.array==i & !is.na(month.array))
}

#Combine covariate array with intercept along 4th dimension
int.array<-array(1, dim(control1.array))
control1.array.int<-abind(int.array,control1.array,control2.array, control3.array, control4.array, along=4 )
dimnames(control1.array.int)[[4]]<-c('Intercept','ACM',"A00","P00", "E00")
control1.array.int2<-abind(control1.array.int,post.index.array, along=4 )
dimnames(control1.array.int2)[[4]]<-c('Intercept','ACM',"A00","P00", "E00","post.index")

#Create variables to be used by JAGS
N.countries=dim(control1.array.int2)[1]
n.times=apply(outcome.array,c(1,2), function(x) sum(!is.na(x))) #Vector (or matrix) giving number of non-missing observations per state
n.times.pre=apply(outcome.array.pre,c(1,2), function(x) sum(!is.na(x))) #Vector (or matrix) giving number of non-missing observations per state

N.states.country<-apply(n.times,1, function(x) sum(x!=0)  )

N.preds=3 #N predictors ( intercept, acm_noj, post-vax time trend)
q=6  #N predictors of slopes (if none, set to 1 for int only)
# I_Sigma<-replicate( N.countries, diag(N.preds) )
I_Sigma<-diag(N.preds) 
I_Omega<-diag(q)
I_Omicron<-diag(11)

#Z needs to be an array i,j,q ; with assignment of covariate based on hdi.index df
z<-array(0, dim=c(N.countries, max(N.states.country),q))
z[,,1]<-1 #intercept for z 
hdi.index$hdi.countryN<- as.numeric(as.factor(hdi.index$country))
for(i in 1:N.countries){
  for(j in 1:N.states.country[i]){
   hdi.countryN<- as.numeric(as.factor(hdi.index$country))
  if( sum(hdi.countryN==i & hdi.index$grp.index==j  & hdi.index$hdi=='Lo')>0){  
      z[i,j,2]<-1
  }else{if(sum(hdi.index$hdi.countryN==i & hdi.index$grp.index==j  & hdi.index$hdi=='Me')>0){
    z[i,j,3]<-1
  }}
  }
}

m<-1 
w<-matrix(NA, nrow=N.countries, ncol=m) 
for(i in 1:N.countries){
  for(k in 1:m){
    w[i,k]<-1
  }
}
offset<-array(0, dim=dim(outcome.array)) #no offset

aux.output<-list('post.index.array'=post.index.array, 'countries'=countries,
                 'n.times.pre'=n.times.pre,'n.times'=n.times,'N.countries'=N.countries,'N.states.country'=N.states.country )
#saveRDS(aux.output, file=paste0(output_directory,"aux.output", ag.select,".rds"), compress=FALSE)

###########################
#Call JAGS Model
source('model.R')
##########################


#Extract output############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
#############################################################################
##### create posterior estimates into data frames (w.true and reg.mean) #####
#############################################################################
x.func <-  function(x){ sub(".*\\[(.*)\\].*", "\\1", x, perl=TRUE)} 

beta1<-posterior_samples[[1]][,grep("^beta.*,1]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
beta1.lab<-x.func(dimnames(beta1)[[2]]) 
beta1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta1.lab.spl)<-c('country','state')


##melt and cast predicted values into 4D array N,t,i,j array
log.pred.mu<-posterior_samples[[1]][,grep("log.pred.mu",dimnames(posterior_samples[[1]])[[2]])]
pred.indices<- x.func(dimnames(log.pred.mu)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')

#predictive distributio
pred1<-rpois(n=log.pred.mu, lambda=exp(log.pred.mu)) #get prediction interval
pred1<-matrix(pred1, nrow=nrow(log.pred.mu))
reg_unbias2<-cbind.data.frame(pred.indices.spl,t(pred1))
reg_unbias_m<-melt(reg_unbias2, id=c('time','country','state'))
reg_unbias_c<-acast(reg_unbias_m, variable~time~country~state)
reg_unbias_c<-reg_unbias_c[,order(as.numeric(dimnames(reg_unbias_c)[[2]])),order(as.numeric(dimnames(reg_unbias_c)[[3]])),order(as.numeric(dimnames(reg_unbias_c)[[4]])), drop=F]
preds.unbias.q<-apply(reg_unbias_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.unbias.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

#fitted values that only account for parameter uncertainty--use this for 2nd stage model
log_reg_mean2<-cbind.data.frame(pred.indices.spl,t(log.pred.mu))
log_reg_mean_m<-melt(log_reg_mean2, id=c('time','country','state'))
log_reg_mean_c<-acast(log_reg_mean_m, variable~time~country~state)
log_reg_mean_c<-log_reg_mean_c[,order(as.numeric(dimnames(log_reg_mean_c)[[2]])),order(as.numeric(dimnames(log_reg_mean_c)[[3]])),order(as.numeric(dimnames(log_reg_mean_c)[[4]])), drop=F]
preds.logregmean.q<-apply(log_reg_mean_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.logregmean.q)[[2]]<- as.numeric(as.character(dimnames(preds.logregmean.q)[[2]]))
preds.logregmean.sd<-apply(log_reg_mean_c,c(2,3,4),sd, na.rm=TRUE)

saveRDS(preds.logregmean.sd, file=paste0(output_directory,"preds.logregmean.sd", ag.select,".rds"), compress=FALSE)
saveRDS(preds.logregmean.q, file=paste0(output_directory,"preds.logregmean.q", ag.select,".rds"), compress=FALSE)
saveRDS(outcome.array.pre, file=paste0(output_directory,"outcome.array.pre", ag.select,".rds"), compress=FALSE)


tiff(paste0(output_directory,'obs.vs.exp.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
matplot(t(preds.unbias.q[,,i,1]), type='l', col='gray', lty=c(2,1,2), bty='l')
  title(countries[i])
  points(outcome.array[i,1,], cex=0.5)
  abline(v=n.times.pre[i], lty=2, col='red')
}
dev.off()
#####################################
#####Pointwise Rate ratio calculation
#####################################
rr<-array(NA, dim=dim(reg_unbias_c))
rd<-array(NA, dim=dim(reg_unbias_c))
for(i in c(1:N.countries)){
rr[,,i,1]<-t(apply(reg_unbias_c[,,i,1],1,function(x) (outcome.array[i,1,]+0.5) / (x+0.5)   )) 
rd[,,i,1]<- t(apply(reg_unbias_c[,,i,1],1,function(x) outcome.array[i,1,] - x   )) 
}



rr.q<-apply(rr,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
rd.q<-apply(rd,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)

log.rr.sd<-apply(log(rr),c(2,3,4),sd,na.rm=TRUE)

dimnames(rr.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))
dimnames(rd.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

log.rr.q<-log(rr.q)

tiff(paste0(output_directory,'rd.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(rd.q[,,i,1]), type='l', col='gray', lty=c(2,1,2), bty='l')
  title(paste0("Rate Difference ",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=sqrt(0), lty=2, col='red')
}
dev.off()


tiff(paste0(output_directory,'rr.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(log.rr.q[,,i,1]), type='l', col='gray', lty=c(2,1,2),ylim=c(-1.0,1.5), bty='l')
  title(paste0("Rate Ratio ",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=0, lty=2, col='red')
}
dev.off()

saveRDS(log.rr.q, file=paste0(output_directory,"log.rr.q_", ag.select,".rds"))
saveRDS(log.rr.sd, file=paste0(output_directory,"log.rr.sd_", ag.select,".rds"))
 }
}

###############NEXT PERFORM THE SMOOTHING META_MODEL

for(ag.select in c('2-59m','2-11m','12-23m','24-59m','<2m')){
  rm(list=ls()[-which(ls() %in% c('ag.select', 'subnat'))]) #for instance 
  if(subnat){
    output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/subnat'   #Directory where results will be saved.
  }else{
    output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/nat'   #Directory where results will be saved.
  }
  if(ag.select=='<2m'){ ag.select<-'u2m'}
  output_directory <- paste(output_directory,  '/',ag.select,'/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
  dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  preds.logregmean.sd<- readRDS( file=paste0(output_directory,"preds.logregmean.sd", ag.select,".rds"))
  preds.logregmean.q<- readRDS( file=paste0(output_directory,"preds.logregmean.q", ag.select,".rds"))
  outcome.array.pre<- readRDS( file=paste0(output_directory,"outcome.array.pre", ag.select,".rds"))
  
  #########################
  #CALL JAGS POOLING MODEL
  #########################
  source('model.pool.R')
  
  
  
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
  tiff(paste0(output_directory,'cp by country.tiff'), width = 7, height = 8, units = "in",res=200)
  cp.plot<-ggplot(data=plot.cp, aes(x=median.cp, y=order2, color=country2)) +
    geom_point(aes(size=inv.var.cp)) +
    scale_size_continuous(range=c(1,15)) +
    theme_bw()+
    guides( size = FALSE)+
    # theme(legend.position = "none")+
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
  print(cp.plot)
  dev.off()
  saveRDS(plot.cp, file=paste0(output_directory,'CP1 pool model.rds'))
  
  
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
  tiff(paste0(output_directory,'country.slope.tiff'), width = 7, height = 8, units = "in",res=200)
  slp.plot<-ggplot(data=plot.slp, aes(x=median.slp, y=order2, color=country2)) +
    geom_point(aes(size=inv.var.slp)) +
    scale_size_continuous(range=c(1,15)) +
    theme_bw()+
    guides( size = FALSE)+
    # theme(legend.position = "none")+
    scale_color_manual(values=c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c','#fdbf6f','#ff7f00','#cab2d6',  '#6a3d9a'))
  print(slp.plot)
  dev.off()
  saveRDS(plot.slp, file=paste0(output_directory,"slope pool model.rds"))
  
  ##########################################
  
  #Relationship between change point location and slope
  tiff(paste0(output_directory,'cp vs slope.tiff'), width = 7, height = 8, units = "in",res=200)
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
  
  tiff(paste0(output_directory,'country.trend.tiff'), width = 7, height = 8, units = "in",res=200)
  par(mfrow=c(5,2), mar=c(4,2,1,1))
  for(i in 1:length(countries)){
    # for(j in 1:N.states[i]){
    plot.data<-t(preds.unbias.q[,,i])
    matplot( ((1:tot_time)-pre.vax.time), plot.data,type='l',yaxt='n', xlim=c(0, max.time.points), xlab='months post-PCV introduction', ylim=c(-0.7,0.4), col='gray', lty=c(2,1,2), bty='l')
    abline(h=0)
    axis(side=2, at=c(-0.7,-0.35,0,0.35,0.7), las=1,labels=round(exp(c(-0.7,-0.35,0,0.35,0.7)),1 ))
    # abline(v=0)
    title(countries[i])
  }
  dev.off()
  saveRDS(preds.unbias.q, file=paste0(output_directory,"reg_mean_with_pooling CP nobias.rds"))
  
}



