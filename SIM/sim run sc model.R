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


input_directory  <- "C:/Users/dmw63/Your team Dropbox/PAHO mortality/Data/" #Directory or URL page containing input data file
file_name="PAHO all age cuts_SubChapters.csv"
data_file <- paste0(input_directory, file_name)
#prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM LOCAL
prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM URL

#Filter to obtain relevant age groups
#prelog_data<-prelog_data[substr(prelog_data$age_group,4,8)==ag.select,]  #Only <12m
prelog_data<-prelog_data[grep(ag.select,prelog_data$age_group ),]
if(ag.select=='2-23m'){ prelog_data<-prelog_data[!grepl('12-23m',prelog_data$age_group ),]}
if(subnat){
  prelog_data<-prelog_data[!grepl('A',prelog_data$age_group ),]  #filter out summary categories
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/subnat'   #Directory where results will be saved.
  
}else{
  prelog_data<-prelog_data[grepl('A',prelog_data$age_group ),]  #filter out summary categories
  output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results/nat'   #Directory where results will be saved.
}
if(ag.select=='<2m'){ ag.select<-'u2m'}
output_directory <- paste(output_directory,  '/',ag.select,'.SIM','/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
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
library(stringr)
prelog_data$interval.date<-NULL
prelog_data$hdi<-word(prelog_data$age_group,start=3)

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
outcome.array<- sim.counts
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

outcome.array.pre<- outcome.array
outcome.array.pre[post.index.array>0]<- NA

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

aux.output1<-list( 'countries'=countries,'month.dummy'=month.dummy,'outcome.array'=outcome.array,'post.index.array'=post.index.array,
                 'n.times.pre'=n.times.pre,'n.times'=n.times,'N.countries'=N.countries,'N.states.country'=N.states.country )
saveRDS(aux.output1, file=paste0(output_directory,"in.data.", ag.select,".rds"), compress=FALSE)


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
rr.dynamic<-array(NA, dim=dim(reg_unbias_c))
rr.0.05<-array(NA, dim=dim(reg_unbias_c))
rd<-array(NA, dim=dim(reg_unbias_c))
scaled.rd<-array(NA, dim=dim(reg_unbias_c))
dynamic.correction<-array(NA, dim=dim(reg_unbias_c[,,,]))

comb.pred.array<- abind(reg_unbias_c, log_reg_mean_c , along=5)
outcome.array.reorder<-aperm(outcome.array, c(3,1,2)) #Change dimension order of array, and add dimension to match comb.pred.array

#str(log_reg_mean_c[k,,,])
#mean.counts<- apply(exp(log_reg_mean_c), c(3), mean, na.rm=T) #Mean counts by country

#Continuity corrections that minimize bias basd on a simulation https://github.com/weinbergerlab/continuity-correction
dynamic.correction[exp(log_reg_mean_c) >=4]<-0.5
dynamic.correction[exp(log_reg_mean_c) <4 & exp(log_reg_mean_c)>=3 ]<-0.4 
dynamic.correction[exp(log_reg_mean_c) <3 & exp(log_reg_mean_c)>=2 ]<-0.3 
dynamic.correction[exp(log_reg_mean_c) <2 & exp(log_reg_mean_c)>=1 ]<-0.2 
dynamic.correction[exp(log_reg_mean_c) <1 & exp(log_reg_mean_c)>=0.5 ]<-0.08 
dynamic.correction[exp(log_reg_mean_c) <0.5 & exp(log_reg_mean_c)>=0.4 ]<-0.04  
dynamic.correction[exp(log_reg_mean_c) <0.4 & exp(log_reg_mean_c)>=0.3 ]<-0.02  
dynamic.correction[exp(log_reg_mean_c) <0.3 & exp(log_reg_mean_c)>=0.2 ]<-0.01 
dynamic.correction[exp(log_reg_mean_c) <0.2 & exp(log_reg_mean_c)>=0.1 ]<-0.001  
dynamic.correction[exp(log_reg_mean_c) <0.1 & exp(log_reg_mean_c)>=0.05 ]<-0.00001
dynamic.correction[exp(log_reg_mean_c) <0.05 & exp(log_reg_mean_c)>=0.01 ]<-0.00001

#Noteif using multiple states, need to modify this
for(k in 1:dim(reg_unbias_c)[1]){
  rr[k,,,]<-(outcome.array.reorder[,,1]+0.5)/(reg_unbias_c[k,,,]+0.5)
  rr.0.05[k,,,]<-(outcome.array.reorder[,,1]+0.05)/(reg_unbias_c[k,,,]+0.05)
  rd[k,,,]<-(outcome.array.reorder[,,1])-(reg_unbias_c[k,,,])
  scaled.rd[k,,,]<- 1- ( reg_unbias_c[k,,,]-outcome.array.reorder[,,1])/exp(log_reg_mean_c[k,,,]) # 1-(exp-obs)/mu1
  rr.dynamic[k,,,]<-(outcome.array.reorder[,,1]+dynamic.correction[k,,])/(reg_unbias_c[k,,,]+dynamic.correction[k,,])
}

 
rr.q<-apply(rr,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
rr.0.05.q<-apply(rr.0.05,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
rr.dynamic.q<-apply(rr.dynamic,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)

scaled.rd.q<-apply(scaled.rd,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
rd.q<-apply(rd,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
par(mfrow=c(1,1))
plot(rr.q[2,,,],scaled.rd.q[2,,,], main='Median RR estimates w/ and w/o continuity correction', bty='l')

log.rr.sd<-apply(log(rr),c(2,3,4),sd,na.rm=TRUE)
log.rr.0.05.sd<-apply(log(rr.0.05.q),c(2,3,4),sd,na.rm=TRUE)
log.rr.dynamic.sd<-apply(log(rr.dynamic.q),c(2,3,4),sd,na.rm=TRUE)

scaled.rd.sd<-apply(scaled.rd,c(2,3,4),sd,na.rm=TRUE)

dimnames(rr.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))
dimnames(rr.0.05.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))
dimnames(rd.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))
dimnames(rr.dynamic)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

log.rr.q<-log(rr.q)
log.rr.0.05.q<-log(rr.0.05.q)
log.rr.dynamic.q<-log(rr.dynamic.q)

tiff(paste0(output_directory,'rd.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(rd.q[,,i,1]), type='l', col='gray', lty=c(2,1,2), bty='l')
  title(paste0("Rate Difference ",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=sqrt(0), lty=2, col='red')
}
dev.off()

tiff(paste0(output_directory,'scaled.rd.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(scaled.rd.q[,,i,1]), type='l', col='gray', lty=c(2,1,2), ylim=c(-1,3) ,bty='l')
  title(paste0("Scaled Rate Difference ",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=1, lty=2, col='red')
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

tiff(paste0(output_directory,'rr.0.05.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(log.rr.0.05.q[,,i,1]), type='l', col='gray', lty=c(2,1,2),ylim=c(-1.0,1.5), bty='l')
  title(paste0("Rate Ratio ",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=0, lty=2, col='red')
}
dev.off()

tiff(paste0(output_directory,'rr.dynamic.',ag.select,'.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(2,2,1,1))
for(i in 1: length(countries)){
  matplot(t(log.rr.dynamic.q[,,i,1]), type='l', col='gray', lty=c(2,1,2),ylim=c(-1.0,1.5), bty='l')
  title(paste0("Rate Ratio With Dynamic Cont Correct.",countries[i]))
  abline(v=n.times.pre[i], lty=2, col='red')
  abline(h=0, lty=2, col='red')
}
dev.off()

saveRDS(log.rr.q, file=paste0(output_directory,"log.rr.q_", ag.select,".rds"))
saveRDS(scaled.rd.q, file=paste0(output_directory,"scaled.rd.q_", ag.select,".rds"))
saveRDS(scaled.rd.sd, file=paste0(output_directory,"scaled.rd.sd_", ag.select,".rds"))
saveRDS(log.rr.sd, file=paste0(output_directory,"log.rr.sd_", ag.select,".rds"))

saveRDS(log.rr.0.05.q, file=paste0(output_directory,"log.rr.0.05.q_", ag.select,".rds"))
saveRDS(log.rr.0.05.sd, file=paste0(output_directory,"log.rr.0.05.sd_", ag.select,".rds"))

saveRDS(log.rr.dynamic.q, file=paste0(output_directory,"log.rr.0.05.q_", ag.select,".rds"))
saveRDS(log.rr.dynamic.sd, file=paste0(output_directory,"log.rr.0.05.sd_", ag.select,".rds"))

