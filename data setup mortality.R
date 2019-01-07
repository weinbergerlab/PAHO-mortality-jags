##Test model with hospitalization data


library(RCurl)
library(reshape2)
library(lubridate)
library(abind)
library(dummies)

max.time.points=48


input_directory  <- "C:/Users/dmw63/Your team Dropbox/PAHO mortality/Data/" #Directory or URL page containing input data file
file_name="PAHO all age cuts_SubChapters.csv"
output_directory <- 'C:/Users/dmw63/Desktop/My documents h/PAHO mortality/jags cp results'   #Directory where results will be saved.
output_directory <- paste(output_directory,  '/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
data_file <- paste0(input_directory, file_name)
#prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM LOCAL
prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM URL

#Filter to obtain relevant age groups
prelog_data<-prelog_data[substr(prelog_data$age_group,4,8)=='2-23m',]  #Only <12m
prelog_data<-prelog_data[which(substr(prelog_data$age_group,10,10)!='A'),]  #filter out summary categories
prelog_data<-prelog_data[,c('age_group', 'monthdate','J12_J18_prim','acm_noj_prim' )]
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
prelog_data_spl<-prelog_data_spl[,c( 'index','country_index','grp.index',"J12_J18_prim","acm_noj_prim","post.index","monthdate")]

#Reshape data into 3D arrays--separate arrays for outcome, controls, time index
ds.m<-melt(prelog_data_spl,id.vars=c( 'index','country_index','grp.index'))
outcome.array<- acast(ds.m[ds.m$variable=="J12_J18_prim",], country_index~grp.index~index )   #Outcome array
control1.array<- acast(ds.m[ds.m$variable=="acm_noj_prim",],  country_index~grp.index~index )  #Control array
control1.array<-log(control1.array+0.5) #Log transform
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
control1.array.int<-abind(int.array,control1.array, along=4 )
control1.array.int2<-abind(control1.array.int,post.index.array, along=4 )

#Create variables to be used by JAGS
N.countries=dim(control1.array.int2)[1]
n.times=apply(outcome.array,c(1,2), function(x) sum(!is.na(x))) #Vector (or matrix) giving number of non-missing observations per state
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

#Extract first change point time:
beta4<-posterior_samples[[1]][,grep("^beta.*,4]",dimnames(posterior_samples[[1]])[[2]])] #Intercept
cp1<-max.time.points*exp(beta4) #Intercept
beta4.lab<-x.func(dimnames(beta4)[[2]]) 
beta4.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(beta4.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(beta4.lab.spl)<-c('country','state')
#for(i in c(1:ncol(beta4))){hist(beta4[,i])}
quant.cp1<-t(max.time.points*exp(apply(beta4,2,quantile, probs=c(0.025,0.5,0.975))))
var.cp1<-apply(beta4,2,var)
cp1.lab<-x.func(dimnames(quant.cp1)[[1]]) 
cp1.lab.spl<-as.data.frame(matrix(as.numeric(as.character(unlist(strsplit(cp1.lab, ',',fixed=TRUE)))), ncol=3, byrow=TRUE))
names(cp1.lab.spl)<-c('country','state')
par(mfrow=c(2,2))
plot(y=1:nrow(quant.cp1), x=quant.cp1[,'50%'], bty='l')
library(ggplot2)
plot.cp<-cbind.data.frame('strata'=1:nrow(quant.cp1),'median.cp'=quant.cp1[,'50%'],'lcl.cp'=quant.cp1[,'2.5%'],'ucl.cp'=quant.cp1[,'97.5%'], 'inv.var.cp'=1/var.cp1,beta4.lab.spl)
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
#Unbias
reg_unbias<-posterior_samples[[1]][,grep("log_rr_estimate",dimnames(posterior_samples[[1]])[[2]])]
pred1<-reg_unbias
pred.indices<- x.func(dimnames(pred1)[[2]]) 
pred.indices.spl<-matrix(unlist(strsplit(pred.indices, ',',fixed=TRUE)), ncol=3, byrow=TRUE)
pred.indices.spl<-as.data.frame(pred.indices.spl)
names(pred.indices.spl)<- c('country','state','time')
reg_unbias2<-cbind.data.frame(pred.indices.spl,t(reg_unbias))
reg_unbias_m<-melt(reg_unbias2, id=c('time','country','state'))
reg_unbias_c<-acast(reg_unbias_m, variable~time~country~state)
reg_unbias_c<-reg_unbias_c[,order(as.numeric(dimnames(reg_unbias_c)[[2]])),order(as.numeric(dimnames(reg_unbias_c)[[3]])),order(as.numeric(dimnames(reg_unbias_c)[[4]]))]
preds.unbias.q<-apply(reg_unbias_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.unbias.q)[[2]]<- as.numeric(as.character(dimnames(preds.unbias.q)[[2]]))

tiff(paste0(output_directory,'subnat.rr.tiff'), width = 7, height = 8, units = "in",res=200)
par(mfrow=c(5,2), mar=c(4,2,1,1))
for(i in 1:length(countries)){
  # for(j in 1:N.states[i]){
  for(j in c(1:3)){
  plot.data<-t(preds.unbias.q[,,i,j])
  if( sum(plot.data, na.rm=T)>0){
  matplot( post.index.array[i,1,]*dim(post.index.array)[3], plot.data,type='l',yaxt='n',add=(j>1), xlim=c(0.1, max.time.points), xlab='months post-PCV introduction', col='gray', lty=c(2,1,2), bty='l')
  abline(h=0)
  axis(side=2, at=c(-0.7,-0.35,0,0.35,0.7), las=1,labels=round(exp(c(-0.7,-0.35,0,0.35,0.7)),1 ))
  # abline(v=0)
  }
  title(countries[i])
  }
}
dev.off()
saveRDS(preds.unbias.q, file=paste0(output_directory,"reg_mean_with_pooling CP nobias.rds"))


