##Test model with hospitalization data


library(RCurl)
library(reshape2)
library(lubridate)
library(abind)
setwd('C:/Users/dmw63/Desktop/My documents h/PAHO mortality/paho mrtality jags')
input_directory  <- "C:/Users/dmw63/Your team Dropbox/PAHO mortality/Data/" #Directory or URL page containing input data file
file_name="PAHO all age cuts_SubChapters.csv"
output_directory <- '../JAGS_results'   #Directory where results will be saved.
output_directory <- paste(output_directory,'_', format(Sys.time(), '%Y-%m-%d-%H%M%S'), '/', sep = '')                     #Adds a subfolder to output directory to organize results by date and time run.
data_file <- paste0(input_directory, file_name)
#prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM LOCAL
prelog_data <- read.csv(data_file, check.names = FALSE)# IF IMPORTING FROM URL
prelog_data<-prelog_data[substr(prelog_data$age_group,4,8)=='2-23m',]  #Only <12m
prelog_data<-prelog_data[which(substr(prelog_data$age_group,10,10)!='A'),]  #filter out summary categories
prelog_data<-prelog_data[,c('age_group', 'monthdate','J12_J18_prim','acm_noj_prim' )]
prelog_data$monthdate<-as.Date(prelog_data$monthdate)
prelog_data$country<-substr(prelog_data$age_group,1,2)
################Set Vax intro date for each country
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
######################
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
prelog_data_spl<-prelog_data_spl[,c( 'index','country_index','grp.index',"J12_J18_prim","acm_noj_prim","post.index")]

ds.m<-melt(prelog_data_spl,id.vars=c( 'index','country_index','grp.index'))
outcome.array<- acast(ds.m[ds.m$variable=="J12_J18_prim",], country_index~grp.index~index )   #Outcome array
control1.array<- acast(ds.m[ds.m$variable=="acm_noj_prim",],  country_index~grp.index~index )  #Control array
control1.array<-log(control1.array+0.5) #Log transform
post.index.array <- acast(ds.m[ds.m$variable=="post.index",],  country_index~grp.index~index ) #time index array
post.index.array[post.index.array<0]<-0 #Pre-vax indices are 0
post.index.array<-post.index.array/dim(post.index.array)[3]
  
#Combine covariate array with intercept along 4th dimension
int.array<-array(1, dim(control1.array))
control1.array.int<-abind(int.array,control1.array, along=4 )
control1.array.int2<-abind(control1.array.int,post.index.array, along=4 )

N.countries=dim(control1.array.int2)[1]
n.times=apply(outcome.array,c(1,2), function(x) sum(!is.na(x))) #Vector (or matrix) giving number of non-missing observations per state
# test.array<-post.index.array
# test.array[is.na(test.array)]<-0
# n.times.test<-apply(test.array, c(1,2), function(x) max(x, na.rm=TRUE))
N.states.country<-apply(n.times,1, function(x) sum(x!=0)  )
#which.states<-apply(n.times,1, function(x) paste(which(x!=0), collapse=',')  )

N.preds=3 #N predictors ( intercept, acm_noj, post-vax time trend)
q=3  #N predictors of slopes (if none, set to 1 for int only)
# I_Sigma<-replicate( N.countries, diag(N.preds) )
I_Sigma<-diag(N.preds) 
I_Omega<-diag(q)
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

source('model.R')