#SETUP Data for simulation study

#for(ag.select in c('2-59m','2-11m','12-23m','24-59m','<2m')){
#for(ag.select in c('2-11m','12-23m','24-59m','<2m', '2-23m')){
 ag.select<-'2-59m'
 subnat='F'
  
    print(ag.select)
    print(subnat)
rm(list=ls()[-which(ls() %in% c('ag.select', 'subnat'))]) #for instance 
set.seed(123)
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
reg_mean2<-cbind.data.frame(pred.indices.spl,t(log.pred.mu))
reg_mean_m<-melt(reg_mean2, id=c('time','country','state'))
reg_mean_c<-acast(reg_mean_m, variable~time~country~state)
preds.q<-apply(reg_mean_c,c(2,3,4),quantile, probs=c(0.025,0.5,0.975),na.rm=TRUE)
dimnames(preds.q)[[2]]<- as.numeric(as.character(dimnames(preds.q)[[2]]))
preds.q<-preds.q[,order(as.numeric(dimnames(preds.q)[[2]])),order(as.numeric(dimnames(preds.q)[[3]])),]


#For each country, have decline, that varies in magnitude
dur.decline<-24 #how long does it take to achieve max effect from change point 1 to change point 2
cp1<-matrix(c(12,20,6,8,15,14,15,17,13,10), nrow=10, ncol=1)/max.time.points
cp2<-matrix(c(12,20,6,8,15,14,15,17,13,10)+dur.decline, nrow=10, ncol=1)/max.time.points
max.rr<-0.8
beta.fix<- log(max.rr)/(dur.decline/max.time.points)

step1<-function(x){
  x>=0
}
vax.effect<- array(NA, dim=dim(post.index.array))
for(i in 1:length(countries)){
  j=1
  for(v in 1: n.times[i,j]){ 
    (vax.effect[i,j,v]<-step1(post.index.array[i,j,v] - cp1[i,j])*(1 - step1(post.index.array[i,j,v] - cp2[i,j]))*beta.fix*(post.index.array[i,j,v] - cp1[i,j]) 
    + step1(post.index.array[i,j,v] - cp2[i,j])*beta.fix*(cp2[i,j] - cp1[i,j]) )
  }
}
matplot(t(vax.effect[,1,]))

log.pred.mu.vax<- preds.q['50%',,] +  t(vax.effect[,1,])
pred.count.vax<-rpois(n=length(log.pred.mu.vax),lambda=exp(log.pred.mu.vax)) #gives error when try to simulate from NA is OK
pred.count.vax<-array(pred.count.vax, dim=dim(log.pred.mu.vax))

#Now need to combine with the real observed data---for pre vax use real data, for post.vax, use 
post.index.array.t<-t(post.index.array[,1,])
outcome.array.t<-t(outcome.array[,1,])
post.index.array.t[is.na(post.index.array.t)]<- -99
combine.counts<-matrix(NA, nrow=nrow(pred.count.vax), ncol=ncol(pred.count.vax))
combine.counts[post.index.array.t==0]<- outcome.array.t[post.index.array.t==0]
combine.counts[post.index.array.t>0]<- pred.count.vax[post.index.array.t>0]
  
matplot(combine.counts, type='l')


combine.counts.array<-array(NA, dim=dim(outcome.array))
combine.counts.array[,1,]<-t(combine.counts)

saveRDS(combine.counts.array, file=paste0(output_directory,"sim.counts.", ag.select,".rds"))



#ALSO create a dataset that can be used in regular SC analysis

#combine.counts
