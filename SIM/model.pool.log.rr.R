###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
library(rjags)

model_string<-"
model{
for(i in 1:n.countries){
for(j in 1:n.states[i]){
for(v in 1:ts.length[i,j]){       
##################
#Observation model
##################
w_hat[v,i,j] ~ dnorm(w_true[i,j, v], w_hat_prec[ v, i,j])

#################################
#Model of 'true' time series data
#################################
w_true[i,j, v] ~ dnorm(reg_mean[i,j,v], w_true_prec_inv[i,j]) 
##################### 
#CHANGE POINT MODEL #
#####################
reg_mean[i,j,v]<-(beta[i,j,1] +
step(time.index[i,j,v] - cp1[i,j])*(1 - step(time.index[i,j,v] - cp2[i,j]))*beta[i,j,2]*(time.index[i,j,v] - cp1[i,j]) 
+step(time.index[i,j,v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]))

reg_unbias[i,j,v]<-(step(time.index[i,j,v] - cp1[i,j])*(1 - step(time.index[i,j,v] - cp2[i,j]))*beta[i,j,2]*(time.index[i,j,v] - cp1[i,j]) 
+step(time.index[i,j,v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]))
}

w_true_prec_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
w_true_sd[i,j] ~ dunif(0, 100)
cp1[i,j]<-exp(beta[i,j,3])
cp2.add[i,j]<-exp(beta[i,j,4])
cp2[i,j]<-cp1[i,j] +cp2.add[i,j]  + t.interval   #ensure Cp2 is at least 1 unit after CP1
###########################################################
#Second Stage Statistical Model
###########################################################
beta[i,j, 1:4] ~ dmnorm(lambda[1:4], Omega_inv[1:4, 1:4])
}
########################################################
#Third Stage Statistical Model
########################################################
#gamma[i, 1:4] ~ dmnorm(lambda[1:4], Sigma_inv[1:4, 1:4])
}
#######################################################
#Remaining Prior Distributions
#######################################################
Omega_inv[1:4, 1:4] ~ dwish(I_Omega[1:4, 1:4], (4 + 1))
Omega[1:4, 1:4]<-inverse(Omega_inv[1:4, 1:4])
# Sigma_inv[1:4, 1:4] ~ dwish(I_Sigma[1:4, 1:4], (4 + 1))
# Sigma[1:4, 1:4]<-inverse(Sigma_inv[1:4, 1:4])
for(j in c(1:4)){
lambda[j] ~ dnorm(0, 1e-4)
}
}
"
#Model Organization
model_jags<-jags.model(textConnection(model_string),
                       data=list('n.countries' = N.countries, 
                                 'n.states' = N.states, 
                                 'w_hat' = log.rr.med,
                                 'w_hat_prec' = log.rr.prec, 
                                 'ts.length' = ts.length_mat,
                                 't.interval'= min(time.index[time.index>0], na.rm=T),
                                 'I_Omega'= I_Omega,
                                 'max.time.points'=max.time.points,
                                 'time.index'=time.index), n.chains=2, n.adapt=2000) 

#Posterior Sampling
update(model_jags, n.iter=10000)  

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean",'reg_unbias' ,'cp1','cp2',"beta",'lambda'),
                                thin=10,
                                n.iter=20000)