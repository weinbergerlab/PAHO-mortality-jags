##############
#Packages
##############
library(rjags)
###################################################
##### JAGS (Just Another Gibbs Sampler) model #####
###################################################
model_string<-"
model{
for(i in 1:c){
for(j in 1:s[i]){
##################
#Observation model
##################
for(v in 1:n[i,j]){     
    Y[i,j,v]~dpois(mu[i,j,v]) #Likelihood

    ##################### 
    #CHANGE POINT MODEL #
    #####################
    log(mu[i,j,v])<-(beta[i,j,1] 
    + step(time.index[i,j,v] - cp1[i,j])*(1 - step(time.index[i,j,v] - cp2[i,j]))*beta[i,j,2]*(time.index[i,j,v] - cp1[i,j]) 
    + step(time.index[i,j,v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]) + beta[i,j,3]*x[i,j,v,2]
    +disp[i,j,v]
)
    reg_unbias[i,j,v]<-(
    step(time.index[i,j,v] - cp1[i,j])*(1 - step(time.index[i,j,v] - cp2[i,j]))*beta[i,j,2]*(time.index[i,j,v] - cp1[i,j]) 
    + step(time.index[i,j,v] - cp2[i,j])*beta[i,j,2]*(cp2[i,j] - cp1[i,j]) + beta[i,j,3]*x[i,j,v,2]
    +disp[i,j,v]
    )

disp[i,j,v]~dnorm(0, tau.disp[i,j])
}
for(k1 in 1:n[i,j]){
for(k2 in 1:n[i,j]){
w_true_cov_inv[i,j,k1,k2]<-ifelse(k1==k2, w_true_var_inv[i,j], 0)
}
}


w_true_var_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
w_true_sd[i,j] ~ dunif(0, 1000)
cp1[i,j]<-exp(beta[i,j,4])
cp2.add[i,j]<-exp(beta[i,j,5])
tau.disp[i,j]<-1/sd.disp[i,j]^2
sd.disp[i,j]~dunif(0,100)
cp2[i,j]<-cp1[i,j] +cp2.add[i,j]  + 1/max.time.points   #ensure Cp2 is at least 1 unit after CP1
###########################################################
#Second Stage Statistical Model
###########################################################
beta[i,j, 1:5] ~ dmnorm(lambda[1:5], Omega_inv[1:5, 1:5])
}
########################################################
#Third Stage Statistical Model
########################################################
#gamma[i, 1:5] ~ dmnorm(lambda[1:5], Sigma_inv[1:5, 1:5])
}
#######################################################
#Remaining Prior Distributions
#######################################################
Omega_inv[1:5, 1:5] ~ dwish(I_Omega[1:5, 1:5], (5 + 1))
Omega[1:5, 1:5]<-inverse(Omega_inv[1:5, 1:5])

for(j in c(1:5)){
lambda[j] ~ dnorm(0, 1e-4)
}

}
"

##############################################################
#Model Fitting
##############################################################
model_jags<-jags.model(textConnection(model_string),
                       data=list('c' = N.countries, 
                                 's' = N.states.country, 
                                 'time.index'=post.index.array,
                                 'max.time.points'= max.time.points , 
                                 'n' = n.times, 
                                 'p' = N.preds, #predictors of outcome
                                 'q' = q,  #predictors of slope
                                 'm' = m,
                                 'Y' = outcome.array,
                                 'O' = offset,
                                 'z' =     z   ,
                                 'x' = control1.array.int,
                                 'w' = w,
                                 'I_Sigma' = I_Sigma,
                                 'I_Omega' = I_Omega)) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("reg_mean",'reg_unbias' 
                                                 ,"beta",'lambda'),
                                
                                thin=1,
                                n.iter=20000)

#post1.summary<-summary(posterior_samples)
# plot(posterior_samples, 
#      ask=TRUE)


