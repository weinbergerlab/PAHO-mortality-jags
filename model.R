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
for(v in 1:n.pre[i,j]){     
    Y[i,j,v]~dpois(mu[i,j,v]) #Likelihood

    ##################### 
    #MODEL #
    #####################
    log(mu[i,j,v])<-(beta[i,j,1] 
    + beta[i,j,2]*x[i,j,v,2]
    # + beta[i,j,3]*x[i,j,v,3]
    # + beta[i,j,4]*x[i,j,v,4]
    # + beta[i,j,5]*x[i,j,v,5]
    + month.dummy[i,j,v,1] * delta[i,j,1] + month.dummy[i,j,v,2] * delta[i,j,2]
    + month.dummy[i,j,v,3] * delta[i,j,3] + month.dummy[i,j,v,4] * delta[i,j,4]
    + month.dummy[i,j,v,5] * delta[i,j,5] + month.dummy[i,j,v,6] * delta[i,j,6]
    + month.dummy[i,j,v,7] * delta[i,j,7] + month.dummy[i,j,v,8] * delta[i,j,8]
    + month.dummy[i,j,v,9] * delta[i,j,9] + month.dummy[i,j,v,10] * delta[i,j,10]
    + month.dummy[i,j,v,11] * delta[i,j,11]
    +disp[i,j,v]
    )
disp[i,j,v]~dnorm(0, tau.disp[i,j])
}

for(u in 1:n[i,j]){     
  ##################### 
  #Prediction#
  #####################
log.pred.mu[i,j,u]<-(beta[i,j,1] 
+ beta[i,j,2]*x[i,j,u,2]
# + beta[i,j,3]*x[i,j,u,3]
# + beta[i,j,4]*x[i,j,u,4]
# + beta[i,j,5]*x[i,j,u,5]
+ month.dummy[i,j,u,1] * delta[i,j,1] + month.dummy[i,j,u,2] * delta[i,j,2]
+ month.dummy[i,j,u,3] * delta[i,j,3] + month.dummy[i,j,u,4] * delta[i,j,4]
+ month.dummy[i,j,u,5] * delta[i,j,5] + month.dummy[i,j,u,6] * delta[i,j,6]
+ month.dummy[i,j,u,7] * delta[i,j,7] + month.dummy[i,j,u,8] * delta[i,j,8]
+ month.dummy[i,j,u,9] * delta[i,j,9] + month.dummy[i,j,u,10] * delta[i,j,10]
+ month.dummy[i,j,u,11] * delta[i,j,11]
+disp.pred[i,j,u]
)
disp.pred[i,j,u]~dnorm(0, tau.disp[i,j])
}




w_true_var_inv[i,j]<-1/(w_true_sd[i,j]*w_true_sd[i,j])
w_true_sd[i,j] ~ dunif(0, 1000)

tau.disp[i,j]<-1/sd.disp[i,j]^2
sd.disp[i,j]~dunif(0,100)
###########################################################
#Second Stage Statistical Model
###########################################################
beta[i,j, 1:2] ~ dmnorm(lambda[1:2], Omega_inv[1:2, 1:2])
delta[i,j, 1:11] ~ dmnorm(theta[1:11], Omicron_inv[1:11, 1:11]) #Seasonal parameters

}
########################################################
#Third Stage Statistical Model
########################################################
#gamma[i, 1:2] ~ dmnorm(lambda[1:2], Sigma_inv[1:2, 1:2])
}
#######################################################
#Remaining Prior Distributions
#######################################################
Omega_inv[1:2, 1:2] ~ dwish(I_Omega[1:2, 1:2], (2 + 1))
Omega[1:2, 1:2]<-inverse(Omega_inv[1:2, 1:2])

Omicron_inv[1:11, 1:11] ~ dwish(I_Omicron[1:11, 1:11], (11 + 1))
Omicron[1:11, 1:11]<-inverse(Omicron_inv[1:11, 1:11])

for(j in c(1:2)){
lambda[j] ~ dnorm(0, 1e-4)
}
for(k in c(1:11)){
theta[k] ~ dnorm(0, 1e-4)
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
                                 'n' = n.times, 
                                 'n.pre' = n.times.pre, 
                                 'p' = N.preds, #predictors of outcome
                                 'q' = q,  #predictors of slope
                                 'm' = m,
                                 'Y' = outcome.array.pre,
                                 'O' = offset,
                                 'z' =     z   ,
                                 'month.dummy'=month.dummy,
                                 'x' = control1.array.int,
                                 'w' = w,
                                 'I_Sigma' = I_Sigma,
                                 'I_Omega' = I_Omega,
                                 'I_Omicron'=I_Omicron
                                 )) 

update(model_jags, 
       n.iter=20000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c( 'mu','log.pred.mu',"beta", 'lambda'),
                                
                                thin=1,
                                n.iter=50000)

#post1.summary<-summary(posterior_samples)
# plot(posterior_samples, 
#      ask=TRUE)


