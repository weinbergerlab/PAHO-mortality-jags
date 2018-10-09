##############
#Packages
##############
library(rjags)

#####################################################################################
#Statistical Model
#####################################################################################
model_string<-"

model{

for(i in 1:c){
for(j in 1:s[i]){
for(t in 1:n[i,j]){

############################################################################
#First Stage Statistical Model
############################################################################
Y[i,j,t] ~ dpois(lambda[i,j,t])
log(lambda[i,j,t]) <- O[i,j,t] + x[i,j,t, 1:p]%*%beta[i,j, 1:p] + phi[i,j,t]
phi[i,j,t] ~ dnorm(0, sigma2_phi_inv)
}

##############################################################
#Second Stage Statistical Model
##############################################################
beta[i,j, 1:p] ~ dmnorm(mu1[i,j, 1:p], Sigma_inv[i, 1:p, 1:p])
for(k in 1:p){
mu1[i,j,k] <- z[i,j, 1:q]%*%gamma[i,k, 1:q]
}
}
for(k in 1:p){

###############################################################
#Third Stage Statistical Model
###############################################################
gamma[i,k, 1:q] ~ dmnorm(mu2[i,k, 1:q], Omega_inv[k, 1:q, 1:q])
for(l in 1:q){
mu2[i,k,l] <- w[i, 1:m]%*%theta[k,l, 1:m]
}
}
Sigma_inv[i, 1:p, 1:p] ~ dwish(I_Sigma[1:p, 1:p], (p + 1))
Sigma[i, 1:p, 1:p] <- inverse(Sigma_inv[i, 1:p, 1:p])
}

#############################################################
#Remaining Prior Distributions
#############################################################
sigma2_phi_inv <- 1/(sigma_phi*sigma_phi)
sigma_phi ~ dunif(0, 1000)

for(k in 1:p){
Omega_inv[k, 1:q, 1:q] ~ dwish(I_Omega[1:q, 1:q], (q + 1))
Omega[k, 1:q, 1:q] <- inverse(Omega_inv[k, 1:q, 1:q])
for(l in 1:q){
for(r in 1:m){
theta[k,l,r] ~ dnorm(0, 0.0001)
}
}
}

}
"

##############################################################
#Model Fitting
##############################################################
model_jags<-jags.model(textConnection(model_string),
                       data=list('c' = N.countries, 
                                 's' = N.states.country, 
                                 'n' = n.times, 
                                 'p' = N.preds, #predictors of outcome
                                 'q' = q,  #predictors of slope
                                 'm' = m,
                                 'Y' = outcome.array,
                                 'O' = offset,
                                 'z' =     z   ,
                                 'x' = control1.array.int2,
                                 'w' = w,
                                 'I_Sigma' = I_Sigma,
                                 'I_Omega' = I_Omega)) 

update(model_jags, 
       n.iter=10000) 

posterior_samples<-coda.samples(model_jags, 
                                variable.names=c("sigma_phi",
                                                 "beta",
                                                 "gamma","theta"),
                                thin=1,
                                n.iter=10000)

#post1.summary<-summary(posterior_samples)
# plot(posterior_samples, 
#      ask=TRUE)


