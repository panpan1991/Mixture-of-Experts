library(pscl)
library(MASS)
library(pgdraw)

#load dataset
dat=read.table("data4hw6")
Y=dat$V2
#construct design matrix
X=matrix(c(rep(1,length(Y)),dat$V1), ncol = 2)

X.append=matrix(c(rep(1,250), seq(1.0001,2,length.out=250)), ncol = 2)
Y.append=Y[1:250]

X=rbind(X,X.append)
Y=c(Y,Y.append)
plot(X[,2],Y)

n=length(Y)
k=3 #number of components in the mixture of regression
  
#Hyperparameter for sigma's square's prior
alpha=1
lamda=1
#Hyperparameter for deltas' prior
B0=diag(100, nrow=2,ncol=2)

#Starting point
sigma.sq=rep(1,k)
delta=matrix(rep(c(1,1),k),ncol = k)
p=exp(X%*%delta)/(1+exp(X%*%delta))

#Z=rbern(n,p)
Z=t(apply(p, 1, rmultinom, n=1, size=1))


nk=apply(Z, 2, sum)
#Xk[[j]] is partition of design matrix that assign to regression j
Xk=list()
for (j in 1:k) {
  Xk[[j]]=matrix(X[Z[,j]==1,],ncol = 2)
}

beta=matrix(rep(c(1,1),k),ncol = k)

####################Iterations###################################
n.iteration=5000
#containers
sigma.sq.record=matrix(NA, nrow = n.iteration, ncol = k)
beta.record=list()
delta.record=list()


error=0
i=1
for(i in 1:5){

#???
individual.density=dnorm(Y,mean=X%*%beta, 
                         matrix(rep(sigma.sq,n), byrow = TRUE, ncol = 3))


mixture=apply(p*individual.density, 1, sum)
mixture=matrix(rep(mixture,3), ncol = 3)
prob=individual.density/mixture
Z=t(apply(prob, 1, rmultinom, n=1, size=1))

###########################################################
#n1=sum(Z)
nk=apply(Z, 2, sum)

for (j in 1:k) {
  Xk[[j]]=matrix(X[Z[,j]==1,],ncol = 2)
}

for (j in 1:k) {
  if(nk[j]>3){
    D=diag(Z[,j])
    sigma.sq[j]=rigamma(1, 
                      sum(Z[,j])/2+alpha+1, 
                      1/2*t(Y-X%*%beta[,j])%*%D%*%(Y-X%*%beta[,j])+
                        lamda+1/2*t(beta[,j])%*%((1/nk[j])*t(Xk[[j]])%*%Xk[[j]])%*%beta[,j])
    
    #betas
    B=diag(Z[,j]/sigma.sq[j])
    m0=nk[j]*sigma.sq[j]*solve(t(Xk[[j]])%*%Xk[[j]])
    m=solve(t(X)%*%B%*%X+solve(m0))
    beta[,j]=mvrnorm(1,m%*%t(X)%*%B%*%Y,m)
  }else{
    error=error+1
  }
}

##################################################################
#W=rpg(500,1,X%*%delta)
for (j in 1:k) {
  
  C=log(exp(X%*%delta[,-j]))
  Cj=apply(C, 1, sum)
  W=pgdraw(1, X%*%delta[,j]-Cj)
  
  K=Z[,j]-1/2
  omega=diag(W)
  mw=solve(t(X)%*%omega%*%X+solve(B0))%*%(t(X)%*%(K+omega%*%Cj))
  vw=solve(t(X)%*%omega%*%X+solve(B0))
  
  delta[,j]=mvrnorm(1,mw,vw)
}
p=exp(X%*%delta)/(1+exp(X%*%delta))

sigma.sq.record[i,]=sigma.sq
beta.record[[i]]=beta
delta.record[[i]]=delta

print( i)

}


#Trace plots
par(mfrow=c(2,2))
plot(beta1.record[,2], main="Beta1 slope")
plot(beta1.record[,1], main="Beta1 intercept")

plot(beta2.record[,2], main="Beta2 slope")
plot(beta2.record[,1], main="Beta2 intercept")
par(mfrow=c(1,2))
plot(sigma.record[,1], main="Sigma 1")
plot(sigma.record[,2], main="Sigma 2")

par(mfrow=c(1,1))
#Fitted Curve
fitted=(t(p.record[1000:n.iteration,])*X%*%t(beta1.record[1000:n.iteration,])+t(1-p.record[1000:n.iteration,])*X%*%t(beta2.record[1000:n.iteration,]))
fit=apply(fitted, 1, mean)
plot(dat$V1,Y, xlab="X", ylab = "Y", main = "Sample points and fitted line")
lines(dat$V1, fit, col="red", type = "l")

par(new=FALSE)

####################
fitted=(t(p.record[1000:n.iteration,])*X%*%t(beta1.record[1000:n.iteration,])+t(1-p.record[1000:n.iteration,])*X%*%t(beta2.record[1000:n.iteration,]))
cred=apply(fitted, 1, prob=c(0.25,0.975), quantile)

lines(dat$V1,cred[1,], col="green", type = "l")
lines(dat$V1,cred[2,], col="orange", type = "l")

