library(LaplacesDemon)
library(pscl)
library(MASS)
library(BayesLogit)

#load dataset
dat=read.table("data4hw6")
Y=dat$V2
#construct design matrix
X=matrix(c(rep(1,length(Y)),dat$V1), ncol = 2)
n=length(Y)

  
#Starting point
alpha=1
lamda=1
sigma1.sq=1
sigma2.sq=1

n1=sum(Z)
X1=matrix(X[Z==1,],ncol = 2)
beta1=c(1,1)

n2=n-sum(Z)
X2=matrix(X[Z==0,],ncol = 2)
beta2=c(1,1)

B0=matrix(c(100,0,0,100),ncol = 2)
delta=c(1,1)
p=exp(X%*%delta)/(1+exp(X%*%delta))

Z=rbern(n,p)

####################Iterations###################################
n.iteration=10000
#containers
sigma.record=matrix(NA, nrow = n.iteration, ncol = 2)
beta1.record=matrix(NA, nrow = n.iteration, ncol = 2)
beta2.record=matrix(NA, nrow = n.iteration, ncol = 2)
p.record=matrix(NA, nrow = n.iteration, ncol = 500)
error=0
for(i in 1:n.iteration){
  

p=(dnorm(Y,mean=X%*%beta1, sigma1.sq)*p)/
  (p*dnorm(Y,mean=X%*%beta1, sigma1.sq)+(1-p)*dnorm(Y,mean=X%*%beta2, sigma2.sq))
Z=rbern(n,p)


###########################################################
n1=sum(Z)

if(n1>3){
  X1=matrix(X[Z==1,],ncol = 2)
  D=diag(Z)
  sigma1.sq=rigamma(1, 
                    sum(Z)/2+alpha+1, 
                    1/2*t(Y-X%*%beta1)%*%D%*%(Y-X%*%beta1)+
                      lamda+1/2*t(beta1)%*%((1/n1)*t(X1)%*%X1)%*%beta1)
  
  #betas
  B=diag(Z/sigma1.sq)
  m0=n1*sigma1.sq*solve(t(X1)%*%X1)
  m=solve(t(X)%*%B%*%X+solve(m0))
  beta1=mvrnorm(1,m%*%t(X)%*%B%*%Y,m)
}
else{
  error=error+1
}

n2=n-sum(Z)

if(n2>3){
  X2=matrix(X[Z==0,],ncol = 2)
  Din=diag(rep(1,n))-diag(Z)

  sigma2.sq=rigamma(1, 
                    (n-sum(Z))/2+alpha+1, 
                    1/2*t(Y-X%*%beta2)%*%Din%*%(Y-X%*%beta2)+
                      lamda+1/2*beta2%*%((1/n2)*t(X2)%*%X2)%*%beta2)
  
  #betas
  B=Din/sigma2.sq
  m0=n2*sigma2.sq*solve(t(X2)%*%X2)
  m=solve(t(X)%*%B%*%X+solve(m0))
  beta2=mvrnorm(1,m%*%t(X)%*%B%*%Y,m)
}
else{
  error=error+1
}

##################################################################
W=rpg(500,1,X%*%delta)
K=Z-1/2
omega=diag(W)
mw=solve(t(X)%*%omega%*%X+solve(B0))%*%(t(X)%*%K)
vw=solve(t(X)%*%omega%*%X+solve(B0))
delta=mvrnorm(1,mw,vw)
p=exp(X%*%delta)/(1+exp(X%*%delta))

sigma.record[i,]=c(sigma1.sq,sigma2.sq)
beta1.record[i,]=beta1
beta2.record[i,]=beta2
p.record[i,]=p

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

