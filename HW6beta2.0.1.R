library(pscl)
library(MASS)
library(pgdraw)
library(matrixcalc)
library(ggcorrplot)
#load dataset
#dat=read.table("data4hw6")
dat=read.table("h3data")
Y=dat$y
#construct design matrix
X=matrix(c(rep(1,length(Y)),dat$x), ncol = 2)
n=length(Y)

# #Three pieces
# X.append=matrix(c(rep(1,250), seq(1.0001,2,length.out=250)), ncol = 2)
# Y.append=Y[1:250]
# X=rbind(X,X.append)
# Y=c(Y,Y.append)
# 
# #four pieces
# X.append=matrix(c(rep(1,250), seq(2.0001,3,length.out=250)), ncol = 2)
# Y.append=Y[251:500]
# X=rbind(X,X.append)
# Y=c(Y,Y.append)
#n=length(Y)

#Histogram of data
plot(X[,2],Y)


k=10 #number of components in the mixture of regression
  
#Hyperparameter for sigma's square's prior
alpha=1
lamda=1
#Hyperparameter for deltas' prior
B0=diag(100, nrow=2,ncol=2)

#Starting point
sigma.sq=rep(1,k)
beta=matrix(rep(c(1,1),k),ncol = k)
delta=matrix(rep(c(1,1),k),ncol = k)

#This is an n*k matrix
p=exp(X%*%delta)/apply(exp(X%*%delta), 1, sum)
#initial indicator matrix
Z=t(apply(p, 1, rmultinom, n=1, size=1))
Z.record=list()

nk=apply(Z, 2, sum)
#Xk[[j]] is partition of design matrix that assign to regression j
Xk=list()
for (j in 1:k) {
  Xk[[j]]=X[Z[,j]==1,]
}

####################Iterations###################################
n.iteration=1500
#containers
sigma.sq.record=matrix(NA, nrow = n.iteration, ncol = k)
beta.record=list()
delta.record=list()


 error=0
# i=1
for(i in 1:n.iteration){

#???
individual.density=dnorm(Y,mean=X%*%beta,
                         sd=matrix(rep(sqrt(sigma.sq),n), byrow = TRUE, ncol = k))


# for (ii in 1:n) {
#   dnorm(Y[ii],mean=X[ii,]%*%beta,  sd=sqrt(sigma.sq))
# }



mixture=apply(p*individual.density, 1, sum)
#mixture=matrix(rep(mixture,3), ncol = 3)
prob=p*individual.density/mixture
Z=t(apply(prob, 1, rmultinom, n=1, size=1))
Z.record[[i]]=Z
###########################################################
nk=apply(Z, 2, sum)

for (j in 1:k) {
  Xk[[j]]=X[Z[,j]==1,]
}

for (j in 1:k) {
  # Since we have two parameters beta_0 and beta_1 to estimate, 
  # we need two observations being allocated to each expert at least
  if(nk[j]>=2){

    #sigmas
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
  
  Cj=log(apply(exp(X%*%delta[,-j]), 1, sum))
  
  #W=pgdraw(1, X%*%delta[,j]-Cj)
  W=pgdraw(1, X%*%delta[,j]-Cj)
  K=Z[,j]-1/2
  omega=diag(W)
  mw=solve(t(X)%*%omega%*%X+solve(B0))%*%(t(X)%*%(K+omega%*%Cj))
  vw=solve(t(X)%*%omega%*%X+solve(B0))
  
  delta[,j]=mvrnorm(1,mw,vw)
}
p=exp(X%*%delta)/apply(exp(X%*%delta), 1, sum)

sigma.sq.record[i,]=sigma.sq
beta.record[[i]]=beta
delta.record[[i]]=delta

print( i)

}


#Trace plots
par(mfrow=c(3,2))
beta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta1.record[i,]=beta.record[[i]][,1]
}

beta2.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta2.record[i,]=beta.record[[i]][,2]
}

beta3.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  beta3.record[i,]=beta.record[[i]][,3]
}


plot(beta1.record[,2], main="Beta1 slope")
plot(beta1.record[,1], main="Beta1 intercept")

plot(beta2.record[,2], main="Beta2 slope")
plot(beta2.record[,1], main="Beta2 intercept")

plot(beta3.record[,2], main="Beta2 slope")
plot(beta3.record[,1], main="Beta2 intercept")

delta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  delta1.record[i,]=delta.record[[i]][,1]
}

plot(delta1.record[,2], main="delta1 slope")
plot(delta1.record[,1], main="delta1 intercept")

delta1.record=matrix(rep(0,2*n.iteration), nrow = n.iteration)
for (i in 1:n.iteration) {
  delta1.record[i,]=delta.record[[i]][,1]
}


par(mfrow=c(1,3))
plot(sigma.sq.record[,1], main="Sigma 1")
plot(sigma.sq.record[,2], main="Sigma 2")
plot(sigma.sq.record[,3], main="Sigma 2")

par(mfrow=c(1,1))

#Fitted Curve
p.record=list()
for (i in 1:n.iteration) {
  p.record[[i]]=exp(X%*%delta.record[[i]])/apply(exp(X%*%delta.record[[i]]), 1, sum)
}
fitted=list()
for (i in 1:n.iteration) {
  fitted[[i]]=apply(p.record[[i]]*(X%*%beta.record[[i]]), 1,sum)
}


summ=rep(0,n)
for (i in 1001:n.iteration) {
  summ=summ+fitted[[i]]
}

fit=summ/(n.iteration-1000)

plot(X[,2],Y, xlab="X", ylab = "Y", main = "Sample points and fitted line")
lines(X[,2], fit, col="red", type = "l",lwd=3)
lines(X[,2],h3(X[,2]), col="green", type = 'l', lwd=3)

# Add a legend
legend(0.8, 3, legend=c("fitted", "h3"),
       col=c("red", "green"), lty=1:1, lwd=3:3, cex=1.5)



#Average Sigma over H
Sigma.sum=diag(n, x=0)
i=1
error=0
for (i in 1:1500) {
Z=Z.record[[i]]
nk=apply(Z, 2, sum)
if(min(nk)>3){
  
  for (j in 1:k) {
    Xk[[j]]=X[Z[,j]==1,]
  }
  H <- matrix(apply(Z, 2, function(x) rep(x, 2)), ncol = ncol(Z)*2)
  X2<- matrix(rep(X,k), ncol = 2*k)
  
  B=hadamard.prod(H, X2)
  Sigma.beta=diag(2*k, x=0.0)
  for (j in 1:k) {
    Sigma.beta[(2*j-1):(2*j), (2*j-1):(2*j)]=nk[j]*solve(t(Xk[[j]])%*%Xk[[j]])
  }
  I=diag(n)
  Sigma=B%*%Sigma.beta%*%t(B)
  Sigma.sum=Sigma.sum+Sigma
  
}
else
{error=error+1}

}

ggcorrplot(Sigma.sum)




