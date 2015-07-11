set.seed(123456)

n <- 10^2

dt <- 0.1
erra <- 0.2
errm <- 10

x0 <- (dt^2)/2

ex <- (x0*erra)^2
ev <- (dt*erra)^2




C <- matrix(data=c(1,0),nrow=1,ncol=2)
acel<-vector(length=n)
Sz <- matrix(data=c(errm),nrow=1,ncol=1)
#Sk <- matrix(data=c(ex^2,ex*ev,ex*ev,ev^2),nrow=2,ncol=2)
Sk <- matrix(data=c(ex^2,0,0,ev^2),nrow=2,ncol=2)


P <- array(dim=c(n+1,2,2))
P[1,,] <- matrix(data=c(20,20,20,20),nrow=2,ncol=2)
est_Psi_k <- matrix(data=c(0,0),nrow=2,ncol=1)



S <- vector(length=n)

Pose <- vector(length = n)

err<-vector(length = n)

S[1]<-1/2+rnorm(1,0,erra)

for(k in 1:n){
  
  if(is.nan(P[k,1,1])){
    break
  }
  
  deltaT <- dt*k
  
  A <- matrix(data=c(1,0,deltaT,1),nrow=2,ncol=2)
  B <- matrix(data=c(deltaT^2/2,deltaT),nrow=2,ncol=1)
    
  
  acel[k] <- rnorm(1,1,erra)
  
  S[k+1] <-(acel[k]*deltaT^2)/2
  
  ex1<- S[k+1]
  
  err[k] <- ex1-(C%*%est_Psi_k)
  
  Kk <- A%*%P[k,,]%*%t(C)%*%solve(C%*%P[k,,]%*%t(C)+Sz)
  
  est_Psi_k1 <- (A%*%est_Psi_k + B%*%acel[k]) + Kk%*%err[k]

  P[k+1,,] <- A%*%P[k,,]%*%t(A) + Sk - 
    A%*%P[k,,]%*%t(C)%*%solve(Sz)%*%C%*%P[k,,]%*%t(A)
  
   
  est_Psi_k <- est_Psi_k1
  Pose[k] <- est_Psi_k1[1]
}
k
plot(Pose,pch=".",xlim=c(0,n),ylim=c(0,max(S)))
lines(Pose,col=1)
lines(S,col=2)
P[1,,]
P[k-1,,]

