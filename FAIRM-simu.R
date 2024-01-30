library(mvtnorm)
library(glmnet)
library(CVXR)
##data generating process for Section 5.1
Data.gen <- function(p, s, h, K.tr, K.te, delta, n.vec, n.test, setting){
  #beta0 <- c(rep(0.5, s), rep(0, p - s))
  if(setting==1){ ##diagonal cov
    beta.S<-seq(0.8,-0.55, length.out=s)
    Sig0=diag(1,p)
  }else if(setting==2){ #equi-cor
    beta.S<-runif(s,-1,1)
    Sig0=diag(0.7,p)+0.3
  }
  X<-matrix(0,ncol=p,nrow=sum(n.vec))
  y<-rep(0, sum(n.vec))
  S=1:s
  Nu<- sample((s+1):p, h, replace=F) #spurious covariates
  tn.vec<-rep(n.test, K.te+K.tr)
  X.test<-NULL
  y.test<-NULL
  for (k in 1:K.tr) {
    #gam.vec<-rnorm(h,0,delta)
    gam.vec<-rnorm(h,-delta,delta)
    ind.k<-ind.set(n.vec,k)
    X[ind.k,] <- rmvnorm(n.vec[k], rep(0, p), Sig0)
    y[ind.k] <-  X[ind.k,S] %*% beta.S+rnorm(n.vec[k]) 
    X[ind.k,Nu]<- y[ind.k]%*%t(gam.vec)+X[ind.k,Nu]
    X.test.k<-rmvnorm(n.test, rep(0, p), Sig0) # test data for the training domains
    y.test.k <-  X.test.k[,S] %*% beta.S+rnorm(n.test) 
    X.test.k[,Nu] <- y.test.k%*%t(gam.vec)+X.test.k[,Nu]
    X.test<-rbind(X.test, X.test.k)
    y.test<-c(y.test, y.test.k)
  }
  
  beta<-rep(0,p)
  beta[S]<-beta.S
  #generate a test model

  for (k in 1:K.te) {
    gam.vec<-rnorm(h,0,delta)
    #gam.vec<-rnorm(h,delta/2,delta)
    ind.k<-ind.set(tn.vec,k)
    X.test.k <- rmvnorm(n.test, rep(0, p), Sig0)
    y.test.k <-  X.test.k[,S] %*% beta.S+rnorm(tn.vec[k]) 
    #X.test.k[,Nu]<- y.test.k%*%t(gam.vec)+X.test.k[,S]%*%matrix(rnorm(s*h,0,delta), nrow=s, ncol=h)+X.test.k[,Nu]
    X.test.k[,Nu]<- y.test.k%*%t(gam.vec)+X.test.k[,Nu]
    X.test<-rbind(X.test, X.test.k)
    y.test<-c(y.test, y.test.k)
  }
  return(list(S=S, Sv=Nu,X.tr=X, y.tr=y, beta=beta, y.test=y.test, X.test=X.test))
}

ind.set <- function(n.vec, k.vec){
  ind.re <- NULL
  for(k in k.vec){
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
  }
  ind.re
}

Maximin<-function(beta.mat,Sigma){
  
  L=ncol(beta.mat)
  Gamma.positive <- diag(0,L)
  for(i in 1:L){
    for(j in 1:L){
      Gamma.positive[i,j]=t(beta.mat[,j])%*%Sigma%*%beta.mat[,i]
    }
  }
  eigen.val<-eigen(Gamma.positive)$values
  if(eigen.val[L]<0){
    Gamma.positive=Gamma.positive+diag(2*eigen.val[L],L)
  }
  v <- Variable(L)
  objective <- Minimize(quad_form(v, Gamma.positive))
  constraints <- list(v >= 0, sum(v) == 1)
  prob.weight <- Problem(objective, constraints)
  result <- solve(prob.weight)
  opt.weight= result$getValue(v)
  maxmin=beta.mat%*%opt.weight
  
  return(list(opt.weight=opt.weight, beta.mm=maxmin))
}


FAIRM.lm<-function(X, y, n.vec, X.til, y.til){
#  X=X1;y=y1; n.vec=n1.vec
  p<-ncol(X)
  K=length(n.vec)
  alpha<-n.vec/sum(n.vec)
  Sig.list<-NULL
  Sig.tr<-diag(0,p)
  Mhat<-matrix(0,nrow=p,ncol=K)
  X.sc<-NULL
  for(k in 1:K){
    X.k<-scale(X[ind.set(n.vec,k),])
    y.k<-y[ind.set(n.vec,k)]
    Sig.list[[k]]<-t(X.k)%*%X.k/n.vec[k]
    Sig.tr<-Sig.tr+alpha[k]*Sig.list[[k]]
    Mhat[,k]<-t(X.k)%*%y.k/n.vec[k]
    X.sc<-rbind(X.sc,X.k)
  }
  sc.X<-apply(X,2,sd)
  Mhat.tr<-t(X.sc)%*%y/length(y)
  M.stat = rowSums(sapply(1:K, function(k)  alpha[k]*(Mhat[, k] - Mhat.tr)^2))
  rho.M = sqrt(mean(y^2))* (log(p)^2 + sqrt(K*log(p))) / sum(n.vec)
  #plot(M.stat);abline(h=rho.M);points(data.all$Sv, M.stat[data.all$Sv],col='red')
  
  #step 1
  Im.hat<-which(M.stat<=rho.M)
  Chat<-list()
  for(j in 1:p){Chat[[j]]=j}
  D.hat<-diag(0,p)
  for(k in 1:K){
    D.hat<-D.hat+alpha[k]*(Sig.list[[k]]-Sig.tr)^2
  }
  rho.D = sqrt(mean(diag(Sig.list[[1]])))* (log(p)^2 + sqrt(K*log(p))) / sum(n.vec)
  for(j in Im.hat){
    Ij.hat<-intersect(which(D.hat[j,]<=rho.D),Im.hat)
    if(max(D.hat[Ij.hat,Ij.hat])<=rho.D){
      Chat[[j]]<-Ij.hat
    }
  }
  ###step 2 & step 3
  Chat<-unique(Chat)
  hbeta.list<-list()
  mse<-rep(0,length(Chat))
  for(j in 1:length(Chat)){
    if(length(Chat[[j]])<=sqrt(length(y))){
      lm1<-lm(y~X[,Chat[[j]]]-1)
      hbeta.list[[j]]<-lm1$coef
    }else{
      cv.re<-cv.glmnet(x = X[,Chat[[j]]], y =y)
      hbeta.list[[j]]<-coef(cv.re,s='lambda.min')[-1]
    }
    if(length(Chat[[j]])==1){
      mse[j]<- mean((y.til-X.til[,Chat[[j]]]*hbeta.list[[j]])^2)
    }else{
      mse[j]<- mean((y.til-X.til[,Chat[[j]]]%*%hbeta.list[[j]])^2)
    }
  }
  j.hat<-which.min(mse)
  CC<-Chat[[j.hat]]
  beta.hat<-rep(0,p)
  beta.hat[CC]<-hbeta.list[[j.hat]]
  list(Shat=CC,beta.hat=beta.hat)
}
#######simulation begins#######
set.seed(123)
p = 400
K.tr = 10
K.te=3
s = 10
n.test=100 #number of test samples
Niter = 200 #number of iterations
setting=2
delta=0.4
h=5
#for(delta in c(0.2, 0.6, 1)){
 for(K.tr in seq(12,4,-2)){
  sse.inv <- matrix(0, ncol = 11, nrow = Niter)
  sse.erm <- matrix(0, ncol=9, nrow=Niter)
  sse.mm <- matrix(0, ncol=9, nrow=Niter)
  sse.oracle <- matrix(0, ncol=9, nrow=Niter)
  for (iter in 1:Niter) {
    n.vec <- rep(100, K.tr)
    n.test=n.vec[1]
    data.all <- Data.gen(p = p, s = s, h = h, K.tr = K.tr, K.te=K.te, delta = delta,
                         n.vec = n.vec, n.test=n.test, setting=setting)
    X.tr <- data.all$X.tr
    y.tr <- data.all$y.tr
    X.test<-data.all$X.test
    y.test<-data.all$y.test
    alpha<-n.vec/sum(n.vec)
    ##sample splitting
    #X1<-NULL;y1<-NULL
    #n1.vec<-rep(0,K.tr)
    #tn.vec<-rep(0,K.tr)
    # X.til<-NULL
    # y.til<-NULL
    # for(k in 1:K.tr){
    #   samp1<-sample(1:n.vec[k], n.vec[k]/2)
    #   X.k<-X.tr[ind.set(n.vec,k)[samp1],]
    #   y.k<-y.tr[ind.set(n.vec,k)[samp1]]
    #   X1<-rbind(X1, X.k)
    #   y1<-c(y1, y.k)
    #   n1.vec[k]<-length(samp1)
    #   tn.vec[k]<-n.vec[k]-n1.vec[k]
    #   X.til<-rbind(X.til, X.tr[ind.set(n.vec,k)[-samp1],])
    #   y.til<-c(y.til,y.tr[ind.set(n.vec,k)[-samp1]])
    # }
    #we realize the method without sample splitting
    FAIRM.re<-FAIRM.lm(X=X.tr,y=y.tr, n.vec=n.vec, X.til=X.tr,y.til=y.tr)
    beta.inv<-FAIRM.re$beta.hat
    Sv.hat<-setdiff(1:p,FAIRM.re$Shat)
    cat(sum(data.all$Sv%in%FAIRM.re$Shat)/length(data.all$Sv),sum((1:s)%in%FAIRM.re$Shat)/s,'\n')
    ###ERM
    cv.erm <- cv.glmnet(scale(X.tr), y = y.tr)
    beta.erm<-coef(cv.erm,s='lambda.1se')[-1]
    beta.erm<-beta.erm/apply(X.tr, 2, sd)

    ###Maximin
    beta.mat<-matrix(nrow=p,ncol=K.tr)
    for(k in 1:K.tr){
      X.sc.k<-scale(X.tr[ind.set(n.vec,k),])
      y.k<-y.tr[ind.set(n.vec,k)]
      cv.lasso<-cv.glmnet(x=X.sc.k, y=y.k)
      beta.hat<-coef(cv.lasso,s='lambda.min')[-1]
      beta.mat[,k]<-beta.hat/apply(X.tr[ind.set(n.vec,k),],2,sd)
    }
    mm.re<-Maximin(beta.mat,t(X.tr)%*%X.tr/sum(n.vec))
    beta.mm<-mm.re$beta.mm
    ###oracle
    beta.oracle<-rep(0,p)
    cv.ora<-cv.glmnet(X.tr[,-data.all$Sv],y.tr)
    beta.oracle[-data.all$Sv]<-coef(cv.ora,s='lambda.min')[-1]
    #plot(beta.inv-beta.oracle)

    ###training errors
    mse.tr<-matrix(0,ncol=K.tr,nrow=4)
    bias.tr<-matrix(0,ncol=K.tr,nrow=4)
    for(k in 1: K.tr){
      X.k<-X.tr[ind.set(n.vec,k),]
      y.k<-y.tr[ind.set(n.vec,k)]
      mse.tr[1,k]<-mean((y.k-X.k%*%beta.inv)^2)
      bias.tr[1,k]<-mean((y.k-X.k%*%beta.inv))
      mse.tr[2,k]<-mean((y.k-X.k%*%beta.erm)^2)
      bias.tr[2,k]<-mean((y.k-X.k%*%beta.erm))
      mse.tr[3,k]<-mean((y.k-X.k%*%beta.mm)^2)
      bias.tr[3,k]<-mean((y.k-X.k%*%beta.mm))
      mse.tr[4,k]<-mean((y.k-X.k%*%beta.oracle)^2)
      bias.tr[4,k]<-mean((y.k-X.k%*%beta.oracle))
    }    
    ###test errors
    tn.vec<-rep(n.test, K.te+K.tr)
    mse.te<-matrix(0,ncol=K.te+K.tr,nrow=4)
    bias.te<-matrix(0,ncol=K.te+K.tr,nrow=4)
    for(k in 1:(K.te+K.tr)){
      X.k<-X.test[ind.set(tn.vec,k),]
      y.k<-y.test[ind.set(tn.vec,k)]
      mse.te[1,k]<-mean((y.k-X.k%*%beta.inv)^2)/mean(y.k^2)
      bias.te[1,k]<-mean(y.k-X.k%*%beta.inv)
      mse.te[2,k]<-mean((y.k-X.k%*%beta.erm)^2)/mean(y.k^2)
      bias.te[2,k]<-mean(y.k-X.k%*%beta.erm)
      mse.te[3,k]<-mean((y.k-X.k%*%beta.mm)^2)/mean(y.k^2)
      bias.te[3,k]<-mean(y.k-X.k%*%beta.mm)
      mse.te[4,k]<-mean((y.k-X.k%*%beta.oracle)^2)/mean(y.k^2)
      bias.te[4,k]<-mean(y.k-X.k%*%beta.oracle)
    }
    sse.inv[iter, ] <- c(K.tr, mean(mse.tr[1,]), max(mse.tr[1,]), mean(mse.te[1,]), max(mse.te[1,]),
                        mean(mse.te[1,-(1:K.tr)]), max(mse.te[1,-(1:K.tr)]), mean(abs(bias.te[1,])), 
                        max(abs(bias.te[1,])), sum(data.all$Sv%in%FAIRM.re$Shat)/length(data.all$Sv),sum((1:s)%in%FAIRM.re$Shat)/s)
    sse.erm[iter, ] <-  c(K.tr,mean(mse.tr[2,]), max(mse.tr[2,]), mean(mse.te[2,]), max(mse.te[2,]),
                          mean(mse.te[2,-(1:K.tr)]), max(mse.te[2,-(1:K.tr)]),mean(abs(bias.te[2,])), max(abs(bias.te[2,])))
    sse.mm[iter, ] <-  c(K.tr,mean(mse.tr[3,]), max(mse.tr[3,]), mean(mse.te[3,]), max(mse.te[3,]),
                         mean(mse.te[3,-(1:K.tr)]), max(mse.te[3,-(1:K.tr)]), mean(abs(bias.te[3,])), max(abs(bias.te[3,])))
    sse.oracle[iter, ]<-c(K.tr,mean(mse.tr[4,]), max(mse.tr[4,]), mean(mse.te[4,]), max(mse.te[4,]),
                          mean(mse.te[4,-(1:K.tr)]), max(mse.te[4,-(1:K.tr)]), mean(abs(bias.te[4,])), max(abs(bias.te[4,])))

    #cat('iter=', iter, sse.inv[iter,], '\n')
    #cat('iter=', iter, sse.erm[iter,], '\n')
    #cat('iter=', iter, sse.oracle[iter,], '\n')
  }
  write.table(sse.inv,file=paste('Simu',setting,'-inv-delta',delta,'.txt',sep=''), append=T, row.names=F, col.names=F)
  write.table(sse.erm,file=paste('Simu',setting,'-erm-delta',delta,'.txt',sep=''),append=T, row.names=F, col.names=F)
  write.table(sse.mm, file= paste('Simu',setting,'-mm-delta',delta,'.txt',sep=''),append=T, row.names=F, col.names=F)
  write.table(sse.oracle,file=paste('Simu',setting,'-ora-delta',delta,'.txt',sep=''),append=T, row.names=F, col.names=F)
 }


