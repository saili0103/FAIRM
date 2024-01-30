library(mvtnorm)
library(glmnet)
library(CVXR) #compute Maximin

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

#Maximin method
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




logit<-function(x){1/(1+exp(-x))}
####FAIRM for GLM used on the color MNIST experiments
FAIRM.glm<-function(X, y, n.vec){
  #  X=X.tr;y=y.tr; n.vec=n.vec
  X<-as.matrix(X)
  X.til<-X
  y.til<-y
  X.til<-as.matrix(X.til)
  p<-ncol(X)
  K=length(n.vec)
  alpha<-n.vec/sum(n.vec)
  Sig.list<-NULL
  Sig.tr<-diag(0,p)
  Mhat<-matrix(0,nrow=p,ncol=K)
  X.sc<-NULL
  for(k in 1:K){
    X.k<-scale(X[ind.set(n.vec,k),])
    exl.set<-which(apply(X[ind.set(n.vec,k),],2,sd)==0)
    X.k[,exl.set] <- 0
    y.k<-y[ind.set(n.vec,k)]
    Sig.list[[k]]<-t(X.k)%*%X.k/n.vec[k]
    Sig.tr<-Sig.tr+alpha[k]*Sig.list[[k]]
    Mhat[,k]<-t(X.k)%*%y.k/n.vec[k]
    X.sc<-rbind(X.sc,X.k)
  }
  X.til<-as.matrix(X.til)
  sc.X<-apply(X,2,sd)
  Mhat.tr<-t(X.sc)%*%y/length(y)
  M.stat = rowSums(sapply(1:K, function(k)  alpha[k]*(Mhat[, k] - Mhat.tr)^2))
  rho.M = var(y^2)* (log(p)^2 + sqrt(K*log(p))) / sum(n.vec)
  #plot(M.stat);abline(h=rho.M);points(data.tr1$bg.ind,M.stat[data.tr1$bg.ind],col='red')
  
  #step 1
  Im.hat<-which(M.stat<=rho.M)
  Chat<-list()
  D.hat<-diag(0,p)
  for(k in 1:K){
    D.hat<-D.hat+alpha[k]*(Sig.list[[k]]-Sig.tr)^2
  }
  rho.D = sqrt((log(p)^2 + sqrt(K*log(p))) / sum(n.vec))
  for(k in 1:length(Im.hat)){
    j<-Im.hat[k]
    Ij.hat<-intersect(which(D.hat[j,]<=rho.D*sqrt(max(D.hat[j,]))),Im.hat)
    #cat(j,max(D.hat[Ij.hat,Ij.hat]),'\n')
    if(max(D.hat[Ij.hat,Ij.hat])<=rho.D){
      Chat[[k]]<-Ij.hat
    }else{
      Chat[[k]]<-j
    }
  }
  #max(unlist(lapply(Chat,length)))
  ###step 2 & step 3
  Chat<-unique(Chat)
  hbeta.list<-list()
  mse<-rep(0,length(Chat))
  for(j in 1:length(Chat)){
      cv.re<-cv.glmnet(x = as.matrix(X[,Chat[[j]]]), y =y, family='binomial')
      hbeta.list[[j]]<-as.numeric(coef(cv.re,s='lambda.min'))
    
    if(length(Chat[[j]])==1){
      mse[j]<- mean(abs(y.til-logit(X.til[,Chat[[j]]]*hbeta.list[[j]][-1]+hbeta.list[[j]][1])))
    }else{
      mse[j]<- mean(abs(y.til-logit(X.til[,Chat[[j]]]%*%hbeta.list[[j]][-1]+hbeta.list[[j]][1])))
    }
  }
  j.hat<-which.min(mse)
  CC<-Chat[[j.hat]]
  beta.hat<-rep(0,p+1)
  beta.hat[c(1,CC+1)]<-hbeta.list[[j.hat]]
  list(Shat=CC,beta.hat=beta.hat)
}


####creating color MNIST data used in Section 5.2
cmnist.gen1<-function(pe,x,y){
  n<-length(y)
  y.til<-y
  rand=rbinom(n,1,pe)
  z=rand*y+(1-rand)*(1-y) #z=y with probability pe
  x.final<-x
  #create frame
  bg.ind=1:784
  bg.coord<-cbind(ceiling(bg.ind/28),bg.ind%%28)
  sub1<-which(bg.coord[,1]<=4 | bg.coord[,1]>=25| bg.coord[,2]<=5|bg.coord[,2]>=23)
  #plot(bg.coord[sub1,])
  bg.ind<-(bg.coord[sub1,1]-1)*28+bg.coord[sub1,2]+28*(bg.coord[sub1,2]==0)
  
  for(i in 1:n){
    #temp1<-which(x[i,]<=20) #background
    x.final[i,bg.ind] <-sapply(x.final[i,bg.ind]+z[i]*100,function(xx) min(xx, 255))
  }
  list(y=y,x=x.final, bg.ind=bg.ind)
}


library(dslabs)
#download mnist data from http://yann.lecun.com/exdb/mnist/
#need to modify the file names according to the warnings
mnist <- read_mnist("~/Downloads")

image(1:28, 1:28, matrix(mnist$train$images[4,], nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
image(1:28, 1:28, matrix(mnist$train$images[4,], nrow=28)[ , 28:1],
      col=hcl.colors(20, "YlOrRd", rev = TRUE), xlab = "", ylab="")
image(1:28, 1:28, matrix(mnist$train$images[4,], nrow=28)[ , 28:1],
      col=hcl.colors(20, "terrain"), xlab = "", ylab="")
image(1:28, 1:28, matrix(mnist$train$images[16,], nrow=28)[ , 28:1],
      col=hcl.colors(20, "YlOrRd", rev = TRUE), xlab = "", ylab="")
image(1:28, 1:28, matrix(mnist$train$images[16,], nrow=28)[ , 28:1],
      col=hcl.colors(20, "terrain"), xlab = "", ylab="")

## the label for this image is: 
mnist$train$labels[4]
########Exp 1: 1 or 7 ########
set.seed(123)
n0=length(mnist$train$labels)
sub0<-which(mnist$train$labels==1 | mnist$train$labels==7) #classify 1 or 7
sub0<-sample(sub0, 1000)###only use first 1k samples
n<-length(sub0)#traning data sample size

#environment 1
y<-as.numeric(mnist$train$labels[sub0]>=5)
data.tr1<-cmnist.gen1(pe=0.9, x=mnist$train$images[sub0,], y=y)
X.tr1<-data.tr1$x
y.tr1<-data.tr1$y
image(1:28, 1:28, matrix(X.tr1[1,], nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
image(1:28, 1:28, matrix(X.tr1[3,], nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")


#environment 2
y<-as.numeric(mnist$train$labels[sub0]>=5)
data.tr2<-cmnist.gen1(pe=0.6, x=mnist$train$images[sub0,], y=y)
X.tr2<-data.tr2$x
y.tr2<-data.tr2$y

cmnist.tr.x<-rbind(X.tr1, X.tr2)
cmnist.tr.y<-c(y.tr1, y.tr2)
n.vec<-rep(length(y.tr1),2)
cor(cmnist.tr.x[,1],cmnist.tr.y)
#test environment 1
test.sub<-which(mnist$test$labels==1 | mnist$test$labels==7)
tn<-length(test.sub)
y<-as.numeric(mnist$test$labels[test.sub]>=5)
data.te1<-cmnist.gen1(pe=0.8, x=mnist$test$images[test.sub,], y=y)
x.te1<-data.te1$x
y.te1<-data.te1$y
cor(x.te1[,1],y.te1)
cor(X.tr1[,1],y.tr1)
#test environment 2
data.te2<-cmnist.gen1(pe=0.5, x=mnist$test$images[test.sub,], y=y)
x.te2<-data.te2$x
y.te2<-data.te2$y

#test 3
data.te3<-cmnist.gen1(pe=0.2, x=mnist$test$images[test.sub,], y=y)
x.te3<-data.te3$x
y.te3<-data.te3$y


cmnist.tr<-cbind(cmnist.tr.x,cmnist.tr.y)
p<-ncol(cmnist.tr) #785
anyNA(cmnist.tr)
write.table(cmnist.tr,file='cmnist-train-dat.txt')
cmnist.te.x<-rbind(x.te1, x.te2, x.te3)
cmnist.te.y=c(y.te1, y.te2, y.te3)
cmnist.te<-cbind(cmnist.te.x,cmnist.te.y)
p<-ncol(cmnist.te)
anyNA(cmnist.te)
write.table(cmnist.te,file='cmnist-test-dat.txt')

#######start training: 1vs 7#######
df1<-read.table('cmnist-train-dat.txt')
p<-ncol(df1)-1
X.tr<-df1[,1:p]
y.tr<-df1[,p+1]
K.tr=2
n.vec<-rep(nrow(df1)/K.tr,K.tr)
n1.vec<-rep(0,K.tr)
tn.vec<-rep(0,K.tr)


X.tr<-as.matrix(X.tr)
bg.ind=data.tr1$bg.ind
#X1<-X.tr[,-data.tr1$bg.ind[-1]]
###FAIRM###
inv.re<-FAIRM.glm(X=X.tr,y=y.tr, n.vec=n.vec)
beta.inv<-inv.re$beta.hat
length(inv.re$Shat)
x0<-rep(100,p)
x0[-bg.ind]<-0
x0[beta.inv[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")

#ERM
cv.erm <- cv.glmnet(x = X.tr, y = y.tr, family='binomial')
beta.erm <- as.numeric(coef(cv.erm, s='lambda.min') )
#beta.erm<-rep(0,p+1)
#beta.erm[-(bg.ind+1)]<-beta.erm0[1:(length(beta.erm0)-1)];beta.erm[bg.ind+1]<-beta.erm0[length(beta.erm0)]
x0<-rep(100,p)
x0[-bg.ind]=0
x0[beta.erm[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
#oracle
X.tr0<-X.tr[,-bg.ind]
cv.ora <- cv.glmnet(x = X.tr0, y = y.tr, family='binomial')
beta.ora0 <- as.numeric(coef(cv.ora, s='lambda.min') )
beta.ora<-rep(0,p+1)
beta.ora[1]<-beta.ora0[1]; beta.ora[-c(1,data.tr1$bg.ind+1)]<-beta.ora0[-1]

x0<-rep(100,p)
x0[-bg.ind]=0
x0[beta.ora[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")

#Maximin
beta.mat<-matrix(0,nrow=p+1,ncol=K.tr)
n.vec<-rep(nrow(df1)/K.tr,K.tr)
for(k in 1:K.tr){
  X.k<-X.tr[ind.set(n.vec,k),]
  y.k<-y.tr[ind.set(n.vec,k)]
  cv.lasso<-cv.glmnet(x=X.k, y=y.k, family='binomial')
  beta.mat[,k]<-as.numeric(coef(cv.lasso,s='lambda.min'))
}
X.tr<-as.matrix(X.tr)
mm.re<-Maximin(beta.mat[-1,],t(X.tr)%*%X.tr/sum(n.vec))
mm.re$opt.weight
beta.mm<-c(sum(beta.mat[1,]*mm.re$opt.weight),mm.re$beta.mm)

x0<-rep(100,p)
x0[-bg.ind]=0
x0[beta.mm[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
####training err & test err-logit####
library(CalibratR)
mse.tr<-matrix(0,ncol=K.tr,nrow=4)
for(k in 1: K.tr){
  X.k<-as.matrix(X.tr[ind.set(n.vec,k),])
  y.k<-y.tr[ind.set(n.vec,k)]
  mse.tr[1,k]<-mean(abs(y.k-logit(X.k%*%beta.inv[-1]+beta.inv[1])))
  mse.tr[2,k]<-mean(abs(y.k-logit(X.k%*%beta.erm[-1]+beta.erm[1])))
  mse.tr[3,k]<-mean(abs(y.k-logit(X.k%*%beta.mm[-1]+beta.mm[1])))
  mse.tr[4,k]<-mean(abs(y.k-logit(X.k%*%beta.ora[-1]+beta.ora[1])))
}    
rowMeans(mse.tr)

###test errors
K.te=3
df2<-read.table('cmnist-test-dat.txt')
X.test<-as.matrix(df2[,1:p])
y.test<-df2[,p+1]
tn.vec<-rep(nrow(df2)/K.te, K.te)
mse.te<-matrix(0,ncol=K.te,nrow=4)
ece<-matrix(0,ncol=K.te,nrow=4)
for(k in 1:K.te){
  X.k<-as.matrix(X.test[ind.set(tn.vec,k),])
  y.k<-y.test[ind.set(tn.vec,k)]
  mse.te[1,k]<-mean(abs(y.k-(logit(X.k%*%beta.inv[-1]+beta.inv[1])>=0.5)))
  mse.te[2,k]<-mean(abs(y.k-(logit(X.k%*%beta.erm[-1]+beta.erm[1])>=0.5)))
  mse.te[3,k]<-mean(abs(y.k-(logit(X.k%*%beta.mm[-1]+beta.mm[1])>=0.5)))
  mse.te[4,k]<-mean(abs(y.k-(logit(X.k%*%beta.ora[-1]+beta.ora[1])>=0.5)))
  ece[1,k]<-getECE(y.k,logit(X.k%*%beta.inv[-1]+beta.inv[1]))
  ece[2,k]<-getECE(y.k,logit(X.k%*%beta.erm[-1]+beta.erm[1]))
  ece[3,k]<-getECE(y.k,logit(X.k%*%beta.mm[-1]+beta.mm[1]))
  ece[4,k]<-getECE(y.k,logit(X.k%*%beta.ora[-1]+beta.ora[1]))
}

mse.te
ece
library(xtable)
xtable(cbind(mse.te[,1],mse.te[,2],mse.te[,3],ece[,1],ece[,2],ece[,3]),digits=c(0,4,4,4,3,3,3))

########Exp2: 0 or 6####
set.seed(123)
n0=length(mnist$train$labels)
sub0<-which(mnist$train$labels==0 | mnist$train$labels==6) #classify 0 or 6
sub0<-sample(sub0,1000)###only use first 1k samples
n<-length(sub0)#traning data sample size
#environment 1
y<-as.numeric(mnist$train$labels[sub0]>=5)
data.tr1<-cmnist.gen1(pe=0.9, x=mnist$train$images[sub0,], y=y)
X.tr1<-data.tr1$x
y.tr1<-data.tr1$y
image(1:28, 1:28, matrix(X.tr1[3,], nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
image(1:28, 1:28, matrix(X.tr1[7,], nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
cor(X.tr1[,1],y.tr1)
#environment 2
y<-as.numeric(mnist$train$labels[sub0]>=5)
data.tr2<-cmnist.gen1(pe=0.4, x=mnist$train$images[sub0,], y=y)
X.tr2<-data.tr2$x
y.tr2<-data.tr2$y
cor(X.tr2[,1],y.tr2)
cmnist.tr.x<-rbind(X.tr1, X.tr2)
cmnist.tr.y=c(y.tr1, y.tr2)
n.vec<-rep(length(y.tr1),2)


#test environment 1
test.sub<-which(mnist$test$labels==0 | mnist$test$labels==6)
tn<-length(test.sub)
y<-as.numeric(mnist$test$labels[test.sub]>=5)
data.te1<-cmnist.gen1(pe=0.8, x=mnist$test$images[test.sub,], y=y)
x.te1<-data.te1$x
y.te1<-data.te1$y
cor(x.te1[,1],y.te1)
cor(X.tr1[,1],y.tr1)
#test environment 2
data.te2<-cmnist.gen1(pe=0.5, x=mnist$test$images[test.sub,], y=y)
x.te2<-data.te2$x
y.te2<-data.te2$y

#test 3
data.te3<-cmnist.gen1(pe=0.2, x=mnist$test$images[test.sub,], y=y)
x.te3<-data.te3$x
y.te3<-data.te3$y


cmnist.tr<-cbind(cmnist.tr.x,cmnist.tr.y)
p<-ncol(cmnist.tr) #785
anyNA(cmnist.tr)
write.table(cmnist.tr,file='cmnist-train-dat06.txt')
cmnist.te.x<-rbind(x.te1, x.te2, x.te3)
cmnist.te.y=c(y.te1, y.te2, y.te3)
cmnist.te<-cbind(cmnist.te.x,cmnist.te.y)
p<-ncol(cmnist.te)
anyNA(cmnist.te)
write.table(cmnist.te,file='cmnist-test-dat06.txt')

#######start training#######
df1<-read.table('cmnist-train-dat06.txt')
p<-ncol(df1)-1
X.tr<-df1[,1:p]
y.tr<-df1[,p+1]

K.tr=2
n.vec<-rep(nrow(df1)/K.tr,K.tr)
#n1.vec<-rep(0,K.tr)
#tn.vec<-rep(0,K.tr)
# X.til<-NULL
# y.til<-NULL
# for(k in 1:K.tr){
#   samp1<-sample(1:n.vec[k], n.vec[k])
#   X.k<-X.tr[ind.set(n.vec,k)[samp1],]
#   y.k<-y.tr[ind.set(n.vec,k)[samp1]]
#   X1<-rbind(X1, X.k)
#   y1<-c(y1, y.k)
#   n1.vec[k]<-length(samp1)
#   tn.vec[k]<-n.vec[k]-n1.vec[k]
#   #X.til<-rbind(X.til, X.tr[ind.set(n.vec,k)[-samp1],])
#   #y.til<-c(y.til,y.tr[ind.set(n.vec,k)[-samp1]])
# }

##FAIRM###
X.tr<-as.matrix(X.tr)
inv.re<-FAIRM.glm(X=X.tr,y=y.tr, n.vec=n.vec)
beta.inv<-inv.re$beta.hat
length(inv.re$Shat)
bg.ind=data.tr1$bg.ind
x0<-rep(100,p)
x0[-bg.ind]=0
x0[beta.inv[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")

#ERM
cv.erm <- cv.glmnet(x = X.tr, y = y.tr, family='binomial')
beta.erm <- as.numeric(coef(cv.erm, s='lambda.min') )

x0<-rep(100,p)
x0[-bg.ind]=0
x0[beta.erm[-1]!=0]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")
#oracle
X.tr0<-X.tr[,-data.tr1$bg.ind]
cv.ora <- cv.glmnet(x = X.tr0, y = y.tr, family='binomial')
beta.ora0 <- as.numeric(coef(cv.ora, s='lambda.min') )
beta.ora<-rep(0,p+1)
beta.ora[1]<-beta.ora0[1]; beta.ora[-c(1,data.tr1$bg.ind+1)]<-beta.ora0[-1]

x0<-rep(100,p)
x0[-bg.ind]=0
x0[which(abs(beta.ora[-1])>cv.ora$lambda.min)]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")

#Maximin
beta.mat<-matrix(0,nrow=p+1,ncol=K.tr)
n.vec<-rep(nrow(df1)/K.tr,K.tr)
for(k in 1:K.tr){
  X.k<-X.tr[ind.set(n.vec,k),]
  y.k<-y.tr[ind.set(n.vec,k)]
  cv.lasso<-cv.glmnet(x=X.k, y=y.k, family='binomial')
  beta.mat[,k]<-as.numeric(coef(cv.lasso,s='lambda.min'))
}
X.tr<-as.matrix(X.tr)
mm.re<-Maximin(beta.mat[-1,],t(X.tr)%*%X.tr/sum(n.vec))
mm.re$opt.weight
beta.mm<-c(sum(beta.mat[1,]*mm.re$opt.weight),mm.re$beta.mm)
x0<-rep(100,p)
x0[-bg.ind]=0
x0[which(abs(beta.mm[-1])!=0)]<-256
image(1:28, 1:28, matrix(x0, nrow=28)[ , 28:1], 
      col = gray(seq(0, 1, 0.05)), xlab = "", ylab="")


####training err & test err-logit####
library(CalibratR)
mse.tr<-matrix(0,ncol=K.tr,nrow=4)
for(k in 1: K.tr){
  X.k<-as.matrix(X.tr[ind.set(n.vec,k),])
  y.k<-y.tr[ind.set(n.vec,k)]
  mse.tr[1,k]<-mean(abs(y.k-logit(X.k%*%beta.inv[-1]+beta.inv[1])))
  mse.tr[2,k]<-mean(abs(y.k-logit(X.k%*%beta.erm[-1]+beta.erm[1])))
  mse.tr[3,k]<-mean(abs(y.k-logit(X.k%*%beta.mm[-1]+beta.mm[1])))
  mse.tr[4,k]<-mean(abs(y.k-logit(X.k%*%beta.ora[-1]+beta.ora[1])))
}    
rowMeans(mse.tr)

###test errors
K.te=3
df2<-read.table('cmnist-test-dat06.txt')
X.test<-as.matrix(df2[,1:p])
y.test<-df2[,p+1]
tn.vec<-rep(nrow(df2)/K.te, K.te)
mse.te<-matrix(0,ncol=K.te,nrow=4)
ece<-matrix(0,ncol=K.te,nrow=4)
for(k in 1:K.te){
  X.k<-as.matrix(X.test[ind.set(tn.vec,k),])
  y.k<-y.test[ind.set(tn.vec,k)]
  mse.te[1,k]<-mean(abs(y.k-(logit(X.k%*%beta.inv[-1]+beta.inv[1])>=0.5)))
  mse.te[2,k]<-mean(abs(y.k-(logit(X.k%*%beta.erm[-1]+beta.erm[1])>=0.5)))
  mse.te[3,k]<-mean(abs(y.k-(logit(X.k%*%beta.mm[-1]+beta.mm[1])>=0.5)))
  mse.te[4,k]<-mean(abs(y.k-(logit(X.k%*%beta.ora[-1]+beta.ora[1])>=0.5)))
  ece[1,k]<-getECE(y.k,logit(X.k%*%beta.inv[-1]+beta.inv[1]))
  ece[2,k]<-getECE(y.k,logit(X.k%*%beta.erm[-1]+beta.erm[1]))
  ece[3,k]<-getECE(y.k,logit(X.k%*%beta.mm[-1]+beta.mm[1]))
  ece[4,k]<-getECE(y.k,logit(X.k%*%beta.ora[-1]+beta.ora[1]))
}

mse.te
ece
library(xtable)
xtable(cbind(mse.te[,1],mse.te[,2],mse.te[,3],ece[,1],ece[,2],ece[,3]),digits=c(0,4,4,4,3,3,3))




