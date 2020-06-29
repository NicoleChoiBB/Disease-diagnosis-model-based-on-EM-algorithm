bootEM<-function (lambda,mu,sigma,samplesize,B = 100)
{ library(MASS)
  library(mixtools)
  k = length(lambda)
  n = samplesize
  p = length(mu)
  j = 0
  lambda.bs = NULL
  mu.bs = NULL
  sigma.bs = NULL
  while (j < B) {
    j = j + 1
    w = rmultinom(n, size = 1, prob = lambda)
    y.sim = t(sapply(1:n, function(i) rmvnorm(1, mu = mu[w[,i] == 1][[1]], sigma = sigma[w[, i] == 1][[1]])))
    em.bs = try(mvnormalmixEM(x = y.sim, k = k, arbmean = T,arbvar = T, lambda = lambda, mu = mu,sigma = sigma), silent = TRUE)
    if (class(em.bs) == "try-error" || em.bs$restarts != 
        0) {
      j = j - 1
    }
    else {
      lambda.bs = cbind(lambda.bs, em.bs$lambda)
      mu.bs1 = as.vector(sapply(em.bs$mu, c))
      mu.bs = cbind(mu.bs, mu.bs1)
      sigma.bs1 = lapply(1:k, function(i) as.vector(em.bs$sigma[[i]]))
      sigma.bs2 = as.vector(sapply(sigma.bs1, c))
      sigma.bs = cbind(sigma.bs, sigma.bs2)
    }
  }
  lambda.se = apply(lambda.bs, 1, sd)
  mu.se1 = apply(mu.bs, 1, sd)
  mu.se = lapply(1:k, function(i) mu.se1[((i - 1) * p + 1):(i * p)])
  sigma.se1 = apply(sigma.bs, 1, sd)
  sigma.se = lapply(1:k, function(i) matrix(sigma.se1[((i-1)*(p^2) + 1):(i * (p^2))], nrow = p, ncol = p))
  bs.list = list(lambda = lambda.bs, lambda.se = lambda.se,mu = mu.bs, mu.se = mu.se, sigma = sigma.bs, sigma.se = sigma.se,se=se)
}

#plot ROC
plot.roc<-function(mu,sigma,mu.sim,sigma.sim,cex=0.8){
  u0<-mu.sim[[1]]
  u1<-mu.sim[[2]]
  sigma0<-as.matrix(sigma.sim[[1]])
  sigma1<-as.matrix(sigma.sim[[2]])
  a<-solve(sigma0+sigma1)%*%(u1-u0)
  AUC_AB.sim<-pnorm(((u1-u0)%*%a)^(1/2))
  ci<-pretty(c(0,200),500)
  FPR<-rep(0,length(ci))
  ROC<-rep(0,length(ci))
  for (i in 1:length(ci)) {
    FPR[i]<-pnorm((-ci[i]+t(a)%*%u0)/(t(a)%*%sigma0%*%a)^(1/2))
    ROC[i]<-pnorm(qnorm(FPR[i])*(t(a)%*%sigma0%*%a)^(1/2)/(t(a)%*%sigma1%*%a)^(1/2)+t(a)%*%(u1-u0)/(t(a)%*%sigma1%*%a)^(1/2))
  }
  FPR.sim<-FPR
  ROC.sim<-ROC
  u0<-mu[[1]]
  u1<-mu[[2]]
  sigma0<-as.matrix(sigma[[1]])
  sigma1<-as.matrix(sigma[[2]])
  a<-solve(sigma0+sigma1)%*%(u1-u0)
  AUC_AB<-pnorm(((u1-u0)%*%a)^(1/2))
  ci<-pretty(c(0,200),500)
  FPR<-rep(0,length(ci))
  ROC<-rep(0,length(ci))
  for (i in 1:length(ci)) {
    FPR[i]<-pnorm((-ci[i]+t(a)%*%u0)/(t(a)%*%sigma0%*%a)^(1/2))
    ROC[i]<-pnorm(qnorm(FPR[i])*(t(a)%*%sigma0%*%a)^(1/2)/(t(a)%*%sigma1%*%a)^(1/2)+t(a)%*%(u1-u0)/(t(a)%*%sigma1%*%a)^(1/2))
  }
  plot(FPR,ROC,'l',col='green',lwd=2,main = 'the ROC curve of simlation and true parameter')
  points(FPR.sim,ROC.sim,'l',col='red',lwd=2)
  legend('bottomright',legend = c('Simulation','True ROC'),col = c('red','green'),lwd = 2,cex=cex)
  AUC<-list(AUC_AB.sim=AUC_AB.sim,AUC_AB=AUC_AB)
}
#relative bias¡¢AUC¡¢roc
showresult<-function(result,lambda,sigma,mu,cex=0.8){
  lamda.bias<-(mean(result$lambda[1,])-lambda[1])/lambda[1]
  mu.bias<-(apply(result$mu, 1, FUN = 'mean')-as.vector(sapply(mu, c)))/as.vector(sapply(mu, c))
  sigma.bias<-(apply(result$sigma, 1, FUN = 'mean')-as.vector(sapply(sigma, c)))/as.vector(sapply(sigma, c))
  sigma.mean<-lapply(1:2, function(i) matrix(apply(result$sigma, 1, FUN = 'mean')[((i-1)*(2^2) + 1):(i * (2^2))], nrow = 2, ncol = 2))
  mu.mean<-lapply(1:2, function(i) apply(result$mu,1,FUN = 'mean')[((i - 1) * 2 + 1):(i * 2)])
  AUC<-plot.roc(mu,sigma,mu.mean,sigma.mean,cex=cex)
  list(lamda.bias=lamda.bias,mu.bias=mu.bias,sigma.bias=sigma.bias,AUC=AUC)
}
boostrap_EMse<-function(data,lambda,mu,sigma,m){
  
  j<-0
  k<-2
  lambda.bs = NULL
  mu.bs = NULL
  sigma.bs = NULL
  library(mixtools)
  library(MASS)
  EM<-mvnormalmixEM(as.matrix(data),arbvar = T,lambda = lambda, mu = mu,sigma = sigma)
  lambda<-EM$lambda
  mu<-EM$mu
  sigma<-EM$sigma
  while (j<m) {
    j<-j+1
    r<-sample(1:length(data[,1]),length(data[,1]),T)
    y.sim<-data[r,]
    em.bs = try(mvnormalmixEM(x = y.sim, k = 2, arbmean = T,arbvar = T, lambda = lambda, mu = mu,sigma = sigma), silent = TRUE)
    if (class(em.bs) == "try-error" || em.bs$restarts != 0) {
      j = j - 1
    }
    else {
      lambda.bs = cbind(lambda.bs, em.bs$lambda)
      mu.bs1 = as.vector(sapply(em.bs$mu, c))
      mu.bs = cbind(mu.bs, mu.bs1)
      sigma.bs1 = lapply(1:k, function(i) as.vector(em.bs$sigma[[i]]))
      sigma.bs2 = as.vector(sapply(sigma.bs1, c))
      sigma.bs = cbind(sigma.bs, sigma.bs2)
    }
  }
  lambda.se = apply(lambda.bs, 1, sd)
  mu.se1 = apply(mu.bs, 1, sd)
  mu.se = lapply(1:k, function(i) mu.se1[((i - 1) * p + 1):(i * p)])
  sigma.se1 = apply(sigma.bs, 1, sd)
  sigma.se = lapply(1:k, function(i) matrix(sigma.se1[((i-1)*(p^2) + 1):(i * (p^2))], nrow = p, ncol = p))
  bs.list = list(lambda = lambda.bs, lambda.se = lambda.se,mu = mu.bs, mu.se = mu.se, sigma = sigma.bs, sigma.se = sigma.se)
}
#simulation 1
mu<-list(c(1,1),c(3,3))
lambda<-c(0.5,0.5)
sigma<-list(matrix(c(1,0,0,1),2,2),matrix(c(2,0,0,2),2,2))
result1.1<-bootEM(lambda,mu,sigma,samplesize=100,B = 100)
result1.2<-bootEM(lambda,mu,sigma,samplesize=200,B = 100)
showresult(result1.1,lambda,sigma,mu)
showresult(result1.2,lambda,sigma,mu)
#simulation 2
mu<-list(c(1,1),c(3,3))
lambda<-c(0.3,0.7)
sigma<-list(matrix(c(1,0,0,1),2,2),matrix(c(2,0,0,2),2,2))
result2.1<-bootEM(lambda1,mu,sigma,samplesize=100,B = 100)
showresult(result2.1,lambda,sigma,mu)

#simulation 3
mu<-list(c(1,1),c(2,2))
lambda<-c(0.5,0.5)
sigma<-list(matrix(c(1,0,0,1),2,2),matrix(c(2,0,0,2),2,2))
result3.1<-bootEM(lambda,mu1,sigma,samplesize=100,B = 100)
showresult(result3.1,lambda,sigma,mu)

#simulrtion 4
mu<-list(c(1,1),c(3,3))
lambda<-c(0.5,0.5)
sigma1<-list(matrix(c(1,0.5,0.5,1),2,2),matrix(c(2,1,1,2),2,2))
sigma2<-list(matrix(c(1,0.7,0.7,1),2,2),matrix(c(2,1.4,1.4,2),2,2))
result4.1<-bootEM(lambda,mu,sigma1,samplesize=100,B = 100)
result4.2<-bootEM(lambda,mu,sigma2,samplesize=100,B = 100)
showresult(result4.1,lambda,sigma1,mu)
showresult(result4.2,lambda,sigma2,mu)
