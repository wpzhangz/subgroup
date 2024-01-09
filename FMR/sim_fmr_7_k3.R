library(SIS)
library(MASS)
library(e1071)
library(fmrs)
seed = seq(10,1000,10)
repeat_times = length(seed)
r_fmr_7_k3 = numeric(repeat_times)
#r_fmr_1_2 = numeric(repeat_times)
MSE_fmr_7_k3 = numeric(repeat_times)
RTO_fmr_7_k3 = numeric(repeat_times)
for (t in 1:repeat_times) {
  set.seed(seed[t])
  options(digits = 2)
  ntr = 300
  p = 200
  b1 = 0.37
  b2 = 1.4
  a = 2
  b = 0.5
  beta1_true = c(-a,0,rep(0,p-2))
  beta2_true = c(0,a,rep(0,p-2))
  beta3_true = c(-b,b,rep(0,p-2))
  epsilon = rnorm(ntr,0,.1)
  sigma = matrix(0.1,p,p)
  diag(sigma) = rep(1,p)
  X = matrix(data = mvrnorm(ntr,rep(0,p),sigma) ,nrow = ntr)
  x1 = X[,1]
  x2 = X[,2]
  orderX = order(x1^2+x2)
  X = X[orderX,]
  mu = I(X[,1]^2+X[,2]<=b1)*(X%*%beta1_true) + 
    I(X[,1]^2+X[,2]>b1)*I(X[,1]^2+X[,2]<=b2)*(X%*%beta2_true) +
    I(X[,1]^2+X[,2] >b2)*(X%*%beta3_true)
  y = mu + epsilon
  class1 = X[,1]^2+X[,2] <= b1
  class2 = X[,1]^2+X[,2] > b2
  beta1 = ifelse(class1, -a, ifelse(class2,-b,0)) 
  beta2 = ifelse(class1, 0,ifelse(class2,b,a))
  num = rep(0,ntr)
  beta_true <- c(beta1, beta2,rep(num,p-2))
  
  beta_true.m = matrix(beta_true,ntr,p)
  
  #求已知组情形的MSE
  group_true = ifelse(class1 == TRUE,1,ifelse(class2 == TRUE,3,2))
  
  nComp = 3
  nObs = ntr
  cor.xy = cor(X,y)
  #p.sis = SIS(X,y,tune='bic',varISIS='cons',seed=10)$ix
  #p.sis = SIS(X,y,tune='bic',varISIS='aggr',seed=10)$ix
  p.sis =  sort(order(abs(cor.xy),decreasing = TRUE)[1:20])
  X.sis = X[,p.sis]
  nCov.sis = ncol(X.sis)
  y.v = as.vector(y)
  delta = rep(1,ntr)
  res.mle.2 <- fmrs.mle(y = y.v, x = X.sis, delta = delta,
                        nComp = nComp, disFamily = "norm",
                        initCoeff = rnorm(nComp*nCov.sis+nComp),
                        initDispersion = rep(1, nComp),
                        initmixProp = rep(1/nComp, nComp))
  res.lam2 <- fmrs.tunsel(y = y.v, x = X.sis, delta = delta,
                          nComp = nComp, disFamily = "norm",
                          initCoeff = c(coefficients(res.mle.2)),
                          initDispersion = dispersion(res.mle.2),
                          initmixProp = mixProp(res.mle.2),
                          penFamily = "lasso")
  res.var2 <- fmrs.varsel(y = y.v, x = X.sis, delta = delta,
                          nComp = ncomp(res.mle.2), disFamily = "norm",
                          initCoeff=c(coefficients(res.mle.2)),
                          initDispersion = dispersion(res.mle.2),
                          initmixProp = mixProp(res.mle.2),
                          penFamily = "lasso",
                          lambPen = slot(res.lam2, "lambPen"))
  beta_fmr2 = round(coefficients(res.var2),5)
  f2 = res.var2@weights
  sub1.2 = which(f2[,1]>=0.5)
  sub2.2 = which(f2[,2]>=0.5)
  sub2.3 = which(f2[,3]>=0.5)
  beta_fmr2_dropi = beta_fmr2[-1,]
  
  s_fmr = c(which(beta_fmr2_dropi[,1]!=0),which(beta_fmr2_dropi[,2]!=0),
            which(beta_fmr2_dropi[,3]!=0))
  s_fmr = s_fmr[!duplicated(s_fmr)]
  rto = length(s_fmr)
  
  
  beta_fmr2_a = matrix(0,ntr,(length(p.sis)+1))
  beta_fmr2_a_t = t(beta_fmr2_a)
  beta_fmr2_a_t[,sub1.2] = beta_fmr2[,1]
  beta_fmr2_a_t[,sub2.2] = beta_fmr2[,2]
  beta_fmr2_a_t[,sub2.3] = beta_fmr2[,3]
  beta_fmr2_a = t(beta_fmr2_a_t)
  
  beta_true_sis = cbind(rep(0,ntr),beta_true.m[,p.sis])
  
  mse_fmr = norm((beta_fmr2_a - beta_true_sis),type = 'f')/ntr
  
  g.e = numeric(ntr)
  g.e[sub1.2] = 1
  g.e[sub2.2] = 2
  g.e[sub2.3] = 3
  
  t_fmr = table(g.e,group_true)
  r = classAgreement(t_fmr)$rand
  
  r_fmr_7_k3[t] = r
  MSE_fmr_7_k3[t] = mse_fmr
  RTO_fmr_7_k3[t] = rto
  #mse_spgr =  norm(beta_spgr-beta_true_spgr,type = 'f')/ntr
}