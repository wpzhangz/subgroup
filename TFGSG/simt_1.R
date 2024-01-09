library(FGSG)
library(MASS)
library(glmnet)
library(SIS)
library(fmrs)
library(e1071)


seed = seq(10,1000,10)
repeat_times = length(seed)
RI_1 = numeric(repeat_times)
MSE_1 = numeric(repeat_times)
RTO_1 = numeric(repeat_times)
for (t in 1:repeat_times) {
  set.seed(seed[t])
  options(digits = 4)
  ntr = 200
  p = 200
  b = 0
  beta1_true = c(-2,-2,0,0,rep(0,p-4))
  beta2_true = c(0,0,2,2,rep(0,p-4))
  epsilon = rnorm(ntr,0,.1)
  sigma = matrix(0.3,p,p)
  diag(sigma) = rep(1,p)
  X = matrix(data = mvrnorm(ntr,rep(0,p),sigma) ,nrow = ntr)
  #X = matrix(data = rnorm(ntr*p,0,1) ,nrow = ntr)
  x1 = X[,1]
  x2 = X[,2]
  x3 = X[,3]
  x4 = X[,4]
  orderX = order(x1+x2+x3+x4)
  X = X[orderX,]
  mu = I(X[,1]+X[,2]+X[,3]+X[,4]<=b)*(X%*%beta1_true) + 
    I(X[,1]+X[,2]+X[,3]+X[,4]>b)*(X%*%beta2_true)
  y = mu + epsilon
  class1 = X[,1]+X[,2]+X[,3]+X[,4] <= b
  beta1 = ifelse(class1, -2, 0) 
  beta2 = ifelse(class1, -2, 0)
  beta3 = ifelse(class1, 0, 2)
  beta4 = ifelse(class1, 0, 2) ## slope
  num = rep(0,ntr)
  beta_true <- c(beta1, beta2,beta3,beta4,rep(num,p-4))
  
  beta_true.m = matrix(beta_true,ntr,p)
  
  group_true = ifelse(class1 == TRUE,1,2)
  
  nComp = 2
  nObs = ntr
  p.sis = SIS(X,y,tune='bic',varISIS='aggr',seed=10)$ix
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
  
  
  #our methods:
  #initial values:
  beta0 = initial_ridge(X.sis,y,lam = 0.001)
  #beta0 = get_local_init_beta2(X.sis,y,kn = 30)
  beta0.0 = matrix(0,ntr,p)
  beta0.0[,p.sis] = beta0
  #利用初值β0得到一组knn:knnmat_beta0
  dist.beta0 = as.matrix(dist(beta0))
  knnmat_beta0 <- matrix(0, ntr, ntr)# 1/0 matrix for knn
  wk = 10
  for(i in 1:ntr){
    knnmati = dist.beta0[i, ]
    knnid = order(knnmati, decreasing=FALSE)[1:(wk+1)]
    knnmat_beta0[i, knnid] = 1
    #knnmat[knnid, i] = 1
  }
  diag(knnmat_beta0) = rep(0, ntr)
  #利用fmr得到的分组得到knn:knnmat_fmr
  w = matrix((rbinom(ntr^2,1,0.05) * runif(ntr^2,0,1)),ntr,ntr)
  for (i in 1:nComp) {
    sub = which(f2[,i] >= 0.5)
    for (j in 1:length(sub)) {
      for (l in 1:length(sub)) {
        w[sub[j],sub[l]] = rbinom(1,1,0.2) * runif(1,1,2)
      }
    }
  }
  diag(w) = 0
  knnmat_fmr = ifelse(w !=0,1,0)
  #利用数据X,Y得到knnmat_xy
  D = createDmat(X.sis,y,wk = 10)
  knnmat_xy = D$knn
  
  #综合得到knnmat_case 得到最终的knnmat,利用此knnmat生成Dmat
  #简单的加起来
  knnmat = knnmat_beta0 + knnmat_fmr + knnmat_xy
  knnmat = ifelse(knnmat >= 2,1,0)
  
  
  graph = c()
  for (i in 1:ntr) {
    for (j in i:ntr) {
      if(knnmat[i,j] != 0){
        graph = append(graph,c(i,j))
      }
    }
  }
  
  X.sis.tilde = matrix(0,ntr,ntr*nCov.sis)
  for (i in 1:ntr) {
    a = (i-1)*nCov.sis + 1
    b = i*nCov.sis
    X.sis.tilde[i,a:b] = X.sis[i,]
  }
  
  
  
  fgsg_1 = ncTFGS(X.sis.tilde,y,graph,0.006,0.03)
  
  fgsg_1_m = matrix(unlist(fgsg_1),nrow = ntr)
  beta_true.m_s = beta_true.m[,p.sis]
  
  mse_fgsg = norm((beta_true.m_s - fgsg_1_m),type = 'f')/ntr
  
  group_fgsg = kmeans(fgsg_1_m,centers = nComp)
  tab_fgsg = table(group_fgsg$cluster,group_true)
  RI_fgsg = classAgreement(tab_fgsg)$rand
 
  
  nonzero = numeric()
  for (i in 1:nComp) {
    center = group_fgsg$centers[i,]
    center_trunc = ifelse(abs(center) < 0.1,0,center)
    nonzero = append(nonzero,which(center_trunc != 0))
  }
  RTO_fgsg = length(unique(nonzero))
  
  RI_1[t] = RI_fgsg
  MSE_1[t] = mse_fgsg
  RTO_1[t] = RTO_fgsg
}