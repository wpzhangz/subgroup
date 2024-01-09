#n=200,p=300, k = 2
options(digits = 4)
ntr = 200
p = 300
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

#求已知组情形的MSE
group_true = ifelse(class1 == TRUE,1,2)
sub1_true = which(class1 == TRUE)
sub2_true = which(class1 == FALSE)
x_sub1 = X[sub1_true,]
y_sub1 = y[sub1_true]
x_sub2 = X[sub2_true,]
y_sub2 = y[sub2_true]
lambda.rw1 = cv.ncvreg(x_sub1,y_sub1,penalty = "SCAD",nfolds = 3)$lambda.min
beta1_scad = ncvfit(x_sub1,y_sub1,penalty = "SCAD",lambda = lambda.rw1)$beta
lambda.rw2 = cv.ncvreg(x_sub2,y_sub2,penalty = "SCAD",nfolds = 3)$lambda.min
beta2_scad = ncvfit(x_sub2,y_sub2,penalty = "SCAD",lambda = lambda.rw2)$beta
scad_mse = sqrt(length(x_sub1) * t(beta1_scad - beta1_true) %*% (beta1_scad - beta1_true) + 
                  length(x_sub2) * t(beta2_scad - beta2_true) %*% (beta2_scad - beta2_true))/ntr

#fmr方法
#sis 筛一部分变量
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

wnz = sum((knnmat + t(knnmat))!=0)/2 
Dmat = matrix(0, ncol=ntr, nrow=wnz)
## nrow of Dmat is determined by # of nonzero in W
rowIndex=1
for(i in 1:(ntr-1)){
  for(j in (i+1):ntr){
    if(knnmat[i,j] ==0){
    }else {
      if(knnmat[j,i] == 0){
        Dmat[rowIndex, i] = 1
        Dmat[rowIndex, j] = -1 
        rowIndex = rowIndex + 1
      }else{
        Dmat[rowIndex, i] = 2
        Dmat[rowIndex, j] = -2 
        rowIndex = rowIndex + 1
      }
    }
  }
}
for (i in 2:ntr) {
  for (j in 1:(i-1)) {
    if(knnmat[i,j] == 0){
    }else{
      if(knnmat[j,i] == 1){
      }else{
        Dmat[rowIndex, i] = -1
        Dmat[rowIndex, j] = 1 
        rowIndex = rowIndex + 1
      }
    }
  }
}
fit = APG_fit_L1(x = X,y = y,beta0 = beta0.0,lam1 = lam1,lam2 = lam2,rho = 1,wk = 10,
                 Dmat = Dmat,W = knnmat,L1multiplier = 40,t = 0.3,maxiter = 1000
                 ,nADMM = 1500)
fit_trunc = ifelse(abs(fit) < 0.2,0,fit)
mse_fit = norm((fit - beta_true),type = 'f')/ntr
group = kmeans(fit_trunc,centers = nComp)
loss = loss_fit(X,y,group)
B = BIC_est(X,y,loss,fit_trunc,nComp)
bic = B$BIC
tab = table(group$cluster,group_true)
RI = classAgreement(tab)$rand
number.v = B$q
coef = group$centers
result =  list(mse_ora = scad_mse,mse_fit = mse_fit,RI = RI,RTO = number.v,bic = bic,beta_est = coef)