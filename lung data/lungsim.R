luad.y = read.csv("luad_y.csv",header = TRUE,sep = ",")
lusc.y = read.csv("lusc_y.csv",header = TRUE,sep = ",")
y1 = luad.y[,2]
y2 = lusc.y[,2]
n1 = length(y1)
n2 = length(y2)
y = c(y1,y2)
ntr = length(y)
group_true = c(rep(1,n1),rep(2,n2))

X = read.csv("X_anova.csv",header = TRUE,sep = ",")
X.df = as.data.frame(X)
label = X[,1]
X.data = matrix(unlist(X[,-1]),nrow = ntr,byrow = FALSE)
X1.data = X.data[1:n1,]
X2.data = X.data[(n1+1):ntr,]
p = ncol(X.data)

#根据方差筛
#策略一：根据方差筛完数据，利用sis得到的数据先进行模拟，得到fit，
#根据kmeans对结果fit进行聚类，同组使用同一初值：y_k=x_k beta_k
#此时的k取2，或较大一数，再进行组合并。
#k取较大时效果不佳（合并效果差）
#想法：k取一个范围，如k=(2,3,4,5,..),基于此k重新进行拟合并计算bic，最后取bic最小值

var.vari = apply(X.data,2,var)
#q1: 分位数
var.s = which(var.vari>quantile(var.vari,q1))
X.sis = X.data[,var.s]
p.s = ncol(X.sis) 
X.sisc = scale(X.sis,center = TRUE,scale = TRUE)
X.datac = scale(X.data,center = TRUE,scale = TRUE)
library(SIS)
p.sis = SIS(X.sisc,y,tune='bic',varISIS='aggr',seed=5)$ix
X.sp = X.sisc[,p.sis]


beta02 = initial_ridge(X.sp,y,lam = 0.0001)
beta0.0 = matrix(0,ntr,p.s)
beta0.0[,p.sis] = beta02
dist.beta0 = as.matrix(dist(beta02))
knnmat_beta0 <- matrix(0, ntr, ntr)# 1/0 matrix for knn
# wk : KNN 
wk = k1
for(i in 1:ntr){
  knnmati = dist.beta0[i, ]
  knnid = order(knnmati, decreasing=FALSE)[1:(wk+1)]
  knnmat_beta0[i, knnid] = 1
  #knnmat[knnid, i] = 1
}
diag(knnmat_beta0) = rep(0, ntr)
library(pheatmap)
pheatmap(knnmat_beta0,cluster_rows = FALSE,cluster_cols = FALSE)
knnmat = knnmat_beta0
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
fit = APG_fit_L1(x = X.sisc,y = y,beta0 = beta0.0,lam1 = lam1,lam2 = lam2,rho = 1,
                 wk = k1,Dmat = Dmat,W = knnmat,L1multiplier = L1,t = t,
                 maxiter = 1000,nADMM = 1000,B=B)

library(e1071)
group = kmeans(fit,centers = nComp)
tab = table(group$cluster,group_true)

RI = classAgreement(tab)$rand

class1 = which(group$cluster == 1)
class2 = which(group$cluster == 2)

group1.xall = X.datac[class1,]
group2.xall = X.datac[class2,]
group1.y = y[class1]
group2.y = y[class2]

library(ncvreg)
lambda.1a = cv.ncvreg(group1.xall,group1.y,penalty = "SCAD",nfolds = 3)$lambda.min
beta.1a = ncvfit(group1.xall,group1.y,penalty = "SCAD",lambda = lambda.1a)$beta

lambda.2a = cv.ncvreg(group2.xall,group2.y,penalty = "SCAD",nfolds = 3)$lambda.min
beta.2a = ncvfit(group2.xall,group2.y,penalty = "SCAD",lambda = lambda.2a)$beta

sis.all = c(which(beta.1a != 0),which(beta.2a != 0))

sis.all = sis.all[!duplicated(sis.all)]

p.a = length(sis.all)

X.sa = X.datac[,sis.all]

beta0.a = initial_ridge(X.sa,y,lam = 0.0001)

beta0.sa = matrix(0,ntr,p)
beta0.sa[,sis.all] = beta0.a
dist.beta0.sa = as.matrix(dist(beta0.sa))
knnmat_beta0.sa <- matrix(0, ntr, ntr)# 1/0 matrix for knn
wk = 10
for(i in 1:ntr){
  knnmati = dist.beta0.sa[i, ]
  knnid = order(knnmati, decreasing=FALSE)[1:(wk+1)]
  knnmat_beta0.sa[i, knnid] = 1
  #knnmat[knnid, i] = 1
}
diag(knnmat_beta0.sa) = rep(0, ntr)
library(pheatmap)
pheatmap(knnmat_beta0.sa,cluster_rows = FALSE,cluster_cols = FALSE)

knnmat.a = knnmat_beta0.sa
wnz.a = sum((knnmat.a + t(knnmat.a))!=0)/2 
Dmat.a = matrix(0, ncol=ntr, nrow=wnz.a)
## nrow of Dmat is determined by # of nonzero in W
rowIndex=1
for(i in 1:(ntr-1)){
  for(j in (i+1):ntr){
    if(knnmat.a[i,j] ==0){
    }else {
      if(knnmat.a[j,i] == 0){
        Dmat.a[rowIndex, i] = 1
        Dmat.a[rowIndex, j] = -1 
        rowIndex = rowIndex + 1
      }else{
        Dmat.a[rowIndex, i] = 2
        Dmat.a[rowIndex, j] = -2 
        rowIndex = rowIndex + 1
      }
    }
  }
}
for (i in 2:ntr) {
  for (j in 1:(i-1)) {
    if(knnmat.a[i,j] == 0){
    }else{
      if(knnmat.a[j,i] == 1){
      }else{
        Dmat.a[rowIndex, i] = -1
        Dmat.a[rowIndex, j] = 1 
        rowIndex = rowIndex + 1
      }
    }
  }
}
fit.a = APG_fit_L1(x = X.datac,y = y,beta0 = beta0.sa,lam1 = lam1,lam2 = lam2,
                  rho = 1,wk = k1,Dmat = Dmat.a,W = knnmat.a,L1multiplier = L1,
                  t = t,maxiter = 1000,nADMM = 500,B=B)


