# x: n-by-p design matrix
# y: n-by-1 response matrix
# beta0: initial estimate n-by-p
# lam1:penalty tuning parameter for pairwise
# lam2:penalty tuning parameter for variable selection
# D:在程序中生成
# rho: penalty parameter
# nADMM: maximum number of iterations for ADMM, 2000 is default.
# admmAbsTol :  absolute tolerance for ADMM, 1e-04 is default.
# L1multiplier: 步长调节,初始步长应稍大一些
# t: 步长收缩幅度，t<1
# maxiter: 步长收缩最大步
# alpha: els net
library(ncvreg)
library(MASS)
library(glmnet)
library(SIS)
library(fmrs)
library(e1071)
library(FGSG)
APG_fit_L1 = function(x, y, beta0, lam1, lam2, rho, wk, nADMM = 2000,admmAbsTol = 1e-04,
                      Dmat,W,L1multiplier,t,maxiter,B){
  
  
  n = nrow(x)
  p = ncol(x)
  #list = createDmat(X = x,y = y,method="adapt", wk=wk, beta0_n = beta0)
  #D = list$Dmat
  D = Dmat
  
  #根据初始值生成的邻接矩阵
  W = W
  #pairwise的初始权重
  w.pair = update_w(beta0,W,B)
  
  nD = nrow(D)
  options(digits = 4)
  beta = beta0
  beta_hat = beta0
  #fx关于beta的导数
  gradfx = matrix(0,n,p)
  
  V = D %*% beta0
  
  Lambda = matrix(0,nD,p)
  
  #变量选择的初始权重
  zeta = update_zeta(beta0)
  #run = numeric(nADMM)
  
  #s.all = numeric(nADMM)
  #ADMM循环迭代
  for (i in 1:nADMM){
    #start = Sys.time()
    
    #李普希兹常数
    s = L1multiplier/eigenvalues(x)
    #生成中间变量β.t和V.t
    beta.trans = beta
    Lambda.trans = Lambda
    beta_hat = beta
    #fx关于β的导数更新
    gradfx = update_gradf(x,y,beta_hat)
    
    u = beta_hat - s*(gradfx + rho*t(D) %*% (D %*% beta_hat - V + Lambda))
    
    #更新beta 
    beta = update_prox_11(u = u,lambda = (lam2*s),weight = zeta)
    g.beta = phi_lam(x, y, beta, rho, D, V, Lambda)
    delta.g = t(gradfx + rho*t(D) %*% (D %*% beta_hat - V + Lambda))%*% (beta_hat - beta)
    g.condition = phi_lam(x,y,beta_hat,rho, D, V, Lambda) - 
      sum(diag(delta.g)) + norm((beta_hat - beta),type = 'f')^2/(2*s)
    for (j in 1:maxiter) {
      if(g.beta > g.condition){
        s = t*s
        u = beta_hat - s*(gradfx + rho*t(D) %*% (D %*% beta_hat - V + Lambda))
        beta = update_prox_11(u = u,lambda = (lam2*s),weight = zeta)
        g.beta = phi_lam(x, y, beta, rho, D, V, Lambda)
        delta.g = t(gradfx + rho*t(D) %*% (D %*% beta_hat - V + Lambda))%*% (beta_hat - beta)
        g.condition = phi_lam(x,y,beta_hat,rho, D, V, Lambda) - 
          sum(diag(delta.g)) + norm((beta_hat - beta),type = 'f')^2/(2*s)
      }
      else{
        break
      }
    }
    #s.all[i] = s
    #if(j == maxiter) {print("WARNING: maxium iteration")}
    
    #加速步
    #tau.t = tau
    #tau = 0.5* (1 + sqrt(1 + 4*tau^2))
    #更新加速步
    #beta_hat = beta + (tau.t - 1)/tau * (beta - beta.trans)
    
    
    v = D %*% beta + Lambda
    #更新V
    V = update_prox_12(u = v,lambda = (lam1/rho), weight = w.pair)
    
    #更新Lambda
    Lambda = Lambda + D %*% beta - V
    #Lambda 的加速
    #Lambda_hat = Lambda + (tau.t - 1)/tau (Lambda - Lambda.trans)
    
    
    #停止条件
    r1 = norm(beta - beta.trans,type = "f")
    r2 = norm(Lambda - Lambda.trans,type = "f")
    if(r1 < admmAbsTol && r2 < admmAbsTol){
      cat("Convergence !", "\n")
      print(i)
      break
    }
    #更新权重
    zeta = update_zeta(beta)
    w.pair = update_w(beta,knnmat = W,B)
    
    
    #w.pair = update_w(beta,W,B = sqrt(p))
    
    #w.pair = update_w_scad(beta,W,B,C)
    #end = Sys.time()
    #runningtime = end - start
    #run[i] = runningtime
  }
  #if(i == nADMM) {print("WARNING: maxium iteration achieves")}
  #return(list(beta = beta,s = run))
  return(beta)
}

#生成D矩阵
#生成pairwise的权wij
#beta0_n: n-by-p matrix
createDmat <- function(X, y, method="adapt", wk){
  ntr = nrow(X)
  p = ncol(X)
  ### Use adaptive method to calculate w, w=i/(bi-bj)
  if(method=="adapt"){
    #beta0_n=beta_local_n
    #beta0_mat = matrix(beta0_n, nrow=ntr, byrow=TRUE)
    distXY = as.matrix(dist(cbind(X,y))) # smallest wk
    knnmat <- matrix(0, ntr, ntr)# 1/0 matrix for knn
    for(i in 1:ntr){
      knnmati = distXY[i, ]
      knnid = order(knnmati, decreasing=FALSE)[1:(wk+1)]
      knnmat[i, knnid] = 1
      #knnmat[knnid, i] = 1
    }
  }
  diag(knnmat) = rep(0, ntr)
  k.upper = knnmat# diag useless
  k.upper[!upper.tri(k.upper,diag = FALSE)] = 0
  #########################################
  wnz = sum(k.upper!=0) 
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
          e = numeric(ntr)
          e[i] = -1
          e[j] = 1
          Dmat = rbind(Dmat,e)
        }
      }
    }
  }
  #nonZeroRow = which(apply(Dmat, 1, FUN=function(x) sum(x!=0)>0))
  #DmatNonZero = Dmat[nonZeroRow, ]; dim(DmatNonZero)
  #return(DmatNonZero)
  return(list(Dmat = Dmat, knn = knnmat))
}

#求损失值
phi_lam = function(x,y,beta,rho,D,V,Lambda){
  n = nrow(x)
  p = ncol(x)
  loss1 = 0
  for (i in 1:n) {
    loss1 = loss1 + (y[i] - x[i,] %*% beta[i,])^2
  }
  loss2 = 0.5*rho*(norm((D %*% beta - V + Lambda),type = 'f'))^2
  loss = loss1 + loss2
  return(loss)
}

#生成变量选择的权
update_zeta = function(beta){
  n = nrow(beta)
  p = ncol(beta)
  zeta = matrix(0,n,p)
  for (i in 1:n) {
    for (j in 1:p) {
      zeta[i,j] = 1/abs(beta[i,j])
    }
  }
  return(zeta)
}


#更新pairwise的权
update_w = function(beta,knnmat,B){
  ntr = nrow(beta)
  beta_diff = as.matrix(dist(beta, method = "euclidean"))
  diag(beta_diff) = rep(1, ntr) ## bi-bj matrix, diag useless
  W = knnmat/pmax(beta_diff,1e-03)
  wnz = sum((knnmat + t(knnmat))!=0)/2 ## 2*number of rows for Dmat
  #knn的权重
  w = numeric(wnz)
  rowIndex=1
  for(i in 1:(ntr-1)){
    for(j in (i+1):ntr){
      if(knnmat[i,j] ==0){
      }else {
        wij_bdd = min(W[i, j],B)
        w[rowIndex] = wij_bdd
        rowIndex = rowIndex + 1
      }
    }
  }
  for (i in 2:ntr) {
    for (j in 1:(i-1)) {
      if(knnmat[i,j] ==0){
      }else {
        if(knnmat[j,i] == 0){
          wij_bdd = min(W[i, j],B)
          w[rowIndex] = wij_bdd
          rowIndex = rowIndex + 1
        }
      }
    }
  }
  return(w)
}


#求李普希兹常数
eigenvalues = function(x){
  n = nrow(x)
  p = ncol(x)
  x.new = matrix(0,n,n*p)
  for (i in 1:n) {
    a = (i-1)*p + 1
    b = i*p
    x.new[i,a:b] = x[i,]
  }
  s = max(eigen(x.new %*% t(x.new))$val)
  return(s)
}

#求fx关于β的导数
update_gradf = function(x,y,beta){
  n = nrow(x)
  p = ncol(x)
  gradfx = matrix(0,n,p)
  for(i in 1:n){
    for(j in 1:p){
      gradfx[i,j] = x[i,j]*(y[i] - x[i,] %*% beta[i,])
    }
  }
  return(-2*gradfx)
}

#取正函数
takeP = function(x){
  if(x > 0){
    return(x)
  }else{
    return(0)
  } 
}

#符号函数
sgn = function(x){
  if(x>0){
    return(1)
  }else{
    if(x<0){
      return(-1)
    }
    else{
      return(0)
    }
  }
}


#L_11范数的prox求解
update_prox_11 = function(u,lambda,weight){
  n = nrow(u)
  p = ncol(u)
  for(j in 1:p){
    for(i in 1:n){
      u[i,j] = sgn(u[i,j]) * takeP(abs(u[i,j]) - lambda*weight[i,j])
    }
  }
  return(u)
}

#L_12范数的prox求解
update_prox_12 = function(u,lambda,weight){
  n = nrow(u)
  p = ncol(u)
  for(i in 1:n){
    m = sqrt(sum(u[i,]^2))
    x = 1 - lambda*weight[i]/m
    
    for(j in 1:p){
      u[i,j] = u[i,j] * takeP(x)
    }
  }
  return(u)
}
#初始值估计
get_local_init_beta2 <- function(X, y, eps=2.5, kn=30){
  ## include y information into distance cal
  ## when finding the neighboring for each point, if the x's are similar
  ## but the y is quite different, then it shold not be treated as neighbr
  ## dist(i,j) = |xi-xj|
  X = as.matrix(X); y = as.matrix(y)
  n = nrow(X);p = ncol(X)
  #X,Y的距离矩阵
  dist_mat_X = as.matrix(dist(X,method="manhattan"))
  dist_mat_y = as.matrix(dist(y,method="manhattan"))
  is_greater_than_eps = ifelse(dist_mat_y<eps, 1, Inf)
  dist_mat = dist_mat_X * is_greater_than_eps
  beta_local_init_mat = matrix(NA, n, p)
  for(i in 1:n){
    kntmp = min(kn+1, n - sum(dist_mat[i, ] == Inf, na.rm = T))
    
    if(kntmp>10){
      subID = order(dist_mat[i, ])[1:kntmp]
      Xsub = X[subID, ]; ysub = y[subID]
      Xsub = scale(Xsub,scale = F)
      ysub = scale(ysub,scale = F)
      #Xsub = scale(Xsub,center = TRUE,scale = FALSE)
      #ysub = scale(ysub,center = TRUE,scale = FALSE)
      lambda.rw1 = cv.ncvreg(Xsub,ysub,penalty = "SCAD",nfolds = 5)$lambda.min
      fit = ncvfit(Xsub,ysub,penalty = "SCAD",lambda = lambda.rw1)
      #lambda.rw1 = cv.glmnet(Xsub,ysub,nfolds = 3)$lambda.min 
      #fit = glmnet(Xsub,ysub,family = "gaussian",alpha = 1,intercept = FALSE,lambda = lambda.rw1)
      #beta_local_init_mat[i, ] = coef(fit)[-1]
      beta_local_init_mat[i, ] = fit$beta
    }else{
      beta_local_init_mat[i, ] = rep(0, p)
    }
  }
  return(beta_local_init_mat)
}

#初始值估计
get_local_init_beta_spare <- function(X.sis,X, y, eps=2.5, kn=30){
  ## include y information into distance cal
  ## when finding the neighboring for each point, if the x's are similar
  ## but the y is quite different, then it shold not be treated as neighbr
  ## dist(i,j) = |xi-xj|
  X = as.matrix(X); y = as.matrix(y)
  n = nrow(X);p = ncol(X)
  #X,Y的距离矩阵
  dist_mat_X = as.matrix(dist(X.sis,method="manhattan"))
  dist_mat_y = as.matrix(dist(y,method="manhattan"))
  is_greater_than_eps = ifelse(dist_mat_y<eps, 1, Inf)
  dist_mat = dist_mat_X * is_greater_than_eps
  beta_local_init_mat = matrix(NA, n, p)
  for(i in 1:n){
    kntmp = min(kn+1, n - sum(dist_mat[i, ] == Inf, na.rm = T))
    
    if(kntmp>10){
      subID = order(dist_mat[i, ])[1:kntmp]
      Xsub = X[subID, ]; ysub = y[subID]
      Xsub = scale(Xsub,scale = F)
      ysub = scale(ysub,scale = F)
      #Xsub = scale(Xsub,center = TRUE,scale = FALSE)
      #ysub = scale(ysub,center = TRUE,scale = FALSE)
      lambda.rw1 = cv.ncvreg(Xsub,ysub,penalty = "SCAD",nfolds = 3)$lambda.min
      fit = ncvfit(Xsub,ysub,penalty = "SCAD",lambda = lambda.rw1)
      #lambda.rw1 = cv.glmnet(Xsub,ysub,nfolds = 3)$lambda.min 
      #fit = glmnet(Xsub,ysub,family = "gaussian",alpha = 1,intercept = FALSE,lambda = lambda.rw1)
      #beta_local_init_mat[i, ] = coef(fit)[-1]
      beta_local_init_mat[i, ] = fit$beta
    }else{
      beta_local_init_mat[i, ] = rep(0, p)
    }
  }
  return(beta_local_init_mat)
}

# 一般情形利用岭回归求初值
initial_ridge = function(x,y,lam){
  n = nrow(x)
  p = ncol(x)
  npair = n*(n-1)/2
  ## change the matrix x to diag tilde.x = diag(x1,...,xn) = x %*% (I_p,...,I_p)
  x.new = matrix(0,n,n*p)
  for (i in 1:n) {
    a = (i-1)*p + 1
    b = i*p
    x.new[i,a:b] = x[i,]
  }
  ## generate matrix D and A 
  D <- NULL
  for (j in 1:(n - 1)) {
    D[[j]] <- matrix(0, nrow = n - j, ncol = n)
    D[[j]][, j] <- 1
    for (k in (j + 1):n) {
      D[[j]][k - j, k] <- -1
    }
  }
  D <- do.call(rbind, D)
  beta.vec = solve(t(x.new) %*% x.new + lam *kronecker(t(D) %*% D,diag(1,p))) %*% t(x.new) %*% y
  beta = matrix(data = beta.vec,n,p,byrow = TRUE)
  return(beta)
}


initial_ridge_homo = function(x,y,lam){
  #x = scale(x)
  n = nrow(x)
  p = ncol(x)
  beta = solve(t(x) %*% x + lam*diag(1,p)) %*% t(x) %*% y
  return(beta)
}

# 高维情形利用岭回归求初值
initial_ridge_hd = function(x,y,lam1,lam2){
  n = nrow(x)
  p = ncol(x)
  npair = n*(n-1)/2
  ## change the matrix x to diag tilde.x = diag(x1,...,xn) = x %*% (I_p,...,I_p)
  x.new = matrix(0,n,n*p)
  for (i in 1:n) {
    a = (i-1)*p + 1
    b = i*p
    x.new[i,a:b] = x[i,]
  }
  ## generate matrix D and A 
  D <- NULL
  for (j in 1:(n - 1)) {
    D[[j]] <- matrix(0, nrow = n - j, ncol = n)
    D[[j]][, j] <- 1
    for (k in (j + 1):n) {
      D[[j]][k - j, k] <- -1
    }
  }
  D <- do.call(rbind, D)
  #t(D) %*% D = n*diag(1,n) + matrix(1,n,n)
  A = lam1 * (n*diag(1,n) + matrix(1,n,n)) + 2*lam2 * diag(1,n)
  B = solve(A)
  #C = kronecker(B,diag(1,p))
  #C = kronecker(solve(A),diag(1,p))
  #B = diag(1,n) + x.new %*% kronecker(solve(A),diag(1,p)) %*% t(x.new) 
  #beta.vec = C -  C %*% t(x.new) %*% solve(B) %*% x.new %*% C %*% t(x.new) %*% y
  beta.vec = (kronecker(B,diag(1,p)) -  kronecker(B,diag(1,p)) %*% t(x.new) %*% solve(diag(1,n) + x.new %*% kronecker(B,diag(1,p)) %*% t(x.new)) %*% x.new %*% kronecker(B,diag(1,p))) %*% t(x.new) %*% y
  beta = matrix(data = beta.vec,n,p,byrow = TRUE)
  return(beta)
}


loss_fit = function(x,y,group){
  k = nrow(group$centers)
  loss = 0
  for (i in 1:k) {
    fit_sub = which(group$cluster == i)
    coef_fit_sub = group$centers[i,]
    error = y[fit_sub] - x[fit_sub,] %*% coef_fit_sub
    loss = loss + t(error) %*% error
  }
  return(loss)
}


BIC_est = function(x,y,loss,beta,k){
  n = nrow(x)
  p = ncol(x)
  C_n = log(log(n*p))
  group = kmeans(beta,centers = k)
  nonzero = numeric()
  for (i in 1:k) {
    center = group$centers[i,]
    center_trunc = ifelse(abs(center) < 0.2,0,center)
    nonzero = append(nonzero,which(center_trunc != 0))
  }
  q = length(unique(nonzero))
  #q = length(which(beta != 0))
  BIC = log(loss/n) + C_n * log(n) * q/n
  return(list(BIC = BIC,q = q))
}