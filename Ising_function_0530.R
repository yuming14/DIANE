library(parallel)
library(IsingSampler)
library(dplyr)
# library(RGCCA)
library(caret)
library(lava)
library(Matrix)

wkd = "~/Ising_Model/code/2302/"
exp_p = 10

Simu_Theta_lr_block = function(p = 5, r = 1){
  r = min(r, p)
  U = matrix(rnorm(p*r), nrow = p) / sqrt(r*p)
  Theta = U%*%t(U)
  diag(Theta) = 0
  return(list(U = U, Theta = Theta))
}

Simu_Theta_lr = function(p = 5, r = 1, blocksize = 5){
  if(length(p)>1){
    if(length(r)==1) r = rep(r, length(p))
    Theta_list = lapply(1:length(p), function(i) return(Simu_Theta_lr(p[i], r[i])))
    Theta = lapply(Theta_list, function(x) x$Theta)
    U = do.call("bdiag", lapply(Theta_list, function(x) x$U))
    return(list(U = U, Theta = Theta))
  }
  n1 = p%/%blocksize
  p2 = p%%blocksize
  if(n1>=1){
    Theta_list = lapply(1:n1, function(o){
      return(Simu_Theta_lr_block(blocksize, r))
    })
    if(p2>1){
      Theta_list[[n1+1]] = Simu_Theta_lr_block(p2, r)
    }
    Theta = lapply(Theta_list, function(x) x$Theta)
    U = Reduce(bdiag, lapply(Theta_list, function(x) x$U))
    return(list(U = U, Theta = Theta))
  }else{
    Theta_list = Simu_Theta_lr_block(p2, r)
    Theta = list(U = Theta_list$U, Theta = list(Theta_list$Theta))
  }
  return(Theta)
}

Simu_Theta_block = function(p = 5, prob = c(0,1)){
  Theta = matrix(sample(0:1, p^2, TRUE, prob = prob), p, p) * rnorm(p^2)
  Theta = pmax(Theta, t(Theta)) / p
  diag(Theta) = 0
  return(Theta)
}

Simu_Theta = function(p = 5, blocksize = 5, prob = c(0, 1)){
  if(length(p)>1){
    return(lapply(p, function(p0) return(Simu_Theta_block(p0, prob))))
  }
  n1 = p%/%blocksize
  p2 = p%%blocksize
  if(n1>=1){
    Theta = lapply(1:n1, function(o){
      return(Simu_Theta_block(blocksize, prob))
    })
    if(p2>1){
      Theta[[n1+1]] = Simu_Theta_block(p2, prob)
    }
  }else{
    Theta = list(Simu_Theta_block(p2, prob))
  }
  return(Theta)
}

Simu_x = function(n = 1000, Theta, Beta = 1, Resp = c(-1L,1L), nIter = 100,
                  method = "MH"){
  # cl = makeCluster(detectCores()-1)
  # clusterExport(cl, c("wkd","n","Beta","nIter","Resp","method"), 
  #               envir = environment())
  x = lapply(Theta, function(theta){
    return(IsingSampler(n, theta, thresholds = rep(0, nrow(theta)), beta = Beta, 
                        nIter = nIter, responses = Resp,
                        method = method))
  })
  # stopCluster(cl)
  x = do.call("cbind", x)
  return(x)
}

sigmoid = function(x){
  return(1/(1+exp(-x)))
}

get_grad_Theta_new = function(x, Theta){
  C = 2 * (x %*% (Theta-diag(diag(Theta))))
  grad_Theta = - (t(x+1) %*% x - 2 * t(sigmoid(C)) %*% x) / nrow(x)
  grad_Theta = grad_Theta + t(grad_Theta)
  return(grad_Theta)
}

get_grad_Theta_slow = function(x, Theta){
  n = nrow(x)
  p = nrow(Theta)
  grad = matrix(0, nrow = p, ncol = p)
  for(i in 1:p){
    theta = matrix(0, nrow = p, ncol = p)
    theta[i,] = Theta[i,]
    theta[,i] = Theta[,i]
    for(l in 1:n){
      C = 1/(exp(x[l,,drop = FALSE]%*%theta%*%t(x[l,,drop=FALSE]))+1)
      tmp = t(x[l,,drop=FALSE])%*%x[l,,drop=FALSE] * c(C)
      grad[i,] = grad[i,] + tmp[i,]
      grad[,i] = grad[,i] + tmp[,i]
    }
  }
  return(-2*grad/n)
}

get_grad_Theta = function(x, Theta, ZeroDiag = TRUE){
  # given Theta(i,j) = Theta(j,i)
  stopifnot(ncol(x)==ncol(Theta))
  stopifnot(nrow(Theta)==ncol(Theta))
  if(ZeroDiag) diag(Theta) = 0
  C = x %*% (Theta + t(Theta)) - x * diag(Theta)
  C = C * x
  x_tilde = x / (exp(C)+1)
  grad_Theta = t(x_tilde) %*% x
  grad_Theta = -(grad_Theta + t(grad_Theta))/nrow(x)*2
  return(grad_Theta)
}

Sing_Value_ST = function(A, lambda, PSD){
  eigen.x = eigen(A)
  if(PSD){
    eigen.value = pmax(eigen.x$values - lambda, 0)
  }else{
    eigen.value = pmax(abs(eigen.x$values) - lambda, 0) * sign(eigen.x$values)
  }
  idx = which(eigen.value>0)
  return(list(Theta = eigen.x$vectors %*% 
                diag(eigen.value) %*% 
                t(eigen.x$vectors),
              U = eigen.x$vectors[,idx] %*% 
                diag(sqrt(eigen.value[idx]))))
}

Sing_Value_ST_old = function(A, lambda = 1, PSD = TRUE){
  svd_A = svd(A)
  if(PSD){
    sig = svd_A$d * (2*(sign(svd_A$u[1,])==sign(svd_A$v[1,]))-1)-lambda
  }else{
    sig = svd_A$d-lambda
  }
  idx = which(sig>0)
  sig = sig[idx]
  if(length(idx)==0){
    return(list(Theta = matrix(0, nrow = nrow(A), ncol = ncol(A)),
                U = matrix(0, nrow = nrow(A), ncol = 1)))
  }else if(length(idx)==1){
    return(list(Theta = t(t(svd_A$u[,idx]))%*%t(svd_A$v[,idx]) * sig,
                U = matrix(svd_A$u[,idx], ncol = 1)*sqrt(abs(sig))))
  }else{
    return(list(Theta = svd_A$u[,idx] %*% diag(sig) %*% t(svd_A$v[,idx]),
                U = svd_A$u[,idx] %*% diag(sqrt(abs(sig)))))
  }
}

est_Theta_convex = function(x, Theta0 = NULL, C = NULL,
                            lambda = 1, eta = 1, 
                            epsilon = 1e-2, maxstep = 100,
                            diag_add = 0.5, PSD = TRUE,
                            Theta_star = NULL, info = TRUE,
                            ZeroDiag = TRUE,
                            timing = TRUE, file = NULL){
  if(timing){
    compute_time = numeric(maxstep)
    t0 = Sys.time()
  } 
  if(is.null(Theta0)){
    Theta_old = matrix(0, nrow = ncol(x), ncol = ncol(x))
  }else{
    Theta_old = Theta0
  }
  
  delta = numeric(maxstep)
  if(!is.null(Theta_star)) err = numeric(maxstep)
  for(step in 1:maxstep){
    grad_Theta = get_grad_Theta(x, Theta_old, ZeroDiag)
    if(is.null(C)){
      A = Theta_old - eta * grad_Theta
    }else{
      A = Theta_old - eta * grad_Theta - eta * C
    }
    diag(A) = diag(Theta_old) - eta * diag_add
    Theta_new = Sing_Value_ST(A = A, lambda = lambda, PSD = PSD)
    U = Theta_new$U
    Theta_new = Theta_new$Theta
    if(timing) compute_time[step] = difftime(Sys.time(), t0, units = "secs")
    tmp = Theta_new-Theta_old
    diag(tmp) = 0
    delta[step] = norm(tmp, 'F')
    if(!is.null(Theta_star)){
      tmp = Theta_new-Theta_star
      diag(tmp) = 0
      err[step] = norm(tmp, 'F')
    }
    if(delta[step]<epsilon) break
    Theta_old = Theta_new
    if(!is.null(file)){
      save(Theta_new, U, delta, step, file = file)
    }
  }
  delta = delta[1:step]
  if(timing) compute_time = compute_time[1:step]
  if(!is.null(Theta_star)) err = err[1:step]
  if(step==maxstep & info) cat("not converge!","\n")
  if(is.null(Theta_star)){
    if(timing){
      return(list(Theta = Theta_new, delta = delta, err = NA,
                  U = U, compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta, err = NA))
    }
  }else{
    if(timing){
      return(list(Theta = Theta_new, delta = delta, err = err,
                  U = U, compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta, err = err,
                  U = U))
    }
  }
}

est_Theta_DC_convex = function(x, Theta_star = NULL, lambda = 1, eta = 1,
                               m = 10, diag_add = 0.5, PSD = TRUE,
                               maxstep = c(100,100), epsilon = c(1e-2,1e-2),
                               ZeroDiag = TRUE,
                               info = TRUE, timing = TRUE){
  if(timing){
    compute_time = numeric(maxstep[1])
    t0 = Sys.time()
  }
  Theta_old = matrix(0, nrow = ncol(x), ncol = ncol(x))
  idxm = createFolds(1:nrow(x), k = m)
  delta = numeric(maxstep[1])
  if(!is.null(Theta_star)) err = numeric(maxstep[1])
  deltabase = NULL
  for(t in 1:maxstep[1]){
    for(i in 1:length(idxm)){
      grad_L_m = get_grad_Theta(x[idxm[[i]],], Theta_old, ZeroDiag)
      if(i==1){
        C = grad_L_m * (1/m - 1)
      }else{
        C = C + grad_L_m
      }
    }
    Theta_new = est_Theta_convex(x = x[idxm[[1]],], Theta0 = Theta_old, 
                                 C = C, lambda = lambda, eta = eta, 
                                 epsilon = epsilon[2], maxstep = maxstep[2],
                                 diag_add = diag_add, PSD = PSD,
                                 ZeroDiag = TRUE,
                                 info = FALSE, timing = FALSE)$Theta
    if(timing) compute_time[t] = difftime(Sys.time(), t0, units = "secs")
    tmp = Theta_new - Theta_old
    diag(tmp) = 0
    delta[t] = norm(tmp, 'F')
    if(!is.null(Theta_star)){
      tmp = Theta_new - Theta_star
      diag(tmp) = 0
      err[t] = norm(tmp, 'F')
    }
    if(delta[t] >= exp_p*delta[1]){
      if(info){ cat("explosure","\n") }
      break
    }
    Theta_old = Theta_new
    if(delta[t] < epsilon[1]) break
  }
  if(timing) compute_time = compute_time[1:t]
  delta = delta[1:t]
  if(t==maxstep[1] & info) cat("not converge!","\n")
  if(is.null(Theta_star)){
    if(timing){
      return(list(Theta = Theta_new, delta = delta, err = NA,
                  compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta, err = NA))
    }
  }else{
    if(timing){
      return(list(Theta = Theta_new, delta = delta, err = err[1:t],
                  compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta, err = err[1:t]))
    }
  }
}

get_L_grad_Theta = function(x, Theta){
  diag(Theta) = 0
  stopifnot(ncol(x)==ncol(Theta))
  stopifnot(nrow(Theta)==ncol(Theta))
  C = x %*% (Theta + t(Theta)) - x * diag(Theta)
  C = C * x
  x_tilde = x / (exp(C)+1)
  grad_Theta = t(x_tilde) %*% x
  grad_Theta = -(grad_Theta + t(grad_Theta))/nrow(x)
  return(grad_Theta)
}

get_Q_grad_UV = function(U, V){
  tmp = t(U)%*%U - t(V)%*%V
  return(list(grad_U = U%*%tmp, grad_V = V%*%(-tmp)))
}

get_UV_Theta = function(Theta, d){
  d = min(d, ncol(Theta), nrow(Theta))
  svd_Theta = svd(Theta, nu = d, nv = d)
  if(d==1){
    return(list(U = matrix(svd_Theta$u,ncol=1) * sqrt(svd_Theta$d[1]),
                V = matrix(svd_Theta$v,ncol=1) * sqrt(svd_Theta$d[1])))
  }else{
    return(list(U = svd_Theta$u %*% diag(sqrt(svd_Theta$d[1:d])),
                V = svd_Theta$v %*% diag(sqrt(svd_Theta$d[1:d]))))
  }
}

est_Theta_nonconvex = function(x, Theta0 = NULL, U = NULL, V = NULL,
                               C = NULL, eta = 0.1, d = 100, 
                               maxstep = 100, epsilon = 1e-2, 
                               diag_add = 0.5, PSD = TRUE,
                               Theta_star = NULL, info = FALSE,
                               lambda = 0.01, maxstep_convex = 10, 
                               timing = TRUE, stopsign = "Theta"){
  if(timing){
    compute_time = numeric(maxstep)
    t0 = Sys.time()
  }
  if(is.null(Theta0)){
    Theta_old = est_Theta_convex(x, lambda = lambda, eta = eta,
                                 epsilon = epsilon[1], maxstep = maxstep_convex, 
                                 diag_add = diag_add, PSD = PSD,
                                 Theta_star = NULL,
                                 ZeroDiag = TRUE,
                                 info = FALSE, timing = FALSE)$Theta
  }else{
    Theta_old = Theta0
  }
  if(is.null(U) | is.null(V)){
    UV = get_UV_Theta(Theta_old, d)
    U = UV$U
    V = UV$V
  }
  
  delta = numeric(maxstep)
  if(!is.null(Theta_star)) err = numeric(maxstep)
  for(step in 1:maxstep){
    gL = get_L_grad_Theta(x, U %*% t(V))
    diag(gL) = 0
    gQ = get_Q_grad_UV(U, V)
    dU = gL %*% V + gQ$grad_U
    dV = gL %*% U + gQ$grad_V
    if(!is.null(C)){
      dU = dU + C %*% V
      dV = dV + t(C) %*% U
    }
    U = U - eta * dU
    V = V - eta * dV
    Theta = U %*% t(V)
    if(timing){
      compute_time[step] = difftime(Sys.time(), t0, units = "secs")
    }
    tmp = Theta - Theta_old
    diag(tmp) = 0
    if(stopsign == "U"){
      delta[step] = norm(dU, 'F')*eta
    }else{
      delta[step] = norm(tmp, 'F')
    }
    if(!is.null(Theta_star)){
      errtmp = Theta-Theta_star
      diag(errtmp) = 0
      err[step] = norm(errtmp, 'F')
    }
    if(delta[step]<epsilon)
      break
    if(delta[step]>=delta[1]*exp_p){
      if(info)
        cat("explosure!","\n")
      break
    }
    Theta_old = Theta
  }
  if(step==maxstep & info)
    cat("not converge!")
  delta = delta[1:step]
  if(timing) compute_time = compute_time[1:step]
  if(is.null(Theta_star)){
    if(timing){
      return(list(Theta = Theta, U = U, V = V, 
                  delta = delta,
                  compute_time = compute_time))
    }else{
      return(list(Theta = Theta, U = U, V = V, 
                  delta = delta))
    }
  }else{
    err = err[1:step]
    if(timing){
      return(list(Theta = Theta, U = U, V = V, 
                  delta = delta, 
                  err = err, compute_time = compute_time))
    }else{
      return(list(Theta = Theta, U = U, V = V, 
                  delta = delta,
                  err = err))
    }
  }
  
}

est_Theta_DC_nonconvex = function(x, Theta0 = NULL, m = 10, eta = 0.1, d = 100, 
                                  maxstep = c(100,100), epsilon = c(1e-2,1e-2),
                                  diag_add = 0.5, PSD = TRUE,
                                  Theta_star = NULL, info = FALSE,
                                  lambda = 0.01, maxstep_convex = 100,
                                  ini_convex = 1, timing = TRUE,
                                  stopsign = "Theta",
                                  file = NULL){
  if(timing){
    compute_time = numeric(maxstep[1])
    t0 = Sys.time()
  }
  idxm = createFolds(1:nrow(x), k = m)
  if(is.null(Theta0)){
    if(ini_convex == 1){
      for(i in 1:length(idxm)){
        tmp = est_Theta_convex(x[idxm[[i]],], lambda = lambda, eta = eta,
                               epsilon = epsilon[1], maxstep = maxstep_convex, 
                               diag_add = diag_add, PSD = PSD, Theta_star = NULL,
                               ZeroDiag = TRUE,
                               info = FALSE)$Theta
        if(i==1){
          Theta_old = tmp
        }else{
          Theta_old = Theta_old + tmp
        }
      }
      Theta_old = Theta_old/length(idxm)
    }else{
      Theta_old = est_Theta_convex(x[idxm[[1]],], lambda = lambda, eta = eta,
                                   epsilon = epsilon[1], maxstep = maxstep_convex, 
                                   diag_add = diag_add, PSD = PSD, Theta_star = NULL,
                                   ZeroDiag = TRUE,
                                   info = FALSE)$Theta
    }
  }else{
    Theta_old = Theta0
  }
  delta = numeric(maxstep[1])
  if(!is.null(Theta_star)) err = numeric(maxstep[1])
  deltabase = NULL
  U = NULL; V = NULL
  for(t in 1:maxstep[1]){
    for(i in 1:length(idxm)){
      tmp = get_L_grad_Theta(x[idxm[[i]],], Theta_old)
      if(i==1){
        C = tmp * (1/m - 1)
      }else{
        C = C + tmp/m
      }
    }
    tmp = est_Theta_nonconvex(x = x[idxm[[1]],], Theta0 = Theta_old, 
                              U = U, V = V, C = C,
                              eta = eta, d = d, maxstep = maxstep[2],
                              diag_add = diag_add, PSD = PSD,
                              epsilon = epsilon[2], timing = FALSE,
                              stopsign = stopsign)
    if(timing) compute_time[t] = difftime(Sys.time(), t0, units = "secs")
    Theta_new = tmp$Theta
    U = tmp$U
    V = tmp$V
    tmp = Theta_new - Theta_old
    diag(tmp) = 0
    delta[t] = norm(tmp, 'F')
    if(!is.null(Theta_star)){
      tmp = Theta_new - Theta_star
      diag(tmp) = 0
      err[t] = norm(tmp, 'F')
    }
    if(delta[t] >= exp_p*delta[1]){
      if(info){ cat("explosure","\n") }
      break
    }
    Theta_old = Theta_new
    if(delta[t] < epsilon[1]) break
    if(!is.null(file)) save(Theta_new, U, V, delta, t, file = file)
  }
  if(t==maxstep[1] & info) cat("not converge!")
  delta = delta[1:t]
  if(timing) compute_time = compute_time[1:t]
  if(is.null(Theta_star)){
    if(timing){
      return(list(Theta = Theta_new, delta = delta,
                  U = U, V = V,
                  compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta))
    }
  }else{
    err = err[1:t]
    if(timing){
      return(list(Theta = Theta_new, delta = delta, err = err,
                  U = U, V = V,
                  compute_time = compute_time))
    }else{
      return(list(Theta = Theta_new, delta = delta, 
                  U = U, V = V, err = err))
    }
  }
}
