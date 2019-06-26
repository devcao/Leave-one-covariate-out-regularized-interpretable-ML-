require(lars)
require(glmnet)


extract_lars_null <- function(x, y){
  n = nrow(x); p = ncol(x)
  
  X_scaled <- scale(x)
  y = y - mean(y)
  
  lars.out <- lars(X_scaled, y, type = "lasso", 
                   intercept = FALSE, use.Gram = FALSE, normalize=TRUE)
  lambda.hat <- sort(c(lars.out$lambda,0),decreasing = FALSE)
  
  return(list(y = y,
              X_scaled = X_scaled,
              lars.out = lars.out,
              lambda.hat = lambda.hat,
              path = "lars",
              n = n,
              p = p))
}

extract_lars_alter <- function(lars_null, whichCov){
  
  X_scaled = lars_null$X_scaled; y = lars_null$y;
  lambda.hat = lars_null$lambda.hat; lars.out = lars_null$lars.out;
  n = lars_null$n; p = lars_null$p
  
  lars.j.out <- lars(X_scaled[, -whichCov], y, type = "lasso", 
                       intercept = FALSE, use.Gram = FALSE, normalize=TRUE)
  lambda.j.hat <- sort(c(lars.j.out$lambda,0), decreasing = FALSE)
    
  ## remove the unoverlapping part (test)
  
  if(lambda.hat[1] != lambda.j.hat[1]){
    leftmost = c(lambda.hat[1], lambda.j.hat[1])
    whichone = which.max(leftmost)
    if(whichone ==1){
      lambda.j.hat = lambda.j.hat[lambda.j.hat >= lambda.hat[1]]
      lambda.j.hat = c(lambda.hat[1], lambda.j.hat)
    }else{
      lambda.hat = lambda.hat[lambda.hat >= lambda.j.hat[1]]
      lambda.hat = c(lambda.j.hat[1], lambda.hat)
    }
  }
  
  union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
  beta.j.hat <- matrix(0, length(union.lambda), p)
  beta.hat <- predict(lars.out, s=union.lambda, type="coef", mode="lambda")$coefficients
  beta.j.hat[, -whichCov] <- predict(lars.j.out, s=union.lambda, type="coef", mode="lambda")$coefficients

  #beta.hat <- coef(lars.out, lam = union.lambda)
  #beta_val <- coef(lars.j.out, lam = union.lambda)

  return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}


extract_path = function(x, y, path_type){
  p = ncol(x);
  if(path_type == "lars"){
    lars_null = extract_lars_null(x, y, intercept = FALSE)
    out = lapply(1:p, extract_lars_alter, lars_null = lars_null)
  }else if(path_type == "enet"){
    enet_null = extract_enet_null(x, y, intercept = FALSE)
    out = lapply(1:p, extract_enet_alter, enet_null = enet_null)
  }
  
  return(out)
}



extract_loco = function(path_null_obj, whichCov, s = 2, t = 2){
  path_type = path_null_obj$path
  
  if(path_type == "lars"){
    out = extract_lars_alter(lars_null = path_null_obj, whichCov = whichCov)
    varimp = locoVarImp_Engine(out, s = s, t = t)
  }else if(path_type == "enet"){
    out = extract_enet_alter(enet_null = path_null_obj, whichCov = whichCov)
    varimp = locoVarImp_Engine(out, s = s, t = t)
  }
  
  return(varimp)
}


loco_pipeline = function(x, y, path_type, s = 2, t = 2, n_threads = -1, ...){
  p = ncol(x)
  if(path_type == "lars"){
    path_null = extract_lars_null(x, y)
  }else if(path_type == "enet"){
    path_null = extract_enet_null(x, y, ...)
  }
  
  if(is.null(n_threads)){
    out = unlist(lapply(1:p, extract_loco, path_null_obj = path_null, s = s, t = t))
  }else{
    if(n_threads == -1){n_threads = detectCores()}
    cl = makeCluster(n_threads, type = "FORK")
    out = unlist(parLapply(cl, 1:p, extract_loco, path_null_obj = path_null, s = s, t = t))
    stopCluster(cl)
  }
  return(out)
}



  
locoVarImp_Engine <- function(obj, s = 2, t = 2){

  M <- length(obj$union.lambda)
  Delta <- obj$beta.j.hat - obj$beta.hat
  Delta_1 <- Delta[-M, ]
  Delta_2 <- Delta[-1, ]
  Lambda <- diff(obj$union.lambda)
  equal_sign_indicator <- Delta_2 * Delta_1 < 0
  Epsilon <- 1/(s+1) * Lambda * abs( (Delta_2 ^ (s+1) - ((-1) ^ equal_sign_indicator) * Delta_1 ^ (s+1)) / (Delta_2 - Delta_1) )

  if(s == t){
    return( (sum(Epsilon, na.rm = TRUE)) ^ (1/t) )
  }else{
    return( sum(rowSums(Epsilon) ^ (t/s), na.rm = TRUE) ^ (1/t) )
  }
  
}


locoVarImp <- function(path_obj, s = 2, t = 2){
  return( unlist(lapply(path_obj,locoVarImp_Engine, s = s, t = t)) )
}



extract_enet_null <- function(x, y, nlambda = 1000, alpha = 1, 
                              family = "gaussian"){
  n = nrow(x); p = ncol(x)
  X_scaled <- scale(x)
  y = y - mean(y)
  net.out <- glmnet(X_scaled, y, 
                    nlambda = nlambda, 
                    lambda.min.ratio = 0.001, alpha = alpha, 
                    family = family, standardize = TRUE, intercept = FALSE)
  lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
  
  return(list(y = y,
              X_scaled = X_scaled,
              net.out = net.out,
              lambda.hat = lambda.hat,
              alpha = alpha,
              nlambda = nlambda,
              family = family,
              path = "enet",
              n = n,
              p = p))
}


extract_enet_alter = function(enet_null, whichCov){
  
  X_scaled = enet_null$X_scaled; y = enet_null$y;
  lambda.hat = enet_null$lambda.hat; net.out = enet_null$net.out
  n = lars_null$n; p = lars_null$p
    
  net.j.out = glmnet(x = X_scaled[, -whichCov], y, 
                     nlambda = enet_null$nlambda, 
                     lambda.min.ratio = 0.001, alpha = enet_null$alpha, 
                     family = enet_null$family, intercept = FALSE, standardize = TRUE)
  
  lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
  minLam <- min(lambda.hat, lambda.j.hat)
  
  lambda.hat = lambda.hat[lambda.hat > minLam]
  
  lambda.j.hat = lambda.hat[lambda.hat > minLam]
  
  union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
  
  beta.j.hat <- Matrix(0, length(union.lambda), p)  # sparse Matrix
  
  beta.hat <- predict(net.out, s=union.lambda, type="coef")
  beta.hat <- t( beta.hat[-1, ] )  # we don't need intercept
  
  
  beta.j.tmp <- predict(net.j.out, s=union.lambda, type="coef")
  beta.j.hat[, -whichCov] <- t( beta.j.tmp[-1, ] ) # we don't need intercep
  
  
  return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}


