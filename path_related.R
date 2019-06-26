require(lars)
require(glmnet)
require(parcor)
require(parallel)
require(mvnfast)
#sourceCpp('LOCO_TS.cpp')


LOCO.TS = function(path_type = "lars", x, y, 
                   whichCov, betaNULL = 0, enet.control = list(),
                   s = 2, t = 2){
  if(path_type == "lars"){
    return( LOCOLars.TS.new(ExtractLars.Path(x, y, whichCov, betaNULL), s, t) )  
  }else if(path_type == "enet"){
    if(length(enet.control) == 0){
      enet.control = list(nlambda = 1000, 
                          alpha = 1, family = "gaussian")
    }
    return( LOCOLars.TS.new(ExtractNet.Path(x = x, y = y, 
                                            whichCov = whichCov, 
                                            betaNULL = betaNULL, 
                                            nlambda = enet.control$nlambda, 
                                            alpha = enet.control$alpha, 
                                            family = enet.control$family), 
                            s, t) ) 
  } 
  
}

ExtractLars.Path <- function(x, y, whichCov, betaNULL = 0){
  # Extract Path info
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, and the LOCO path
  #
  #
  n = nrow(x)
  p = ncol(x)
  ################# Here we go ! ! ! #######################################################  
  ##### extract the path #######
  
  if (length(whichCov) != length(betaNULL)){
    stop("Length of variables being tested must equal to # of hypothesis")
  }else{
    
    multiTest <- length(whichCov) > 1
    
  }
  
  # simultanoues testing 
  if(multiTest){
    
    #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
    adjust.X <- rowSums( t( t(x) * betaNULL  ) )    
    adjust.Y <- y - adjust.X #- mean(y)
    
    X_scaled <- scale(x)
    
    
    lars.out <- lars(X_scaled, adjust.Y, type = "lasso", intercept = FALSE, use.Gram = FALSE, normalize=TRUE)
    lambda.hat <- sort(c(lars.out$lambda,0),decreasing=FALSE)
    
    
    
    lars.j.out <- lars(X_scaled[, -whichCov], adjust.Y, type = "lasso", intercept = FALSE, use.Gram = FALSE,normalize=TRUE)
    lambda.j.hat <- sort(c(lars.j.out$lambda,0),decreasing=FALSE)
    
    
    
    
    # indivdual testing   
  }else if(!multiTest){  #indivdual test
    
    adjust.Y <- y - betaNULL * x[,whichCov] #- mean(y)
    X_scaled <- scale(x)
    
    lars.out <- lars(X_scaled, adjust.Y , type = "lasso", intercept = FALSE, use.Gram = FALSE, normalize=TRUE)
    lambda.hat <- sort(c(lars.out$lambda,0),decreasing = FALSE)
    
    
    lars.j.out <- lars(X_scaled[, -whichCov], adjust.Y, type = "lasso", intercept = FALSE, use.Gram = FALSE, normalize=TRUE)
    lambda.j.hat <- sort(c(lars.j.out$lambda,0), decreasing = FALSE)
    
    
  }else{
    stop("wrong input of multiTest, must be boolean")
  } 
  
  
  #### remove the unoverlapping part (test)
  
  #if(lambda.hat[1] != lambda.j.hat[1]){
  #  leftmost = c(lambda.hat[1], lambda.j.hat[1])
  #  whichone = which.max(leftmost)
  #  if(whichone ==1){
  #    lambda.j.hat = lambda.j.hat[lambda.j.hat >= lambda.hat[1]]
  #    lambda.j.hat = c(lambda.hat[1], lambda.j.hat)
  #  }else{
  #    lambda.hat = lambda.hat[lambda.hat >= lambda.j.hat[1]]
  #    lambda.hat = c(lambda.j.hat[1], lambda.hat)
  #  }
  #}
  
  #minLam = 0
  #lambda.hat = lambda.hat[lambda.hat > minLam]
  #lambda.j.hat = lambda.hat[lambda.hat > minLam]
  
  
  union.lambda <- sort(unique(c(lambda.hat,lambda.j.hat)), decreasing = FALSE)
  
  beta.j.hat <- matrix(0, length(union.lambda), p)
  
  beta.hat <- predict(lars.out, s=union.lambda, type="coef", mode="lambda")$coefficients
  
  beta.j.hat[, -whichCov] <- predict(lars.j.out, s=union.lambda, type="coef", mode="lambda")$coefficients
  
  #beta.hat <- coef(lars.out, lam = union.lambda)
  #beta_val <- coef(lars.j.out, lam = union.lambda)
  
  
  
  return(list(union.lambda = union.lambda, beta.hat = beta.hat, beta.j.hat = beta.j.hat)) 
  
}





LOCOLars.TS <- function(obj){
  # Calculate PATH statistic exactly  
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, the LOCO path, and Test Statistic
  #
  #  
  
  M <- length(obj$union.lambda)
  
  Delta <- obj$beta.j.hat - obj$beta.hat
  
  Delta_1 <- Delta[-M, ]
  
  Delta_2 <- Delta[-1, ]
  
  Lambda <- diff(obj$union.lambda)
  
  Epsilon <- 1/3 * Lambda * (Delta_1 * Delta_1 + Delta_1 * Delta_2 + Delta_2 * Delta_2)
  
  return(sum(Epsilon))
  
}




LOCOLars.TS.new <- function(obj, s = 2, t = 2){
  # Calculate PATH statistic exactly  
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, the LOCO path, and Test Statistic
  #
  #  
  
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







## NOW, GLMNET ! @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

#######################################################
#######################################################
################## GLMNET #############################
################## GLMNET #############################
################## GLMNET #############################
#######################################################
#######################################################



ExtractNet.Path <- function(x, y, whichCov, betaNULL = 0, nlambda = 1000,  alpha = 1, family = "gaussian"){
  # Extract Path info
  #
  # Args:
  # X,Y: design matrix and response vector
  # which.covariate: if is a vector, indicating which covariate we will be computing; if is a list: then do multiple testing.
  # 
  # normalize: argguments of lars 
  # betaNULL: same size and same data type with which.covariate, specify the null hypothesis H0: beta = betaNULL. 
  # Returns:
  # A list of lambda vector, original path, and the LOCO path
  #
  #
  n <- nrow(x)
  p <- ncol(x)
  
  ################# Here we go ! ! ! #######################################################  
  
  ##### extract the path #######
  
  if (length(whichCov) != length(betaNULL)){
    
    stop("Length of variables being tested must equal the length of their Null hypothesis")
    
  }else{
    
    multiTest <- length(whichCov) > 1
    
  }
  
  # simultanoues testing 
  if(multiTest){
    
    #### ajust Y by using Y - X_j * beta_j_NULL - X_k * beta_k_NULL - ...
    adjust.X <- rowSums( t( t(x) * betaNULL  ) )    
    adjust.Y <- as.vector( y - adjust.X - mean(y))
    X_scaled = scale(x)
    
    ## TODO
    max_lam = NA
    
    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -whichCov], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
  
    # indivdual testing   
  }else if(!multiTest){  #indivdual test
    
    adjust.Y <- as.vector( y - betaNULL * x[,whichCov] - mean(y) )
    X_scaled <- scale(x)
    
    ## TODO
    net.out <- glmnet(X_scaled, adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.hat <- sort(net.out$lambda, decreasing = FALSE)
    
    net.j.out <- glmnet(X_scaled[, -whichCov], adjust.Y, nlambda = nlambda, lambda.min.ratio = 0.001, alpha = alpha, family = family, intercept = FALSE, standardize = TRUE)
    lambda.j.hat <- sort(net.j.out$lambda, decreasing = FALSE)
    
  
  }else{
    stop("wrong input of multiTest, must be boolean")
  } 
  

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




