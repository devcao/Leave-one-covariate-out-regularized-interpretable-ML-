
depenDesign = function(n, p, beta, rho){
  if(rho == 'equl'){  # equi corr
    Sigma = matrix(rep(0.8,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu,sigma = Sigma)
    Y <- X %*% beta + rnorm(n,0,1)   
  }else if (rho == "weak_equl"){
    Sigma = matrix(rep(0.5,p*p),p,p)
    diag(Sigma) = rep(1,p)
    Mu=rep(0,p)
    X=rmvn(n, mu = Mu, sigma = Sigma)
    Y <- X %*% beta + rnorm(n,0,1)   
  }else if (rho > 0){  # toeplitz matrix
    Sigma = toeplitz( rho^(0:(p - 1)) )
    Mu = rep(0,p)
    #X = rmvnorm(n,mean = Mu,sigma = Sigma)
    X = rmvn(n, mu = Mu,sigma = Sigma)
    Y <- X %*% beta + rnorm(n,0,1)  
  }else if (rho == 0){  # independent
    Sigma = diag(rep(1,p))
    Mu=rep(0,p)
    X <- matrix(rnorm(n*p), n, p)
    Y <- X %*% beta + rnorm(n,0,1)
  }
  return(list(X = X, Y = Y))   
  
}

power.loco.test = function(n, p, beta, rho, 
                           iter = 500, B = 500,
                           whichCov, betaNULL){
  if(n > p){
    
    obj = (lapply(1:iter, power.loco.engine, 
                  n = n, p = p, beta = beta, rho = rho,
                  B = B, whichCov = whichCov, betaNULL = betaNULL))
    loco_pval = unlist(lapply(obj, function(obj){obj$loco_obj$pval}))
    
    lm_pval = unlist(lapply(obj, function(obj){obj$lm_pval}))
    
    
    loco_power = c(); lm_power = c();
    loco_power[1] = mean(loco_pval <= 0.2) 
    loco_power[2] = mean(loco_pval <= 0.1)
    loco_power[3] = mean(loco_pval <= 0.05)
    loco_power[4] = mean(loco_pval <= 0.01)
    
    lm_power[1] = mean(lm_pval <= 0.2) 
    lm_power[2] = mean(lm_pval <= 0.1)
    lm_power[3] = mean(lm_pval <= 0.05)
    lm_power[4] = mean(lm_pval <= 0.01)
    
    return(list(lm_power = lm_power,
                lm_pval = lm_pval,
                loco_power = loco_power, 
                loco_pval = loco_pval, 
                obj = obj))
    
  }else{
  obj = (lapply(1:iter, power.loco.engine, 
                       n = n, p = p, beta = beta, rho = rho,
                       B = B, whichCov = whichCov, betaNULL = betaNULL))
  pval = unlist(lapply(obj, function(obj){obj$loco_obj$pval}))
  power = c()
  power[1] = mean(pval <= 0.2) 
  power[2] = mean(pval <= 0.1)
  power[3] = mean(pval <= 0.05)
  power[4] = mean(pval <= 0.01)
  
  return(list(power = power, pval = pval, obj = obj))
  }
}

power.loco.engine = function(ind, n, p, beta, rho, 
                             B = 500,
                             whichCov, 
                             betaNULL){
  data = depenDesign(n = n, p = p, beta = beta, rho = rho)
  if(n>p){
    lm_pval = (summary(lm(y~., 
                         data = data.frame(y = data$Y, x = data$X)))$coefficients)[whichCov+1, 4] 
  }else{
    lm_pval = NULL
  }
  return(list(lm_pval = lm_pval, 
              loco_obj = loco_resample(path_type = "lars", 
                         x = data$X, y = data$Y,
                         whichCov = whichCov, B = B, 
                         betaNULL = betaNULL, 
                         n_threads = -1)))
  
}



loco.power.curve = function(n, p, s, rho, 
                            iter=500, B=500, 
                            whichCov = 1, betaNULL = 0){
 
  #betaList = seq(-1,1, by = 0.2)
  betaList = 0
  
  pwr = matrix(NA, length(betaList), 4)
  i = 1
  for(b in betaList){
    
    beta = c(rep(1,s),rep(0,p-s))
    beta[1] = b
    
    pwr[i, ] = power.loco.test(n = n, p = p, beta = beta, rho = rho,
                               iter = iter, B = B, 
                               whichCov = 1, betaNULL = 0)$power
    cat("Round: ",i,"Now runing for: ", b, "power is: ", pwr[i,], "\n")
    
    i = i + 1
    
    
  }
  return(pwr)
  
}




#pwr_ld_indep = loco.power.curve(n=100, p = 100, s = 3, rho = 0)
#pwr_ld_wk = loco.power.curve(n=100, p = 100, s = 3, rho = "weak_equl")
#pwr_ld_eq = loco.power.curve(n=100, p = 100, s = 3, rho = "equl")
#pwr_ld_ar09 = loco.power.curve(n=100, p = 100, s = 3, rho = 0.9)
#pwr_ld_ar05 = loco.power.curve(n=100, p = 100, s = 3, rho = 0.5)






data = depenDesign(n=100, p=12, beta=c(0,1,1,rep(0,9)), 0)
x = data$X; y = data$Y

all_result = power.loco.test(n=100, p=12, beta=c(0,1,1,rep(0,9)), rho=0, 
                           iter = 100, B = 500,
                           whichCov=1, betaNULL=0)

all_result$loco_power

all_result$lm_power

#loco_resample(path_type = "lars", x, y, 1, B = 500, betaNULL = 0, n_threads = -1)

#power.loco.engine(ind = 1,n=100, p=12, beta=c(0,1,1,rep(0,9)), 
 #                 rho=0, which=1,betaNULL=0)

