
loco_resample = function(path_type = "lars", x, y, 
                         whichCov, betaNULL = 0, 
                         s = 2, t = 2,
                         enet.control = list(),
                         B = 500, n_threads = -1,
                         plot_null = FALSE){
  
  n = nrow(x); p = ncol(x)
  TS = LOCO.TS(path_type = path_type, x = x, y = y, 
               whichCov = whichCov, betaNULL = betaNULL, 
               enet.control = enet.control,
               s = s, t = t)
  
  rej = c()
    
  if(p >= n){ 
    bhat = adalasso(X = x, y = y, k = 10, use.Gram = FALSE, 
                      both = TRUE, intercept = FALSE)$coefficients.adalasso
  }else{
    bhat = ginv(t(x) %*% x) %*% t(x) %*% y
      #lm.fit(scale(x),y-mean(y))$coefficient#lm.fit(x, y)$coefficients
  }
  residual = y - x%*%bhat
  b.Null = bhat; b.Null[whichCov] = betaNULL; 
  
  if(is.null(n_threads)){
    TS_null=unlist(lapply(X=1:B, loco_boot_engine, 
                          path_type = path_type,
                          x_mat = x, b.Null = b.Null, 
                          residual = residual, betaNULL = betaNULL, 
                          whichCov = whichCov, enet.control = enet.control,
                          s = s, t = t))
  }else{
    if(n_threads == -1){n_threads = detectCores()}
    cl = makeCluster(n_threads, type = "FORK")
    TS_null=unlist(parLapply(cl, X=1:B, loco_boot_engine, 
                             path_type = path_type,
                             x_mat = x, b.Null = b.Null, 
                             residual = residual, betaNULL = betaNULL, 
                             whichCov = whichCov, enet.control = enet.control,
                             s = s, t = t))
    stopCluster(cl)
  }  
  
  rej[1] = TS > quantile(TS_null,0.8)
  rej[2] = TS > quantile(TS_null,0.9)
  rej[3] = TS > quantile(TS_null,0.95)
  rej[4] = TS > quantile(TS_null,0.99)
  pval = mean(TS_null > TS)
  if(plot_null){hist(TS_null, probability = TRUE, xlim = c(0, max(TS_null, TS))); abline(v = TS)}
  return(list(TS = TS, TS_null = TS_null, rej = rej, pval = pval))
}


  
loco_boot_engine = function(ind, path_type, 
                            x_mat, b.Null, residual, 
                            whichCov, betaNULL, 
                            enet.control,
                            s, t){
  
  n = nrow(x_mat)
  ind = sample(1:n,replace = TRUE)
  boot_residual = residual[ind]
  y = x_mat %*% b.Null + boot_residual
  
  return(LOCO.TS(path_type = path_type,
                 x = x_mat, y = y, 
                 whichCov = whichCov, enet.control = enet.control,
                 betaNULL = betaNULL, s = s, t = t))
} 
  
loco_boot_run = function(x, y, whichCov, betaNULL, multiTest, 
                            path_tyep = "enet", B = 500, n_threads = -1){
  
  if(is.null(n_threads)){
    out = unlist(lapply(1:p, extract_loco, path_null_obj = path_null, s = s, t = t))
  }else{
    if(n_threads == -1){n_threads = detectCores()}
    cl = makeCluster(n_threads, type = "FORK")
    out = unlist(parLapply(cl, 1:p, extract_loco, path_null_obj = path_null, s = s, t = t))
    stopCluster(cl)
  }
  
  
}
