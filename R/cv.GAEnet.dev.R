cv.GAEnet.dev = function(mydata, k, penalty_id, adWeight, alpha, lambda, maxiter, eps){
  
  n = dim(mydata[,,1])[1]
  p = dim(mydata[,,1])[2]-1
  D = dim(mydata)[3]
  
  # set up foldid  
  size = round(n/k)
  foldid.order = NULL
  for (i in 1:(k-1)){
    foldid.order = c(foldid.order,rep(i,size))
  }
  foldid.order = c(foldid.order, rep(k,(n-(k-1)*size)))
  set.seed(1)
  foldid = sample(foldid.order)
  
  # L1 and L2 penalty
  L1 = lambda*adWeight*alpha*penalty_id
  L2 = lambda*(1-alpha)*penalty_id
  
  # start cross validation
  dev = NA
  dev.ifold = matrix(0,k,1)
  coef0.train = matrix(NA, nrow = 1, ncol = D)
  coef.train = matrix(NA, nrow = p, ncol = D)
  logL = matrix(NA, nrow = 1, ncol = D)
  dev.test = matrix(NA, nrow = 1, ncol = D)
  
  ## iterate k folds
  for (i in 1:k){
    sel = which(foldid == i)
    train = mydata[-sel,,]
    test = mydata[sel,,]
    
    model.train = cv.GAEnet.sol(train,
                                L1 = L1, 
                                L2 = L2, 
                                maxiter = 1000, 
                                eps = 1e-8)
    coef.train[,] = model.train$coefficients[-1,]
    coef0.train[1,] = model.train$coefficients[1,]
    
    for (d in 1:D){
      xb = test[,1:p,d]%*%coef.train[,d]
      xb = xb+coef0.train[1,d]
      logL[,d] = mean(test[,p+1,d]*(xb) - log(1 + exp(xb)))
      dev.test[,d] = -2*(logL[,d]) #-2*(logL_model - logL_saturate), and logL_saturate=log(1)=0
    }
    dev.ifold[i,1] = mean(dev.test)
  }
  dev = mean(dev.ifold)
  
  # output
  return(list(dev=dev))
}