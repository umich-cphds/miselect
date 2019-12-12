cv.GAEnet.sol = function(mydata, L1, L2, maxiter, eps){
  
  n = dim(mydata[,,1])[1]
  p = dim(mydata[,,1])[2]-1
  D = dim(mydata)[3]
  
  # standardize data
  stdMydata = stdData.GAEnet(mydata)
  X = stdMydata$x
  Y = stdMydata$y
  meanx = stdMydata$meanx
  normx = stdMydata$normx
  
  # coefficient update 
  dif.b = 1
  iter = 0
  b = matrix(1, nrow = p, ncol = D)
  b0 = matrix(1, nrow = 1, ncol = D)
  Y.tilda = Y
  compSet = seq(1:p)
  ita = matrix(0, nrow = n, ncol = D)
  pi  = matrix(0 ,nrow = n, ncol = D)
  res = matrix(0, nrow = n, ncol = D)
  z = matrix(0, nrow = D, ncol = 1)
  
  while ((dif.b > eps) & (iter < maxiter)){
    iter = iter+1
    b.old = b
    b0.old = b0
    
    for (d in 1:D){
      ita[,d] = X[,,d] %*% b[,d]
      ita[,d] = ita[,d] +b0[1,d]
      pi[,d] = exp(ita[,d])/(1+exp(ita[,d]))
      pi[,d][abs(1-pi[,d])<=1e-5] = 1
      pi[,d][abs(pi[,d] <= 1e-5)] = 0
      Y.tilda[,1,d] = ita[,d] + (Y[,1,d] - pi[,d])/0.25
      res[,d] = (Y.tilda[,1,d] - ita[,d])
    }
    
    # update intercept
    for (d in 1:D){
      b0[1,d] = mean(Y.tilda[,1,d])
    }
    # update res
    for (d in 1:D){
      res[,d] = res[,d] - (b0[1,d]-b0.old[1,d])
    }
    
    # update b
    for (j in compSet){
      # update z
      for (d in 1:D){
        z[d,1] = (t(X[,j,d]) %*% res[,d] + n*b[j,d])
      }
      
      #update beta
      b[j,] = t(S_func.GAEnet(z = (0.25/n)*z, L1 = L1, L2 = L2, j = j))
      
      #update res
      for(d in 1:D){
        res[,d] = res[,d] - X[,j,d] * (b[j,d]-b.old[j,d])
      }
    }
    compSet = which(b[,1] != 0)
    dif.b = max(abs(b-b.old))
  }
  
  # naive coefficients
  b.naive = b
  b0.naive = b0
  
  # tranform back into scale
  coefficients.naive = matrix(NA, nrow = (p+1), ncol = D)
  for (d in 1:D){
    mean_d_norm = meanx[,,d]/normx[,,d]
    coefficients.naive[1,d] = b0.naive[,d] - mean_d_norm%*%b.naive[1:p,d] # intercept
    coefficients.naive[2:(p+1),d] = b.naive[1:p,d]/normx[,,d]
  }
  
  # coefficients
  b = b*(1+2*L2[1])
  b0= b0
  
  # tranform back into scale
  coefficients = matrix(NA, nrow = (p+1), ncol = D)
  for (d in 1:D){
    mean_d_norm = meanx[,,d]/normx[,,d]
    coefficients[1,d] = b0[,d] - mean_d_norm%*%b[1:p,d] # intercept
    coefficients[2:(p+1),d] = b[1:p,d]/normx[,,d]
  }
  
  
  return(list(coefficients.naive = coefficients.naive,
              coefficients = coefficients))
}