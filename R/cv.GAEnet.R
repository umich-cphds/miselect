cv.GAEnet = function(mydata, k, penalty_id, adWeight, alphaVec, n.lambda, maxiter, eps){
  
  n = dim(mydata[,,1])[1]
  p = dim(mydata[,,1])[2]-1
  D = dim(mydata)[3]

  # adaptive weight power
  v = log(p*D)/log(n*D)
  gamma = ceiling(2*v/(1-v))+1
  adWeight.power = (gamma+1)/2
  
  # calculate dev for seq alphaVec and lamVec
  dev = matrix(NA, nrow = length(alphaVec), ncol = n.lambda)
  lambda.max = lambda.max.GAEnet(mydata = mydata, adWeight = adWeight, penalty_id = penalty_id, alpha = alphaVec)$lambda.max
  if (all(adWeight == rep(1,p))){
    lamVec = exp(seq(log(lambda.max*1e-3), log(lambda.max), length.out = n.lambda))
  } else {
    lamVec = exp(seq(log(lambda.max*(n*D)^(-adWeight.power)), log(lambda.max*(n*D)^(-0.5)), length.out = n.lambda))
  }
  
  for (i in 1:length(alphaVec)) {
    alpha = alphaVec[i]
    for (j in 1:length(lamVec)){
      lambda = lamVec[j]
      cv.model = cv.GAEnet.dev(mydata = mydata,
                               k = k,
                               penalty_id = penalty_id,
                               adWeight = adWeight,
                               alpha = alpha,
                               lambda = lambda,
                               maxiter = 1000, 
                               eps = 1e-8)
      dev[i,j] = cv.model$dev
    }
  }
  
  # find dev.min
  dev.min = min(dev)
  opt.dev.min = which (dev == min(dev), arr.ind = TRUE)
  alpha.min = alphaVec[opt.dev.min[1]]
  lambda.min = lamVec[opt.dev.min[2]]
  
  # find dev.1se
  dev.min = min(dev)
  dev.sd = sd(dev)
  dev.1se = dev.min+dev.sd
  set.1se = which(dev <= dev.1se, arr.ind = TRUE)
  new.1se = cbind(set.1se, alphaVec[set.1se[,1]]*lamVec[set.1se[,2]])
  row.id = which(new.1se[,3] == max(new.1se[,3]))
  alpha.1se = alphaVec[set.1se[row.id,1]]
  lambda.1se = lamVec[set.1se[row.id,2]]
  if (length(alpha.1se)!= 1){
    opt = which(dev == min(diag(dev[set.1se[row.id,][,1],set.1se[row.id,][,2]])), arr.ind = TRUE)
    alpha.1se = alphaVec[opt[1]]
    lambda.1se = lamVec[opt[2]]
  }
  
  return(list(dev=dev,
              alpha.min = alpha.min,
              lambda.min = lambda.min,
              alpha.1se = alpha.1se,
              lambda.1se = lambda.1se))
}