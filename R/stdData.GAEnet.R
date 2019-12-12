stdData.GAEnet = function(mydata)
{
  n <- dim(mydata[,, 1])[1]
  p <- dim(mydata[,, 1])[2] - 1
  D <- dim(mydata)[3]
  x <- array(0, dim = c(n, p, D))
  y <- array(0, dim = c(n, 1, D))
  meanx <- array(0, dim = c(1, p, D))
  normx <- array(1, dim = c(1, p, D))
  for (d in 1:D) {
    x[,, d] <- mydata[,1:p,d]
    meanx[,, d] <- apply(x[,,d], 2, mean)
    x[,, d] <- scale(x[,,d], center = meanx[,, d], FALSE)
    normx[,, d] <- sqrt(apply(x[,, d] ^ 2, 2, sum) / n)
    x[,, d]   <- scale(x[,, d], center = FALSE, scale = normx[,, d])
    y[, 1, d] <- mydata[, p + 1, d]
  }
  return(list(x = x, y = y, meanx = meanx, normx = normx))
}
