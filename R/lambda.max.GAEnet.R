lambda.max.GAEnet <- function(mydata, adWeight, penalty_id, alpha)
{
    n <- dim(mydata)[1]
    p <- dim(mydata)[2] - 1
    D <- dim(mydata)[3]

    stdMydata <- stdData.GAEnet(mydata)
    X <- stdMydata$x
    Y <- stdMydata$y
    meanx <- stdMydata$meanx
    normx <- stdMydata$normx

    Y_X = matrix(NA, D, p)
    for (d in 1:D)
    Y_X[d,] <- t(Y[, 1, d]) %*% X[,, d]

    norm <- sqrt(apply(Y_X ^ 2, 2, sum))

    lambda.seq <- norm / (n * alpha * adWeight)
    lambda.max <- max(lambda.seq * penalty_id)

    return(list(lambda.max = lambda.max))
}
