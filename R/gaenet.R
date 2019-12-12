GAEnet <- function(mydata, penalty_id, adWeight, alpha, lambda, maxit, eps)
{
    # data info
    n <- dim(mydata[,, 1])[1]
    p <- dim(mydata[,, 1])[2] - 1
    D <- dim(mydata)[3]

    # standardize data
    stdMydata <- stdData.GAEnet(mydata)
    X <- stdMydata$x
    Y <- stdMydata$y
    meanx <- stdMydata$meanx
    normx <- stdMydata$normx

    # L1 and L2 penalty
    L1 <- lambda * adWeight * alpha * penalty_id
    L2 <- lambda * (1 - alpha) * penalty_id

    # coefficient update
    delta <- 1
    beta  <- matrix(1,   nrow = p, ncol = D)
    beta0 <- matrix(1,  nrow = 1, ncol = D)
    ita <- matrix(0, nrow = n, ncol = D)
    pi  <- matrix(0, nrow = n, ncol = D)
    res <- matrix(0, nrow = n, ncol = D)
    z <- matrix(0,   nrow = D, ncol = 1)
    Y.tilde <- Y

    it <- 0
    comp.set <- 1:p
    while (delta > eps && it < maxit)  {
        it <- it + 1
        beta.old  <- beta
        beta0.old <- beta0

        for (d in 1:D) {
            ita[, d] <- X[,, d] %*% beta[, d]
            ita[, d] <- ita[, d] + beta0[1, d]

            pi[, d] <- exp(ita[, d]) / (1 + exp(ita[, d]))
            pi[, d][abs(1 - pi[, d]) <= 1e-5] = 1
            pi[, d][abs(pi[, d]) <= 1e-5] <- 0

            Y.tilde[, 1, d] <- ita[, d] + 4 * (Y[, 1, d] - pi[, d])
            res[, d] = (Y.tilde[, 1, d] - ita[, d])
        }

        # update intercept, res
        for (d in 1:D) {
            beta0[1, d] <- mean(Y.tilde[, 1, d])
            res[, d] <- res[, d] - (beta0[1, d] - beta0.old[1, d])
        }

        # update b
        for (j in comp.set) {
            # update z
            for (d in 1:D)
                z[d, 1] <- t(X[, j, d]) %*% res[, d] + n * beta[j, d]

            #update beta
            beta[j,] <- t(threshold.gaenet(z / (4 * n), L1[j], L2[j]))

            #update res
                for (d in 1:D)
                    res[, d] <- res[, d] - X[, j, d] * (beta[j, d] - beta.old[j, d])
        }

        comp.set <- which(b[,1] != 0)
        delta    <- max(abs(beta - beta.old))
    }

    # naive coefficients
    beta.naive <- beta
    beta0.naive <- beta0

    # tranform back into scale
    coef.naive <- matrix(NA, nrow = p + 1, ncol = D)
    for (d in 1:D){
        mean_d_norm <- meanx[,, d] / normx[,, d]
        # intercept
        coef.naive[1, d] <- beta0.naive[, d] - mean_d_norm %*% beta.naive[1:p, d]
        coef.naive[2:(p + 1), d] <- beta.naive[1:p, d] / normx[,, d]
    }

    # coef
    beta  <- beta * (1 + 2 * L2[1])
    beta0 <- beta0

    # tranform back into scale
    coef = matrix(NA, nrow = (p+1), ncol = D)
    for (d in 1:D) {
        mean_d_norm <- meanx[,, d] / normx[,, d]
        # intercept
        coef[1, d] <- beta0[,d] - mean_d_norm %*% beta[1:p, d]
        coef[2:(p + 1), d] <- beta[1:p, d] / normx[,, d]
    }

    list(coef.naive = coef.naive, coef = coef)
}
