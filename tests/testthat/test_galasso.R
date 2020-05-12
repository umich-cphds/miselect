context("Test galasso")


test_that("galasso works", {
library(mice)
mids <- mice(miselect.df, m = 5, printFlag = FALSE)
dfs <- lapply(1:5, function(i) complete(mids, action = i))

x <- list()
y <- list()
for (i in 1:5) {
    x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
    y[[i]] <- dfs[[i]]$Y
}

pf       <- rep(1, 20)
adWeight <- rep(1, 20)


expect_silent({
    fit <- galasso(x, y, pf, adWeight, nlambda = 50)})
})


test_that("cv.galasso works", {
library(mice)
mids <- mice(miselect.df, m = 5, printFlag = FALSE)
dfs <- lapply(1:5, function(i) complete(mids, action = i))

x <- list()
y <- list()
for (i in 1:5) {
    x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
    y[[i]] <- dfs[[i]]$Y
}

pf       <- c(0, rep(1, 19))
adWeight <- c(0, rep(1, 19))


expect_silent({
    fit <- cv.galasso(x, y, pf, adWeight, nlambda = 50, family = "binomial")})
})
