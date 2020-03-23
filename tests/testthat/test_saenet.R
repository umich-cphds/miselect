context("Test saenet")


test_that("saenet works", {
mids <- mice(miselect.df, m = 5, printFlag = FALSE)
dfs <- lapply(1:5, function(i) complete(mids, action = i))

x <- list()
y <- list()
for (i in 1:5) {
    x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
    y[[i]] <- dfs[[i]]$Y
}

weights  <- 1 - rowMeans(is.na(miselect.df))
pf       <- rep(1, 20)
adWeight <- rep(1, 20)


expect_silent({
    fit <- saenet(x, y, pf, adWeight, weights, family = "binomial")})
})


test_that("cv.saenet works", {
mids <- mice(miselect.df, m = 5, printFlag = FALSE)
dfs <- lapply(1:5, function(i) complete(mids, action = i))

x <- list()
y <- list()
for (i in 1:5) {
    x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
    y[[i]] <- dfs[[i]]$Y
}

weights  <- 1 - rowMeans(is.na(miselect.df))
pf       <- rep(1, 20)
adWeight <- rep(1, 20)


expect_silent({
    fit <- cv.saenet(x, y, pf, adWeight, weights, family = "binomial")})
})
