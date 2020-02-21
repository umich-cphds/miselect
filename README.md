<!-- Badges -->
[![CRAN
Version](https://img.shields.io/cran/v/miselect?style=flat-square&color=blue&label=CRAN)](https://cran.r-project.org/package=miselect)
[![GitHub
Release](https://img.shields.io/github/v/release/umich-cphds/miselect?include_prereleases&label=Github&style=flat-square&color=blue)](https://github.com/umich-cphds/miselect)
[![Travis
CI](https://img.shields.io/travis/umich-cphds/miselect?style=flat-square)](https://travis-ci.org/umich-cphds/mianet)

Variable Selection for Multiply Imputed Data
============================================

Penalized regression methods, such as lasso and elastic net, are used in
many biomedical applications when simultaneous regression coefficient
estimation and variable selection is desired. However, missing data
complicates the implementation of these methods, particularly when
missingness is handled using multiple imputation. Applying a variable
selection algorithm on each imputed dataset will likely lead to
different sets of selected predictors, making it difficult to ascertain
a final active set without resorting to ad hoc combination rules.
'miselect' presents Stacked Adaptive Elastic Net (saenet) and Grouped
Adaptive LASSO (galasso) for continuous and binary outcomes. They, by
construction, force selection of the same variables across multiply
imputed data.

Installation
------------

`miselect` can installed from Github via

    # install.packages("devtools")
    devtools::install_github("umich-cphds/miselect", build_opts = c())

The Github version may contain bug fixes not yet present on CRAN, so if
you are experiencing issues, you may want to try the Github version of
the package.

Example
-------

Here is how to use cross validated `saenet`. A nice feature of `saenet`
is that you can cross validate over `alpha` without having to use
`foldid`.

    library(miselect)
    library(mice)
    #> Loading required package: lattice
    #> 
    #> Attaching package: 'mice'
    #> The following objects are masked from 'package:base':
    #> 
    #>     cbind, rbind

    set.seed(48109)
    # Using the mice defaults for sake of example only.
    mids <- mice(miselect.df, m = 5, printFlag = FALSE)
    dfs <- lapply(1:5, function(i) complete(mids, action = i))

    # Generate list of imputed design matrices and imputed responses
    x <- list()
    y <- list()
    for (i in 1:5) {
        x[[i]] <- as.matrix(dfs[[i]][, paste0("X", 1:20)])
        y[[i]] <- dfs[[i]]$Y
    }

    # Calculate observational weights
    weights  <- 1 - rowMeans(is.na(miselect.df))
    pf       <- rep(1, 20)
    adWeight <- rep(1, 20)
    alpha    <- c(.5 , 1)

    # Since 'Y' is a binary variable, we use 'family = "binomial"'
    fit <- cv.saenet(x, y, pf, adWeight, weights, family = "binomial",
                     alpha = alpha)

    # By default 'coef' returns the betas for (lambda.min , alpha.min)
    coef(fit)
    #>  (Intercept)           X1           X2           X3           X4           X5 
    #>  0.075336069  1.465012227  0.773562992 -0.132892848  1.814095899  0.090605078 
    #>           X6           X7           X8           X9          X10          X11 
    #>  0.000000000  1.803555612  0.002901045 -0.053939183 -0.090949859  1.314651703 
    #>          X12          X13          X14          X15          X16          X17 
    #> -0.001398250  0.000000000  0.000000000 -0.210676071  0.057351962  0.340355076 
    #>          X18          X19          X20 
    #>  0.071108301 -0.269948553  0.000000000

You can supply different values of `lambda` and `alpha`. Here we use the
`lambda` and `alpha` selected by the one standard error rule

    coef(fit, lambda = fit$lambda.1se, alpha = fit$alpha.1se)
    #> (Intercept)          X1          X2          X3          X4          X5 
    #>  0.08216228  1.09429258  0.41249142  0.00000000  1.29535531  0.00000000 
    #>          X6          X7          X8          X9         X10         X11 
    #>  0.00000000  1.32200685  0.00000000  0.00000000  0.00000000  0.93730624 
    #>         X12         X13         X14         X15         X16         X17 
    #>  0.00000000  0.00000000  0.00000000 -0.01425428  0.00000000  0.11562370 
    #>         X18         X19         X20 
    #>  0.00000000 -0.11133637  0.00000000

Bugs
----

If you encounter a bug, please open an issue on the
[Issues](https://github.com/umich-cphds/miselect/issues) tab on Github
or send us an email.

Contact
-------

For questions or feedback, please email Jiacong Du at
<jiacong@umich.edu> or Alexander Rix <alexrix@umich.edu>.

<!-- ## References -->
