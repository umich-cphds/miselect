Stacked and Grouped Penalized Regression For Multiply Imputed Data
==================================================================

Penalized regression methods, such as lasso and elastic net, are used in
many biomedical applications when simultaneous regression coefficient
estimation and variable selection is desired. However, missing data
complicates the implementation of these methods, particularly when
missingness is handled using multiple imputation. Applying a variable
selection algorithm on each imputed dataset will likely lead to
different sets of selected predictors, making it difficult to ascertain
a final active set without resorting to ad hoc combination rules.
'mianet' presents Stacked Adaptive Elestic Net (saenet) and Grouped
Adaptive LASSO (galasso) for continuous and binary outcomes, which by
construction force selection of the same variables across multiply
imputed data.

Installation
------------

`mianet` can installed from Github via

    # install.packages("devtools")
    devtools::install_github("umich-cphds/mianet", build_opts = c())

The Github version may contain bug fixes not yet present on CRAN, so if
you are experiencing issues, you may want to try the Github version of
the package.

Example
-------

TODO

Bugs
----

If you encounter a bug, please open an issue on the
[Issues](https://github.com/umich-cphds/mianet/issues) tab on Github or
send us an email.

Contact
-------

For questions or feedback, please email Jiacong Du at
<jiacong@umich.edu> or Alexander Rix <alexrix@umich.edu>.

<!-- ## References -->
