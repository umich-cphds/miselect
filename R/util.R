# soft threshold function
S <- function(x, t)
{
    if (x > t)
        x - t
    else if (x < -t)
        x + t
    else
        0
}

# Overcomes a technical limitation of R that drops almost all attributes from
# objects when subsetted using `[`
subset_scaled_matrix <- function(x, rows)
{
    .x <- x[rows, , drop = F]
    attr(.x, "scaled:center") <- attr(x, "scaled:center")
    attr(.x, "scaled:scale")  <- attr(x, "scaled:scale")
    .x
}
