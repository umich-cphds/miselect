S <- function(x, t) sign(x) * max((abs(x) - t), 0)

threshold.gaenet <- function(z, L1, L2)
{
    # spectral norm
    normz <- norm(z, type = "2")

    S(normz, L1) / (0.25 + 2 * L2) * (z / normz)
}

threshold.waenet <- function(z, z2, L1, L2)
{
    S(z, L1) / (z2 + 2 * L2)
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
