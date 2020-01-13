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
