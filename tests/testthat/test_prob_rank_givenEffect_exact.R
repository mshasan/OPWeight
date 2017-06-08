
#===============================================================================
prob_rank_givenEffect_exact <- function(k, et, ey, nrep = 10000, m0, m1,
                                        effectType = c("binary", "continuous"))
{
    k0 <- 1:k
    fun.k0 <- function(k0)
    {
        t <- rnorm(nrep, et, 1)
        p0 <- pnorm(-t)

        if(effectType == "binary"){p1 <- pnorm(ey - t)
        } else { if(ey == 0){p1 <- pnorm(ey - t)
            } else {
                m = m0 + m1
                a = ey - 1
                b = ey
                xb = b - t
                xa = a - t
                p1 = (xb*pnorm(xb) - xa*pnorm(xa) + dnorm(xb) - dnorm(xa))/(b-a)
            }
        }

        E.T <- ifelse(et == 0, mean(dbinom(k0-1, m0-1, p0)*dbinom(k-k0, m1, p1)),
                      mean(dbinom(k0-1, m0, p0)*dbinom(k-k0, m1-1, p1)))
        return(E.T)
    }
    prob <- sum(sapply(k0,fun.k0))
    return(prob)
}
