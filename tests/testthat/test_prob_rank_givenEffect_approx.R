
#===============================================================================
prob_rank_givenEffect_approx <- function(k, et, ey, nrep = 10000, m0, m1,
                                  effectType = c("binary", "continuous"))
{
    t <- rnorm(nrep, et, 1)
    p0 <- pnorm(-t, mean = 0, sd = 1, lower.tail = TRUE)

    if(effectType == "binary"){p1 <- pnorm(ey - t, mean = 0, sd = 1, lower.tail = TRUE)
    } else { if(ey == 0){p1 <- pnorm(ey - t, mean = 0, sd = 1, lower.tail = TRUE)
        } else {
            m = m0 + m1
            a = ey - 1
            b = ey
            xb = b - t
            xa = a - t
            p1 = (xb*pnorm(xb) - xa*pnorm(xa) + dnorm(xb) - dnorm(xa))/(b-a)
        }
    }

    mean0 <- (m0 - 1)*p0 + m1*p1 + 1
    mean1 <- m0*p0 + (m1 - 1)*p1 + 1
    var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
    var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)
    prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
                   mean(dnorm(k, mean1, sqrt(var1))))
    return(prob)
}




