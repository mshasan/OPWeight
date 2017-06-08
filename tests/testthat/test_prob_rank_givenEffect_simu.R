
#===============================================================================
prob_rank_givenEffect_simu <- function(s, ey, e.one, m0, m1,
                                    effectType = c("binary", "continuous"))
    {
        m = m0 + m1
        ey0 <- rep(0, m0)

        if(effectType == "binary"){ey1 <- rep(ey, m1)
            } else {if(ey == 0){ey1 <- rep(0, m1)
                } else {ey1 <- runif(m1, ey-1, ey)
                    }
                }
        Ey <- c(ey1, ey0)
        Ey[1] <- e.one
        t01 <- rnorm(m, Ey, 1)
        r1 <- rank(-t01)[1]
        r0 <- rank(-t01)[m1+1]
        cbind(r0, r1)
    }


