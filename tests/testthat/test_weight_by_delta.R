
#===============================================================================

# function to compute weight for delta, the lagrange constant
#---------------------------------------------------------------

weight_by_delta <- function(delta, alpha = .05, et, m, m1, tail = 1L, ranksProb,
                            effectType = c("continuous", "binary"))
    {
        if(effectType == "continuous"){
            weight_per_delta <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta/(alpha*ranksProb)),
                                                   lower.tail = FALSE)
        } else {
            weight_per_delta <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta*m/(alpha*m1*ranksProb)),
                                 lower.tail = FALSE)
        }

        sumWeight_per_delta <- sum(weight_per_delta, na.rm = TRUE)

        return(sumWeight_per_delta)
    }



