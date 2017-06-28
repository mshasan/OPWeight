#' @title Find sum of weights for the LaGrange multiplier
#'
#' @description Compute sum of weights for a given value of the LaGrange multiplier
#'
#' @param delta Numeric value of the LagRange multiplier
#' @param alpha Numeric, significance level of the hypothesis test
#' @param et Numeric, mean effect size of the test statistics
#' @param m Integer, totoal number of hypothesis test
#' @param m1 Integer, number of true alternative tests
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' default is right-tailed test.
#' @param ranksProb Numeric vector of the ranks probability of the filter
#' statistics given the effect size
#' @param effectType Character ("continuous" or "binary"), type of effect sizes
#'
#' @details
#' To obtain the normalized weight, and to make sure that the sum of the weights is
#' equal to the number of tests and the weights are positive, an optimal value of the
#' LaGrange multiplier \code{delta} needed. This function will compute the weights for a given
#' value of the LaGrange multiplier and provide the sum of the weights in return.
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#'
#' @export
#'
#' @return \code{sumWeight_per_delta} sum of weights per delta value
#'
#' @examples
#'
#' # generate a sequence of delta
#' delta <- seq(0, 1, .0001)
#'
#' # compute probability fiven effect
#' filters = runif(100, min = 0, max = 2.5)
#' probs <- dnorm(filters, mean = 0, sd = 1)
#'
#' # compute the sum of weights for each delta
#' weightSum_by_delta <- sapply(delta, weight_by_delta, m = 100, m1 = 50, et = 2,
#'                             ranksProb = probs, effectType = "continuous")
#'
#===============================================================================
# function to compute weight from p(rank=k|filterEffect=ey)
#------------------------------------------------
# Input:-----
# delta = value of the LagRange multiplier
# alpha = significance level of the hypothesis test
# et = mean effect size of the test statistics
# m = totoal number of hypothesis test
# m1 = number of true alternative tests
# tail = right-tailed or two-tailed hypothesis test. default is right-tailed test
# ranksProb = probability of the filter statistics given the effect size
# effectType = type of effect sizes; c("continuous", "binary")

# internal parameters:-----
# weight_per_delta = weight vector per delta value

# output:-----
# sumWeight_per_delta = sum of weights per delta value
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



