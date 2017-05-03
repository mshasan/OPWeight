#' @title Probbaility of rank of test given effect size by normal approximation
#'
#' @description A normal approximation to comnpute the probbaility of rank of a
#' test being higher than any other test given the effect size from external
#' information.
#' @param k rank of the tests
#' @param et actual data test effect for importance sampling
#' @param ey filter test efffect from external information
#' @param nrep number of replications for importance sampling
#' @param m0 number of true null hypothesis
#' @param m1 number of true alternative hypothesis
#' @param effectType type of effect size c("binary","continuous")
#'
#' @details If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i > 0,}
#' then \code{et}  and \code{ey} should be mean of the test and filter effect sizes,
#' respectively. This is called hypothesis testing for the continuous effect sizes.
#' If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i = epsilon,}
#' then \code{et} and \code{ey} should be median or any discrete value of the
#' test and filter effect sizes. This is called hypothesis testing for the Binary
#' effect sizes.
#' @author Mohamad S. Hasan, mshasan@uga.edu
#' @export
#' @import stats
#' @seealso \code{\link{dnorm}} \code{\link{pnorm}} \code{\link{rnorm}}
#' @return \code{prob} probability of the rank of the test
#' @examples
#' # compute the probability of rank for the third test
#' prob <- prob_rank_givenEffect_approx(k=3, et=0, ey=0, nrep=10000, m0=50, m1=50,
#'                          effectType = "continuous")
#'
#' # compute the probabilities of rank for 1 to 100 tests
#' ranks <- 1:100
#' prob <- sapply(ranks, prob_rank_givenEffect_approx, et=2, ey=2, nrep=10000,
#'                      m0=50, m1=50, effectType = "binary")
#'
#' # plot
#' plot(ranks, prob)
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by normal approximation
# we used only uniform effects for continuous case.

# Input:-----
# k = rank of the tests
# et = actual data test effect for importance sampling
# ey = filter test efffect from external information
# nrep = number of replications for importance sampling
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
# effectType = type of effect size c("binary","continuous")

# internal parameters:-----
# t = generate test statistics for target test with effect size et
# p0 = prob of null test having higher test stat value than t
# m = total number of tests
# a = lower limit of the uniform distribution
# b = upper limit of the uniform distribution
# el = vector of uniform effect sizes
# p1 = prob of alt test having higher test stat value than t

# output:-----
# prob = p(rank=k|effect=ey)
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




