#' @title Probbaility of rank of test given effect size by exact method
#'
#' @description An exact method to comnpute the probbaility of rank of a
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
#' prob <- prob_rank_givenEffect_exact(k=3, et=0, ey=0, nrep=10000, m0=50, m1=50,
#'                                 effectType= "continuous")
#'
#' # compute the probabilities of rank for 1 to 100 tests
#' # do not compute probability for a large number of tests
#' ranks <- 1:100
#' prob <- sapply(ranks, prob_rank_givenEffect_exact, et=1, ey=1, nrep=10000,
#'                                  m0=50, m1=50, effectType = "continuous")
#'
#' # plot
#' plot(ranks, prob)
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by exact method
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
# k0 = ranks under null model
# fun.k0 = input=null rank; output=prob of specific combo of r0 and r1
# t = generate test statistics for target test with effect size et
# p0 = prob of null test having higher test stat value than t
# m = total number of tests
# a = lower limit of the uniform distribution
# b = upper limit of the uniform distribution
# el = vector of uniform effect sizes
# p1 = prob of alt test having higher test stat value than t
# E.T = does importance sampling for the integration over t

# output:-----
# prob = p(rank = k | effect = ey)
# prob is obtained by summing all possible combinations of k0 + k1 = k + 1
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
