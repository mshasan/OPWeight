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
#' prob <- prob_rank_givenEffect(k=3,et=0,ey=0,nrep=10000,m0=50,m1=50)
#'
#' # compute the probabilities of rank for 1 to 100 tests
#' ranks <- 1:100
#' prob <- sapply(ranks,prob_rank_givenEffect,et=1,ey=1,nrep=10000,m0=50,m1=50)
#'
#' # plot
#' plot(ranks,prob)
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by normal approximation

# Input:-----
# k = rank of the tests
# et = actual data test effect for importance sampling
# ey = filter test efffect from external information
# nrep = number of replications for importance sampling
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis

# internal parameters:-----
# t = generate test statistics for target test with effect size et
# p0 = prob of null test having higher test stat value than t
# p1 = prob of alt test having higher test stat value than t

# output:-----
# prob = p(rank=k|effect=ey)
#===============================================================================
prob_rank_givenEffect <- function(k, et, ey, nrep = 10000, m0, m1)
{
    t <- rnorm(nrep, et, 1)
    p0 <- pnorm(-t)
    p1 <- pnorm(ey - t)

    mean0 <- (m0 - 1)*p0 + m1*p1 + 1
    mean1 <- m0*p0 + (m1 - 1)*p1 + 1

    var0 <- (m0 - 1)*p0*(1 - p0) + m1*p1*(1 - p1)
    var1 <- m0*p0*(1 - p0) + (m1 - 1)*p1*(1 - p1)

    prob <- ifelse(et == 0, mean(dnorm(k, mean0, sqrt(var0))),
                   mean(dnorm(k, mean1, sqrt(var1))))
    return(prob)
}








