#' @title Weight for the Binary effect sizes
#'
#' @description Compute weight from the proability of rank given the effect size
#' for the binary effect size
#' @param alpha significance level of the hypotheis test
#' @param et effect size from the actual data (test effect)
#' @param m totoal number of hypothesis test
#' @param m1 number of true alternative hypothesis
#' @param tail one-tailed or two-tailed hypothesis test. default is two-tailed test
#' @param delInterval interval between the delta values of a sequence
#' @param prob probability of the rank given the effect size
#'
#' @details
#' If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i = epsilon,}
#' then \code{et} and \code{ey} should be median or any discrete value of the test
#' and filter effect sizes. This is called hypothesis testing for the Binary effect
#' sizes. m1 can be estimated using \code{qvalue} from \code{qvalue}
#' in bioconductor or \code{propTrueNull} from \code{limma} package
#' @author Mohamad S. Hasan, mshasan@uga.edu
#' @export
#' @import OPWeight prob_rank_givenEffect
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_continuous}}
#' @return \code{weight} normalized weight of the tests for the binary case
#' @examples
#' # compute the probabilities of rank for 1 to 100 tests
#' ranks <- 1:100
#' prob <- sapply(ranks,prob_rank_givenEffect,et=1,ey=1,nrep=10000,m0=50,m1=50)
#'
#' # plot the prooabbility
#' plot(ranks,prob)
#'
#' # compute weight for the binary case
#' weight_bin <- weight_binary(alpha=.05,et=1,m=100,m1=50,tail=2,delInterval=.0001,
#' prob=prob)
#'
#' # plot the weight
#' plot(ranks,weight_bin)
#'
#===============================================================================
# function to compute weight from p(rank=k|filterEffect=ey)

# Input:-----
# alpha = significance level of the hypotheis test
# et = effect size from the actual data (test effect)
# m = totoal number of hypothesis test
# m1 = number of true alternative hypothesis
# tail = one-tailed or two-tailed hypothesis test
# delInterval =  interval between the delta values of a sequesnce
# prob =probability of rank given the effect size

# internal parameters:-----
# delta = sequene of delta (lagrange multiplier) values
# findDelta = function to compute sum of weight for each dleta
# deltaOut = optimal delta value
# sumWeight = sum of the weights
# normWeight = normalized weight when necessary


# output:-----
# # Weight.out = weight without normalization
#===============================================================================

# function to compute weight binary case
#--------------------------------------
weight_binary <- function(alpha, et, m, m1, tail = 2L, delInterval = .0001, prob)
{
    prob <- prob/sum(prob, na.rm = T)
    delta <- seq(0, 1, delInterval)
    findDelta <- function(delta)
    {
        weight <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(delta*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
        return(sum(weight, na.rm = TRUE))
    }
    weightSumVec <- vapply(delta, findDelta, 1)
    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}


