#' @title Weight for the Binary effect sizes
#'
#' @description Compute weight from the probability of the rank given the effect size
#' for the binary effect size
#' @param alpha significance level of the hypotheis test
#' @param et mean effect size of the test statistics
#' @param m totoal number of hypothesis test
#' @param m1 number of true alternative hypothesis
#' @param tail right-tailed or two-tailed hypothesis test. default is right-tailed test
#' @param delInterval interval between the \code{delta} values of a sequence. Note that,
#' \code{delta} is a LaGrange multiplier, necessary to normalize the weight
#' @param ranksProb probability of the tests given the effect size
#'
#' @details
#' If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i = epsilon,}
#' then \code{et} and \code{ey} should be median or any discrete value of the test
#' and filter effect sizes, respectively. This is called hypothesis testing for
#' the Binary effect sizes. \code{m1} can be estimated using \code{qvalue} from
#' a bioconductor package \code{qvalue}.
#'
#' @author Mohamad S. Hasan, mshasan@uga.edu
#' @export
#' @import OPWeight prob_rank_givenEffect
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_continuous}}
#' \code{\link{qvalue}}
#' @return \code{weight} normalized weight of the tests for the binary case
#' @examples
#'
#' # compute the probabilities of the ranks of a test being rank 1 to 100 if the
#' # targeted test effect is 2 and the overall mean filter effect is 1.
#' ranks <- 1:100
#' prob2 <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 1, nrep = 10000,
#'                               m0 = 50, m1 = 50)
#' # plot the prooabbility
#' plot(ranks, prob2)
#'
#' # compute weight for the binary case
#' weight_bin <- weight_binary(alpha = .05, et = 1, m = 100, m1 = 50, tail=1,
#'                              delInterval = .0001, ranksProb = prob2)
#'
#' # plot the weight
#' plot(ranks, weight_bin)
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
# ranksProb = probability of the tests given the effect size

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
weight_binary <- function(alpha, et, m, m1, tail = 1L, delInterval = .0001, ranksProb)
{
    prob <- ranksProb/sum(ranksProb, na.rm = T)
    delta <- seq(0, 1, delInterval)

    weightSumVec <- sapply(delta, weight_by_delta, alpha = alpha, et = et, m = m,
                           m1 = m1, tail = tail, ranksProb = prob,
                           effectType = "binary")

    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, .0001, deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut*m/(alpha*m1*prob)),
                                       lower.tail = FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}


