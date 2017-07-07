#' @title Weight for the continuous effect sizes
#'
#' @description Compute weight from the probability of the rank given the effect size
#' for the continuous effect size
#' @param alpha Numeric, significance level of the hypothesis test
#' @param et Numeric, mean effect size of the test statistics
#' @param m Integer, totoal number of hypothesis test
#' @param tail Integer (1 or 2), right-tailed or two-tailed hypothesis test.
#' default is right-tailed test.
#' @param delInterval Numeric, interval between the \code{delta} values of a
#' sequence. Note that, \code{delta} is a LaGrange multiplier, necessary to
#' normalize the weight
#' @param ranksProb Numeric vector of ranks probability of the tests given
#' the effect size
#'
#' @details
#' If one wants to test \deqn{H_0: epsilon_i = 0 vs. H_a: \epsilon_i > 0,}
#' then \code{et} and \code{ey} should be mean value of the test and filter
#' effect sizes, respectively. This is called hypothesis testing for the continuous
#' effect sizes.
#'
#' @author Mohamad S. Hasan, shakilmohamad7@gmail.com
#'
#' @export
#'
#' @import OPWeight prob_rank_givenEffect
#'
#' @seealso \code{\link{prob_rank_givenEffect}} \code{\link{weight_binary}}
#'
#' @return \code{weight} Numeric vector of normalized weight of the tests
#' for the continuous case
#'
#' @examples
#'
#' # compute the probabilities of the ranks of a test being rank 1 to 100 if the
#' # targeted test effect is 2 and the overall mean filter effect is 1.
#' ranks <- 1:100
#' prob2 <- sapply(ranks, prob_rank_givenEffect, et = 2, ey = 1, nrep = 10000,
#'                               m0 = 50, m1 = 50)
#'
#' # plot the prooabbility
#' plot(ranks, prob2)
#'
#' # compute weight for the continuous case
#' weight_cont <- weight_continuous(alpha = .05, et = 1, m = 100, tail = 1,
#'                                      delInterval = .0001, ranksProb = prob2)
#'
#' # plot the weight
#' plot(ranks, weight_cont)
#'
#===============================================================================
# function to compute weight from p(rank=k|filterEffect=ey)
#-----------------------------------------------------
# internal parameters:-----
# delta = sequene of delta (lagrange multiplier) values
# findDelta = function to compute sum of weight for each dleta
# deltaOut = optimal delta value
# sumWeight = sum of the weights
# normWeight = normalized weight when necessary
#===============================================================================

weight_continuous <- function(alpha, et, m, tail = 1L, delInterval = .001, ranksProb)
{
    prob <- ranksProb/sum(ranksProb, na.rm = TRUE)
    delta <- seq(0, 1, delInterval)

    weightSumVec <- sapply(delta, weight_by_delta, alpha = alpha, et = et, m = m,
                           m1 = NULL, tail = tail, ranksProb = prob,
                           effectType = "continuous")

    deltaOut <- delta[min(abs(weightSumVec - m)) == abs(weightSumVec - m)]
    deltaOut <- ifelse(length(deltaOut) > 1, deltaOut[1], deltaOut)
    weight.out <- tail*(m/alpha)*pnorm(et/2 + 1/et*log(deltaOut/(alpha*prob)),
                                       lower.tail=FALSE)
    sumWeight <- sum(weight.out, na.rm = TRUE)
    normWeight <- if(sumWeight == 0) {rep(1, m)} else {weight.out/sumWeight*m}
    return(normWeight)
}
