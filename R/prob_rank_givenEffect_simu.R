#' @title Probability of rank of test given effect size by simulations
#'
#' @description A simulation approach to comnpute the probability of rank of a
#' test being higher than any other test given the effect size from the external
#' information.
#' @param s number of samples of test statistics composed of null and alternative
#'  tests
#' @param ey filter test efffect from the external information
#' @param e.one one test effect that will vary across all tests
#' @param m0 number of true null hypothesis
#' @param m1 number of true alternative hypothesis
#' @param effectType type of effect sizes, c("binary","continuous")
#'
#' @details If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i > 0,}
#' then \code{ey} should be mean of the filter effect sizes,
#' This is called hypothesis testing for the continuous effect sizes.\cr
#'
#' If one wants to test \deqn{H_0: epsilon_i=0 vs. H_a: epsilon_i = epsilon,}
#' then \code{ey} should be median or any discrete value of the
#' filter effect sizes. This is called hypothesis testing for the Binary
#' effect sizes.\cr
#'
#' This is a simulation approach to compute the probability of the rank,
#' P(rank | effect = ey) to verify the actual P(rank | effect = ey).
#' Suppose, we have a vector of m = m1+m0 observations, where the first m1
#' observations are
#' from the true alternative and second m0 are from the true null models. If we
#' pick two tests one from the first position and the other from the (m0+1)-th position,
#' then we would expect that the first observation's rank is greater than m0,
#' and (m1+1)-th observation's rank is less than or equal to m1.
#' However, this is not always true, especially when the effect size of the test
#' statistics is low, but the above scenerio become obvious as the the effect
#' size increases. \code{m1} and \code{m0} can be estimated using \code{qvalue} from
#' a bioconductor package \code{qvalue}.
#'
#' @author Mohamad S. Hasan and Paul Schliekelman
#' @export
#' @import stats
#' @seealso \code{\link{runif}} \code{\link{rnorm}} \code{\link{qvalue}}
#' @return \code{r0} rank of the null test statistic\cr
#'         \code{r1} rank of the alternative test statistic
#' @examples
#' # total number of sample generated (use sample size at least 1,000,000)
#' sampleSize = 10000
#' m0 = 50
#' m1 = 50
#' m = m0 +m1
#'
#' # compute rank of the tests
#' rank <- sapply(1:sampleSize, prob_rank_givenEffect_simu, ey = 1, e.one = 1,
#'                           m0 = m0, m1 = m1, effectType = "continuous")
#'
#' # rank may generate missing valuse because of the large effcet size,
#' # therefore, to make a matplot one needs vector of equal size. This procedure
#' # will replace the missing value to make the equal sized vectors
#' # probability of the rank of a null test
#' prob0 <- rep(NA, m)
#' prob0_x <- tapply(rank[1,], rank[1,], length)/sampleSize
#' prob0[as.numeric(names(prob0_x))] <- as.vector(prob0_x)
#'
#' # probability of the rank of an alternative test
#' prob1 <- rep(NA, m)
#' prob1_x <- tapply(rank[2,], rank[2,], length)/sampleSize
#' prob1[as.numeric(names(prob1_x))] <- as.vector(prob1_x)
#'
#' # plot
#' matplot(1:m, cbind(prob0, prob1), type = "l")
#'
#===============================================================================
# function to compute p(rank=k|filterEffect=ey) by simulation
# we used only uniform effects for continuous case.

# Input:-----
# s = number samples of test statistics composed of null and alternative tests
# ey = filter test efffect from the external information
# e.one one test effect that will vary across all tests
# m0 = number of true null hypothesis
# m1 = number of true alternative hypothesis
# effectType = type of effect size c("binary","continuous")

# internal parameters:-----
# m = total number of test
# ey0 = vector of effects of the null tests
# ey1 = vector of effects of the alterantive tests
# t01 = generate test statistics
# r0 = rank of the null test statistic
# r1 = rank of the alternative test statistic

# output:-----
# rank = pair of null and alternative test rank
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


