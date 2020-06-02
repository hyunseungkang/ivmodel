\name{ivmodel-bransonKeele-internal}
\title{Internal \code{ivmodel} functions for Branson and Keele (2020)}
\alias{permuteData.biasedCoin}
\alias{getBlockPerm}
\alias{getCompletePerms.meanDiffs}
\alias{getCompletePerms.balance}
\alias{getCompletePerms.md}
\alias{getCompletePerms.absBias}
\alias{getBlockPerms.md}
\alias{getBlockPerms.absBias}
\alias{getBernoulliPerms.md}
\alias{getBernoulliPerms.absBias}

\description{Internal \code{ivmodel} functions for Branson and Keele (2020)}
\usage{
	permuteData.biasedCoin(N, probs)
	getBlockPerm(subclassIndicatorTable)
	getCompletePerms.meanDiffs(X, indicator, perms = 1000)
	getCompletePerms.balance(X, indicator, perms = 1000)
	getCompletePerms.md(X, indicator, perms = 1000)
	getCompletePerms.absBias(X, D = NULL, Z = NULL, perms = 1000)
	getBlockPerms.md(X, indicator, subclass, perms = 1000)
	getBlockPerms.absBias(X, D = NULL, Z = NULL, subclass = NULL, perms = 1000)
	getBernoulliPerms.md(X, indicator, perms = 1000)
	getBernoulliPerms.absBias(X, D = NULL, Z = NULL, perms = 1000)
}
\author{Zach Branson and Luke Keele}
\details{
	\code{permuteData.biasedCoin} permutes the treatment indicator according to biased-coin randomization (i.e., Bernoulli trials).
	\code{getBlockPerm} permutes an indicator (instrument or exposure) within a subclass.
	\code{getCompletePerms.meanDiffs} returns the covariate mean differences across many permutations of an indicator.
	\code{getCompletePerms.balance} returns the standardized covariate mean differences across many permutations of an indicator.
	\code{getCompletePerms.md} returns the Mahalanobis distance across many permutations of an indicator.
	\code{getCompletePerms.absBias} returns the sum of absolute biases across many permutations of an indicator.
	\code{getBlockPerms.md} returns the Mahalanobis distance across many block permutations of an indicator.
	\code{getBlockPerms.absBias} returns the sum of absolute biases across many block permutations of an indicator.
	\code{getBernoulliPerms.md} returns the Mahalanobis distance across many Bernoulli-trial permutations of an indicator.
	\code{getBernoulliPerms.absBias} returns the sum of absolute biases across many Bernoulli-trial permutations of an indicator.
}
\keyword{internal}