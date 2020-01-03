#' Diagnostics of instrumental variable analysis
#'
#' @param Y A numeric vector of outcomes.
#' @param D A vector of endogenous variables.
#' @param Z A vector of instruments.
#' @param X A vector, matrix or data frame of (exogenous) covariates.
#'
#' @return a list or data frame
#' \describe{
#' \item{x.mean1}{Mean of X under Z = 1 (reported if Z is binary)}
#' \item{x.mean0}{Mean of X under Z = 0 (reported if Z is binary)}
#' \item{coef}{OLS coefficient of X ~ Z (reported if Z is not binary)}
#' \item{se}{Standard error of OLS coefficient (reported if Z is not binary)}
#' \item{p.val}{p-value of the independence of Z and X (Fisher's test if both are binary, logistic regression if Z is binary, linear regression if Z is continuous)}
#' \item{stand.diff}{Standardized difference (reported if Z is binary)}
#' \item{bias.ratio}{Bias ratio}
#' \item{bias.amplify}{Amplification of bias ratio}
#' \item{bias.ols}{Bias of OLS}
#' \item{bias.2sls}{Bias of two stage least squares)}
#' }
#'
#' @author Qingyuan Zhao
#' @references
#' \itemize{
#' \item{Baiocchi, M., Cheng, J., & Small, D. S. (2014). Instrumental variable methods for causal inference. Statistics in Medicine, 33(13), 2297-2340.}
#' \item{Jackson, J. W., & Swanson, S. A. (2015). Toward a clearer portrayal of confounding bias in instrumental variable applications. Epidemiology, 26(4), 498.}
#' \item{Zhao, Q., & Small, D. S. (2018). Graphical diagnosis of confounding bias in instrumental variable analysis. Epidemiology, 29(4), e29--e31.}
#' }
#'
#' @importFrom stats fisher.test
#' @importFrom stats glm
#'
#' @export
#'
#' @examples
#' n <- 10000
#' Z <- rbinom(n, 1, 0.5)
#' X <- data.frame(matrix(c(rnorm(n), rbinom(n * 5, 1, 0.5)), n))
#' D <- rbinom(n, 1, plogis(Z + X[, 1] + X[, 2] + X[, 3]))
#' Y <- D + X[, 1] + X[, 2] + rnorm(n)
#' print(output <- iv.diagnosis(Y, D, Z, X))
#' iv.diagnosis.plot(output)
#'
#' Z <- rnorm(n)
#' D <- rbinom(n, 1, plogis(Z + X[, 1] + X[, 2] + X[, 3]))
#' Y <- D + X[, 1] + X[, 2] + rnorm(n)
#' print(output <- iv.diagnosis(Y, D, Z, X)) ## stand.diff is not reported
#' iv.diagnosis.plot(output)
#'
iv.diagnosis <- function(Y, D, Z, X) {

    if (length(dim(X)) == 2) { ## if x is a matrix/data.frame
        output <- data.frame(t(sapply(1:ncol(X), function(j) unlist(iv.diagnosis(Y, D, Z, X[, j])))))
        rownames(output) <- colnames(X)
        return(output)
    }

    if (length(unique(X)) == 2 & length(unique(Z)) == 2) {
        p.val <- fisher.test(X, Z)$p.value
    } else if (length(unique(X)) == 2) {
        p.val <- summary(glm(X ~ Z, family = "binomial"))$coefficients[2, 4]
    } else {
        p.val <- summary(lm(X ~ Z))$coefficients[2, 4]
    }
    if (length(unique(Z)) == 2) {
        if (!all.equal(sort(unique(Z)), c(0, 1))) {
            stop("Please convert Z to {0, 1}.")
        }
        stand.diff <- (mean(X[Z==1])-mean(X[Z==0]))/(sqrt((var(X[Z==1])+var(X[Z==0]))/2))
    }
    prev.diff.ratio <- (lm(X ~ Z)$coef[2])/(lm(X ~ D)$coef[2]);
    bias.ratio <- prev.diff.ratio / (lm(D ~ Z)$coef[2])
    bias.amplify <- lm(Y ~ X + D)$coef[2]
    bias.ols <- bias.amplify * (lm(X ~ D)$coef[2])
    bias.2sls <- bias.amplify * (lm(X ~ Z)$coef[2]) / (lm(D ~ Z)$coef[2])
    names(bias.ratio) <- NULL
    names(bias.amplify) <- NULL
    names(bias.ols) <- NULL
    names(bias.2sls) <- NULL
    ## stand.diff.high.treatment <- (mean(x[d==1])-mean(x[d==0]))/sqrt((var(x[d==1])+var(x[d==0]))/2);
    if (length(unique(X)) == 2 & length(unique(Z)) == 2) {
        output <- list(x.mean1 = mean(X[Z == 1]),
                       x.mean0 = mean(X[Z == 0]),
                       coef = NA,
                       se = NA,
                       p.val = p.val,
                       stand.diff = stand.diff,
                       bias.ratio = bias.ratio,
                       bias.amplify = bias.amplify,
                       bias.ols = bias.ols,
                       bias.2sls = bias.2sls)
    } else {
        output <- list(x.mean1 = NA,
                       x.mean0 = NA,
                       coef = unname(lm(X ~ Z)$coef[2]),
                       se = summary(lm(X ~ Z))$coefficients[2, 2],
                       p.val = p.val,
                       stand.diff = NA,
                       bias.ratio = bias.ratio,
                       bias.amplify = bias.amplify,
                       bias.ols = bias.ols,
                       bias.2sls = bias.2sls)
    }
    return(output)
}

#' @describeIn iv.diagnosis IV diagnostic plot
#'
#' @param output Output from \code{iv.diagnosis}.
#' @param bias.ratio Add bias ratios (text) to the plot?
#' @param base_size size of the axis labels
#' @param text_size size of the text (bias ratios)
#'
#' @export
#' @import reshape2
#' @import ggplot2
#'
iv.diagnosis.plot<- function(output, bias.ratio = TRUE, base_size = 15, text_size = 5) {

    output <- data.frame(output)
    output$var <- factor(rownames(output), rownames(output)[order(abs(output$bias.ols))])
    rownames(output) <- NULL
    df <- data.frame(output[, c("var", "bias.ols", "bias.2sls")])
    df <- melt(df, id.vars = "var")
    colnames(df) <- c("var", "method", "bias")
    levels(df$method) <- c("ols", "2sls")
    df$method <- factor(df$method, levels = c("2sls", "ols"))

    p <- ggplot() + geom_bar(aes_string(x = "var", y = "bias", fill = "method", linetype = "method"), stat = "identity", color = "black", position=position_dodge(), data = df) + coord_flip() + theme_bw(base_size) + xlab("") + scale_fill_discrete(breaks = c("ols","2sls")) + scale_linetype_discrete(breaks = c("ols","2sls")) + expand_limits(y=c(0, range(df$bias) * 1.2))

    if (bias.ratio) {
        p <- p + geom_text(mapping = aes(y = max(df$bias) * 1.1, x = var, label = round(bias.ratio, 2)), data = data.frame(output), size = text_size)
    }
    p
}
