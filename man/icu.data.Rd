\name{icu.data}
\alias{icu.data}
\docType{data}
\title{Pseudo-data based on Branson and Keele (2020)}
\description{
Data sampled with replacemenet from the original data from the (SPOT)light study used in Branson and Keele (2020). Also see Keele et al. (2018) for more details about the variables in this dataset.}
\usage{data(icu.data)}
\format{
  A data frame with 13011 observations on the following 18 variables.
  \describe{
    \item{\code{age}}{Age of the patient in years.}
    \item{\code{male}}{Whether or not the patient is male; 1 if male and 0 otherwise.}
    \item{\code{sepsis_dx}}{Whether or not the patient is diagnosed with sepsis; 1 if so and 0 otherwise.}
    \item{\code{periarrest}}{Whether or not the patient is diagnosed with peri-arrest; 1 if so and 0 otherwise.}
    \item{\code{icnarc_score}}{The Intensive Care National Audit and Research Centre physiological score.}
    \item{\code{news_score}}{The National Health Service national early warning score.}
    \item{\code{sofa_score}}{The sequential organ failure assessment score.}
    \item{\code{v_cc1}}{Indicator for level of care at assessment (Level 0, normal ward care).}
    \item{\code{v_cc2}}{Indicator for level of care at assessment (Level 1, normal ward care).}
    \item{\code{v_cc4}}{Indicator for level of care at assessment (Level 2, care within a high dependency unit).}
    \item{\code{v_cc5}}{Indicator for level of care at assessment (Level 3, ICU care).}
    \item{\code{v_cc_r1}}{Indicator for recommended level of care at assessment (Level 0, normal ward care).}
    \item{\code{v_cc_r2}}{Indicator for recommended level of care after assessment (Level 1, normal ward care).}
    \item{\code{v_cc_r4}}{Indicator for recommended level of care after assessment (Level 2, care within a high dependency unit).}
    \item{\code{v_cc_r5}}{Indicator for recommended level of care after assessment (Level 3, ICU care).}
    \item{\code{open_bin}}{Binary instrument; 1 if the available number of ICU beds was less than 4, and 0 otherwise.}
    \item{\code{icu_bed}}{Binary treatment; 1 if admitted to an ICU bed.}
    \item{\code{site}}{ID for the hospital that the patient attended.}
  }
}
\references{
Keele, L. et al. (2018). Stronger instruments and refined covariate balance in an observational study of the effectiveness of prompt admission to intensive care units. Journal of the Royal Statistical Society: Series A (Statistics in Society).

Branson, Z. and Keele, L. (2020). Evaluating a Key Instrumental Variable Assumption Using Randomization Tests. American Journal of Epidemiology. To appear.
}
\examples{
data(icu.data)
}
\keyword{datasets}
