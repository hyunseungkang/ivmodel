# ivmodel

ivmodel is a comprehensive R package for linear instrumental variables analysis with one endogenous variable (e.g. confounded treatment/policy/intervention). The package can handle large data, data from matched sets, and is designed to be robust to many data types (or mixtures of data types). The package contains the usual instrumental variables methods, methods for weak instruments, sensitivity analysis, and power calculations. Please see our paper Jiang, Kang, and Small (2015) for details. 

To install this package in R from GitHub, run the following commands:

```R
install.packages("devtools")
library(devtools) 
install_github("hyunseungkang/ivmodel")
```

To install this package in R from CRAN, run the following command:

```R
install.packages("ivmodel")
```

## Examples

```R
library(ivmodel)

### Use Card (1995) data to estimate the effect of education on log earnings ###
# n: sample size; L: # of IVs; p: # of covariates 
# Y: n by 1 vector of outcomes (must be continuous)
# D: n by 1 vector of treatments (continuous or discrete)
# Z: n by L vector of instruments (continuous or discrete)
# X: n by L vector of instruments (continuous or discrete)
data(card.data)

# One Instrument Anaylsis with Proximity to 4yr College as IV#
Y=card.data[,"lwage"]
D=card.data[,"educ"] 
Z=card.data[,"nearc4"]
Xname=c("exper", "expersq", "black", "south", "smsa", "reg661",
"reg662", "reg663", "reg664", "reg665", "reg666", "reg667",
"reg668", "smsa66")
X=card.data[,Xname]

# Run command #
card.model1IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model1IV

# Multiple IV Analysis with Proximito to 2yr and 4yr Colleges as IVs#
Z = card.data[,c("nearc4","nearc2")]
card.model2IV = ivmodel(Y=Y,D=D,Z=Z,X=X)
card.model2IV
```

## References 
Card, D. (1995). "Using Geographic Variations in College Proximity to Estimate the Return
to Schooling.” In LN Christofides, EK Grant, R Swidinsky (eds.), Aspects of Labor Market
Behaviour: Essays in Honour of John Vanderkamp. University of Toronto Press. 

Kang, H., Jiang, Y., Zhao, Q., and Small, D. S. (2015). <a href="http://pages.stat.wisc.edu/~hyunseung/">ivmodel: An R Package for Inference and Sensitivity Analysis of Instrumental Variables Models with One Endogenous Variable.</a> Technical Report.

