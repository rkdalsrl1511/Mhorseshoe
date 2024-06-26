# Mhorseshoe

## Overview

Mhorseshoe is a package for a high-dimensional Bayesian linear modeling 
algorithm using a horseshoe prior. A feature of this package is that it 
implements approximate MCMC algorithm from Johndrow et al. (2020) and provides 
a horseshoe estimator that can effectively reduce computational costs for 
high-dimensional sparse data. This package provides two different algorithm 
functions :

-`exact_horseshoe()` Run the horseshoe estimator.

-`approx_horseshoe()` Run the horseshoe estimator with the approximate algorithm applied.

## Installation

```r
install.package("Mhorseshoe")
```

## Usage

The following linear model assumptions are made.

$$L(y\ |\ x, \beta, \sigma^2) = (\frac{1}{\sqrt{2\pi}\sigma})^{-N/2}exp
\{ -\frac{1}{2\sigma^2}(y-X\beta)^T(y-X\beta)\},\\ X \in 
\mathbb{R}^{N \times p},\ y \in \mathbb{R}^{N},\ \beta \in \mathbb{R}^{p}$$

- $X \in \mathbb{R}^{N \times p}$ : Matrix of covariates.
- $y \in \mathbb{R}^{N}$ : Response variable.

```r
# Run functions from the Mhorseshoe package
ex_result <- exact_horseshoe(y, X, burn = 5000, iter = 10000)
ap_result <- approx_horseshoe(y, X, burn = 5000, iter = 10000)

# posterior mean of beta
ex_betahat <- ex_result$BetaHat
ap_betahat <- ap_result$BetaHat

# 95% posterior credible intervals
ex_LeftCI <- ex_result$LeftCI
ex_RightCI <- ex_result$RightCI
ap_LeftCI <- ap_result$LeftCI
ap_RightCI <- ap_result$RightCI
```

## References

Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020). Scalable Approximate 
MCMC Algorithms for the Horseshoe Prior. In Journal of Machine Learning 
Research (Vol. 21).

If you would like to discuss this package, please email leehuimin115@g.skku.edu
