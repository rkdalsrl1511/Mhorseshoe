# Mhorseshoe

## Overview

Mhorseshoe is a package for a high-dimensional Bayesian linear modeling 
algorithm using a horseshoe prior. A feature of this package is that it 
implements approximate MCMC algorithm from Johndrow et al. (2020) and provides 
a horseshoe estimator that can effectively reduce computational costs for 
high-dimensional sparse data. This package provides three different algorithm 
functions :

-`exact_horseshoe()` Run the horseshoe estimator assuming a linear model.
-`approx_horseshoe()` Run the horseshoe estimator with the approximate algorithm applied.
-`mapprox_horseshoe()` Run the horseshoe estimator, which estimates and updates the threshold of the approximate algorithm.

mapprox_horseshoe updates parameters in the same way as approx_horseshoe, but while 
approx_horseshoe uses a fixed value as the threshold of the approximation algorithm, 
mapprox_horseshoe adds a new threshold update process. More information about these 
can be found in `vignette("Mhorseshoe")`.

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
# Run functions from the Mhorseshoe package with default settings
result <- exact_horseshoe(X, y, iteration = 5000)
result <- approx_horseshoe(X, y, iteration = 5000)
result <- mapprox_horseshoe(X, y, iteration = 5000)

# posterior mean(burn-in = 1000)
post_mean <- apply(result$BetaSamples[1001:5000, ], MARGIN = 2, mean)

# 95% posterior credible intervals(burn-in = 1000)
post_leftCI <- apply(result$BetaSamples[1001:5000, ], MARGIN = 2, quantile, probs = 0.025)
post_rightCI <- apply(result$BetaSamples[1001:5000, ], MARGIN = 2, quantile, probs = 0.975)
```

## References

Johndrow, J., Orenstein, P., & Bhattacharya, A. (2020). Scalable Approximate 
MCMC Algorithms for the Horseshoe Prior. In Journal of Machine Learning 
Research (Vol. 21).

If you would like to discuss this package, please email leehuimin115@g.skku.edu
