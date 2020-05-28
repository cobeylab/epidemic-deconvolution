# Example code for deconvolving epidemic curves

Ed Baskerville & Katie Gostic<br>
27 May 2020

This repository contains example code for deconvolving a known delay distribution from observed epidemic time series.

None of the code here is ready for prime time, but we’re sharing it in case it’s a useful starting point for other people.

There is code here for two methods, both with attendant warnings:

* The Richardson-Lucy-type deconvolution used in [Goldstein et al. 2009 PNAS](https://doi.org/10.1073/pnas.0902958106), implemented straightforwardly in unoptimized R code. Either a bug in this code or an inherent property of the algorithm makes it extremely sensitive to initial guess.
* An analogous Bayesian model, implemented straightforwardly in Stan. There are two versions of the model, one (`uncorrelated`) where unobserved states are given identical independent priors, and one (`randomwalk`) where log(unobserved states) follow a Brownian motion. Neither of these implementations have been used for real work, but initial runs seem promising for doing things this way.
