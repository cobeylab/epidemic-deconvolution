# Examples of deconvolving epidemic curves

Ed Baskerville & Katie Gostic<br>
27 May 2020

This repository contains examples of deconvolving epidemic curves using two methods:

* (soon) The Richardson-Lucy-type deconvolution used in [Goldstein et al. 2009 PNAS](https://doi.org/10.1073/pnas.0902958106), implemented straightforwardly in unoptimized R code.
* An analogous Bayesian model, implemented straightforwardly in Stan. There are two versions of the model, one where unobserved states are given identical independent priors, and one where log(unobserved states) follow a random walk.
