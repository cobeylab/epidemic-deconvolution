#' Convolve a signal with a delay distribution.
#' 
#' This function takes a signal and a delay distribution and produces
#' a signal with the delay distribution applied.
#' 
#' With `x[t]` the epidemic incidence at time `t`, and `pmf[i + 1]` the
#' fraction of incident cases that are observed at time `t + i`, then
#' the returned `y[s]` will be the total number of cases observed
#' at time `s`.
#' 
#' (Of course, due to the algebraic properties of convolution, it could
#' be used/interpreted in other ways, but was written for this purpose.)
#' 
#' Main limitation: this function assumes nonnegative delays.
#' For real-world situations, you probably will need something more
#' general-purpose.
#' 
#' This is a naive, direct implementation of the convolution sum,
#' and would be more efficient on large data if it just used
#' convolve() directly, since convolve() uses the FFT.
#' (I didn't feel like getting a handle on the indexing conventions of
#' `convolve()` when needing this.)
#' 
#' @param x A signal, treated as a numeric vector.
#' @param p The delay distribution, where `pmf[i + 1]` is the fraction of
#' the signal delayed by `i` timesteps. (That is, `pmf[1]` contains the
#' undelayed fraction.)
#' 
#' @return The delayed signal `y`, a numeric vector with length
#' `length(x) + length(pmf) - 1`.
#' 
#' @author Ed Baskerville
convolve_delay <- function(x, pmf) {
  stopifnot(is.numeric(x))
  stopifnot(is.numeric(pmf))
  stopifnot(abs(sum(pmf) - 1) < 1e-10)
  
  nx <- length(x)
  np <- length(pmf)
  y <- rep(0, nx + np - 1)
  for(i in 1:nx) {
    y[i:(i + np - 1)] <- y[i:(i + np - 1)] + x[i] * pmf
  }
  y
}
