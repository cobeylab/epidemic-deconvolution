# This would be more efficient if it just used convolve() directly
convolve_delay <- function(x, p) {
  nx <- length(x)
  np <- length(p)
  y <- rep(0, nx + np)
  for(i in 1:nx) {
    y[(i + 1):(i + np)] <- y[(i + 1):(i + np)] + x[i] * p
  }
  y
}
