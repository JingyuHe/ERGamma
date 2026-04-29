local_lib <- file.path(getwd(), "r-lib")
if (dir.exists(local_lib)) {
  .libPaths(c(local_lib, .libPaths()))
}

library(GIGrvg)
library(gtools)

gamma_const <- -digamma(1)

# Draw n samples from PIG(c) using the JASA parameterization:
# PIG(c) = sum_k GIG(-3/2, 2 c^2, 1/(2 k^2)).
rPIG_ERGamma <- function(n, cc, N) {
  output <- rep(0, n)

  if (N > 0) {
    for (i in 1:N) {
      if (cc == 0) {
        output <- output + 1 / rgamma(n, shape = 3 / 2) / (4 * i^2)
      } else {
        output <- output + rgig(
          n = n,
          lambda = -3 / 2,
          chi = 1 / (2 * i^2),
          psi = 2 * cc^2
        )
      }
    }
  }

  if (cc == 0) {
    temp1 <- trigamma(1 + N)
    temp2 <- psigamma(1 + N, deriv = 2)
    rate <- -4 * temp1 / temp2
    shape <- 0.5 * temp1 * rate
  } else {
    temp1 <- digamma(1 + N + cc) - digamma(1 + N)
    temp2 <- temp1 - cc * trigamma(1 + N + cc)
    rate <- 2 * cc^2 * temp1 / temp2
    shape <- rate / (2 * cc) * temp1
  }

  output + rgamma(n, shape, rate = rate)
}

# Draw one sample from the power-truncated normal kernel
# x^(m - 1) exp(-a x^2 + b x), x > 0.
rH <- function(m, a, b, max_try = 5000) {
  tau <- ifelse(
    b > 0,
    0.5 + sqrt(1 / 4 + 2 * a * m / b^2),
    -0.5 + sqrt(1 / 4 + 2 * a * m / b^2)
  )
  v1 <- tau * abs(b) - b
  v2 <- tau * abs(b) / (2 * a)

  for (i in 1:max_try) {
    x <- rgamma(1, m, rate = v1)
    if (runif(1) <= exp(-a * (x - v2)^2)) {
      return(x)
    }
  }
  NA_real_
}

# Multinomial-Dirichlet Gibbs sampler used by the opioid script.
Gibbs_MD <- function(data, tau = 0, beta = 0, niter = 1000, N_truncations = 300) {
  M <- dim(data)[1]
  K <- dim(data)[2]
  p <- matrix(1, M, K)
  w <- matrix(1, M, K)
  alpha <- matrix(1, niter, K)

  for (i in 2:niter) {
    if (i %% 100 == 0) {
      cat("sampling ", i, "\n")
    }

    for (m in 1:M) {
      p[m, ] <- pmax(rdirichlet(1, alpha[i - 1, ] + data[m, ]), 1e-9)
    }

    eta <- rgamma(M, shape = sum(alpha[i - 1, ]))

    for (k in 1:K) {
      w[, k] <- rPIG_ERGamma(M, alpha[i - 1, k], N_truncations)
    }

    a <- apply(w, 2, sum)
    b <- M * gamma_const - beta + sum(log(eta)) + apply(log(p), 2, sum)
    for (k in 1:K) {
      alpha[i, k] <- rH(M + tau, a[k], b[k])
    }
  }

  alpha
}

# Expectation-only approximation for large Multinomial-Dirichlet runs.
Gibbs_MD_E <- function(data, tau = 0, beta = 0, niter = 1000) {
  M <- dim(data)[1]
  K <- dim(data)[2]
  alpha <- matrix(1, niter, K)

  for (i in 2:niter) {
    if (i %% 100 == 0) {
      cat("sampling ", i, "\n")
    }

    p_e <- t(t(data) + alpha[i - 1, ])
    log_p_e <- digamma(p_e) - digamma(apply(p_e, 1, sum))
    log_eta_e <- digamma(sum(alpha[i - 1, ]))
    w_e <- (gamma_const + digamma(alpha[i - 1, ] + 1)) / (2 * alpha[i - 1, ])

    a <- M * w_e
    b <- M * gamma_const - beta + M * log_eta_e + apply(log_p_e, 2, sum)
    for (k in 1:K) {
      alpha[i, k] <- rH(M + tau, a[k], b[k])
    }
  }

  alpha
}

ltarget_md <- function(data, alpha) {
  res <- sum(apply(data, 1, function(x) {
    lgamma(sum(alpha)) - lgamma(sum(alpha + x)) + sum(lgamma(alpha + x) - lgamma(alpha))
  }))
  res - sum(alpha)
}

MH_MD <- function(data, ltarget = ltarget_md, tau = 1, step, init = NULL, niter = 1000) {
  K <- dim(data)[2]
  x <- matrix(0, niter, K)
  x[1, ] <- if (is.null(init)) abs(rnorm(K, sd = tau)) else init

  for (i in 2:niter) {
    propose <- x[i - 1, ] * exp(rnorm(K, sd = step))
    ratio <- exp(ltarget(data, propose) - ltarget(data, x[i - 1, ]))
    x[i, ] <- if (runif(1) < ratio) propose else x[i - 1, ]
  }

  x
}
