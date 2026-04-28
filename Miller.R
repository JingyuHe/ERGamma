tol = 10^(-10)
n = 1000
a = 3 
mu = 5
x = rgamma(n, shape = a, rate = a / mu)



R = sum(log(x))
S = sum(x)
T = S / mu - R + n * log(mu) - n

A = a0 + n / 2
B = b0 + T

a0 = 1
b0 = 1

M = 1000

for(i in 1:M){
    a = A / B
    A = a0 - n * a + n * a^2 * trigamma(a)
    B = b0 + (A - a0) / a - n * log(a) + n * digamma(a) + T

    if( abs(a / (A / B) - 1) < tol ){
        break
    }
}


# Compare approximation
aindex = seq(from = 0.1, to = 10, length.out = 1000)

posterior = rep(0, 1000)

for(i in 1:1000){
    # likelihood * prior
    posterior[i] = sum(log(dgamma(x, shape = aindex[i], rate = aindex[i] / mu))) + log(dgamma(aindex[i], shape = a0, rate = b0))
}
posterior = exp(posterior - max(posterior))
posterior = posterior / sum(posterior)


approx = dgamma(aindex, shape = A, rate = B)
approx = approx / sum(approx)

par(mfrow = c(1,2))
plot(aindex, posterior, type = "l", lwd = 2, main = "PDF")
lines(aindex, approx, col = "red", lwd = 2, lty = 2)
plot(aindex, cumsum(posterior), type = "l", lwd = 2, main = "CDF")
lines(aindex, cumsum(approx), col = "red", lwd = 2, lty = 2)
