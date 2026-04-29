library(EBarrays)

## EM algorithm 
## Lognormal-Normal Demo

## mu10,sigma2,tau are parameters in the LNNB model; pde is the
## proportion of differentially expressed genes; n is the
## total number of genes; nr1 and nr2 are the number of replicate
## arrays in each group.

lnnb.sim <- function(mu10, sigmasq, tausq, pde, n, nr1, nr2)
{
    de <- sample(c(TRUE, FALSE), size = n, replace = TRUE, prob = c(pde, 1 - pde))
    x <- matrix(NA, n, nr1)
    y <- matrix(NA, n, nr2)
    mu1 <- rnorm(n, mu10, sqrt(tausq))
    mu2.de <- rnorm(n, mu10, sqrt(tausq))
    mu2 <- mu1
    mu2[de] <- mu2.de[de]
    for(j in 1:nr1) {
        x[, j] <- rnorm(n, mu1, sqrt(sigmasq))
    }
    for(j in 1:nr2) {
        y[, j] <- rnorm(n, mu2, sqrt(sigmasq))
    }
    outmat <- exp(cbind(x, y))
    list(mu1 = mu1, mu2 = mu2, outmat = outmat, de = de)
}

## simulating data with
##  mu_0 = 2.33, sigma^2 = 0.1, tau^2 = 2
##  P(DE) = 0.2

sim.data1 <- lnnb.sim(2.33, 0.1, 2, 0.2, 2000, nr1 = 3, nr2 = 3)
de.true1 <- sim.data1$de ## true indicators of differential expression

sim.data2 <- lnnb.sim(1.33, 0.01, 2, 0.2, 2000, nr1 = 3, nr2 = 3)
de.true2 <- sim.data2$de ## true indicators of differential expression

testdata <- rbind(sim.data1$outmat,sim.data2$outmat)

hypotheses <- ebPatterns(c("1 1 1 1 1 1", "1 1 1 2 2 2")) 

em.out <- emfit(testdata, family = "LNN", hypotheses,
                cluster = 1:5,
                type = 2,
                verbose = TRUE,
                num.iter = 10)

em.out

post.out <- postprob(em.out, testdata)

table(post.out$pattern[, 2] > .5, c(de.true1,de.true2))
table((post.out$cluster[, 2] > .5)+1, c(rep("Cluster 1",2000),rep("Cluster 2",2000)))

plotMarginal(em.out,testdata)
par(ask=TRUE)
plotCluster(em.out,testdata)
par(ask=FALSE)

lnnmv.em.out <- emfit(testdata, family = "LNNMV", hypotheses, groupid=c(1,1,1,2,2,2),
                verbose = TRUE,
                num.iter = 10,
                p.init = c(0.95, 0.05))

lnnmv.em.out
post.out <- postprob(lnnmv.em.out, testdata, groupid=c(1,1,1,2,2,2))

table(post.out$pattern[, 2] > .5, c(de.true1,de.true2))
