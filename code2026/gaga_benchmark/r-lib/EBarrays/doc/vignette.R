### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: setup
###################################################
options(width=50)


###################################################
### code chunk number 2: init
###################################################
library(EBarrays)


###################################################
### code chunk number 3: readpattern
###################################################
pattern <- ebPatterns(c("1, 1, 1, 1, 1, 1, 1, 1, 1, 1",
                        "1, 1, 1, 1, 1, 2, 2, 2, 2, 2"))
pattern


###################################################
### code chunk number 4: readdata
###################################################
data(gould)


###################################################
### code chunk number 5: read2patterns
###################################################
pattern <- ebPatterns(c("1,1,1,0,0,0,0,0,0,0","1,2,2,0,0,0,0,0,0,0"))

pattern


###################################################
### code chunk number 6: emfit
###################################################

gg.em.out <- emfit(gould, family = "GG", hypotheses = pattern, 
                   num.iter = 10)
gg.em.out

gg.post.out <- postprob(gg.em.out, gould)$pattern

gg.threshold <- crit.fun(gg.post.out[,1], 0.05)

gg.threshold

sum(gg.post.out[,2] > gg.threshold)

lnn.em.out <- emfit(gould, family = "LNN", hypotheses = pattern, 
                    num.iter = 10)
lnn.em.out

lnn.post.out <- postprob(lnn.em.out, gould)$pattern

lnn.threshold <- crit.fun(lnn.post.out[,1], 0.05)

lnn.threshold

sum(lnn.post.out[,2] > lnn.threshold)

lnnmv.em.out <- emfit(gould, family = "LNNMV", hypotheses = pattern, 
                      groupid = c(1,2,2,0,0,0,0,0,0,0), num.iter = 10)
lnnmv.em.out

lnnmv.post.out <- postprob(lnnmv.em.out, gould, 
                           groupid = c(1,2,2,0,0,0,0,0,0,0))$pattern

lnnmv.threshold <- crit.fun(lnnmv.post.out[,1], 0.05)

lnnmv.threshold

sum(lnnmv.post.out[,2] > lnnmv.threshold)

sum(gg.post.out[,2] > gg.threshold & lnn.post.out[,2] > lnn.threshold)

sum(gg.post.out[,2] > gg.threshold & lnnmv.post.out[,2] > lnnmv.threshold)

sum(lnn.post.out[,2] > lnn.threshold & lnnmv.post.out[,2] > lnnmv.threshold)

sum(gg.post.out[,2] > gg.threshold & lnn.post.out[,2] > lnn.threshold
    & lnnmv.post.out[,2] > lnnmv.threshold)



###################################################
### code chunk number 7: vignette.Rnw:577-582
###################################################
pattern4 <- ebPatterns(c("1, 1, 1, 1, 1, 1, 1, 1, 1, 1", 
                         "1, 2, 2, 2, 2, 2, 2, 2, 2, 2", 
                         "1,2,2,1,1,1,1,1,2,2", 
                         "1,1,1,1,1,1,1,1,2,2"))
pattern4


###################################################
### code chunk number 8: fourgroups
###################################################
gg4.em.out <- emfit(gould, family = "GG", pattern4, 
                    num.iter = 10)

gg4.em.out

gg4.post.out <- postprob(gg4.em.out, gould)$pattern

gg4.threshold2 <- crit.fun(1-gg4.post.out[,2], 0.05)

sum(gg4.post.out[,2] > gg4.threshold2)

lnn4.em.out <- emfit(gould, family="LNN", pattern4,
                     num.iter = 10)

lnn4.em.out

lnn4.post.out <- postprob(lnn4.em.out, gould)$pattern

lnn4.threshold2 <- crit.fun(1-lnn4.post.out[,2], 0.05)

sum(lnn4.post.out[,2] > lnn4.threshold2)

lnnmv4.em.out <- emfit(gould, family="LNNMV", pattern4,
                       groupid = c(1,2,2,3,3,3,3,3,4,4), num.iter = 10)

lnnmv4.em.out

lnnmv4.post.out <- postprob(lnnmv4.em.out, gould,
                            groupid = c(1,2,2,3,3,3,3,3,4,4))$pattern

lnnmv4.threshold2 <- crit.fun(1-lnnmv4.post.out[,2], 0.05)

sum(lnnmv4.post.out[,2] > lnnmv4.threshold2)

sum(gg4.post.out[,2] > gg4.threshold2 & lnn4.post.out[,2] > lnn4.threshold2)

sum(gg4.post.out[,2] > gg4.threshold2 & lnnmv4.post.out[,2] > lnnmv4.threshold2)

sum(lnn4.post.out[,2] >  lnn4.threshold2 & lnnmv4.post.out[,2] > lnnmv4.threshold2)

sum(gg4.post.out[,2] > gg4.threshold2 & lnn4.post.out[,2] > lnn4.threshold2 
    & lnnmv4.post.out[,2] > lnnmv4.threshold2) 


###################################################
### code chunk number 9: checkccv
###################################################
checkCCV(gould[,1:3])


###################################################
### code chunk number 10: qqplotgamma
###################################################
print(checkModel(gould, gg.em.out, model = "gamma", nb = 50))


###################################################
### code chunk number 11: qqplotlnorm
###################################################
print(checkModel(gould, lnn.em.out, model = "lognormal", nb = 50))


###################################################
### code chunk number 12: qqplotlnnmv
###################################################
print(checkModel(gould, lnnmv.em.out, model = "lnnmv", nb = 50,
		 groupid = c(1,2,2,0,0,0,0,0,0,0)))


###################################################
### code chunk number 13: ggmarg
###################################################
print(plotMarginal(gg.em.out, gould))


###################################################
### code chunk number 14: lnnmarg
###################################################
print(plotMarginal(lnn.em.out, gould))


###################################################
### code chunk number 15: varsqqplotlnnmv
###################################################
print(checkVarsQQ(gould, groupid = c(1,2,2,0,0,0,0,0,0,0)))


###################################################
### code chunk number 16: varsmarglnnmv
###################################################
print(checkVarsMar(gould, groupid = c(1,2,2,0,0,0,0,0,0,0), 
                   nint=2000, xlim=c(-0.1, 0.5)))


###################################################
### code chunk number 17: gg4marg
###################################################
print(plotMarginal(gg4.em.out, data = gould))


###################################################
### code chunk number 18: lnn4marg
###################################################
print(plotMarginal(lnn4.em.out, data = gould))


###################################################
### code chunk number 19: varsmarglnnmv4
###################################################
print(checkVarsMar(gould, groupid = c(1,2,2,3,3,3,3,3,4,4), 
      	           nint=2000, xlim=c(-0.1, 0.5)))


