
source(file.path("src","model2","v1","prob_dens.R"))

alpha=1
beta  =10
n <- 1000

x    <- rinvgamma(n, shape= alpha, rate=beta)

est <- ml_inversegamma(x)

cat('alpha = ',est[1],' (true=',alpha,")\n")
cat('beta = ',est[2],' (true=',beta,")\n")