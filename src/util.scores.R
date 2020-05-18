#####################
#CRPS functions
#####################
##
# evalutes the CRPS for a vector y under normal assumption
# y     - (n x 1) observations
# mu    - (n x 1) predicted mean for each observations
# sigma - (n x 1) predicted standard deviation for each observations
##
CRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) + 0.5 * Exx(mu, sigma))
}

CRPS.sharpness <- function(mu, sigma)
{
  return(-0.5 * Exx(mu, sigma))
}

E_CRPS <- function(mu, sigma, mu_hat, sigma_hat)
{
  return(-EExy(mu_hat, sigma_hat, mu, sigma) + 0.5 * Exx(mu_hat, sigma_hat))
}

#####################
#SCRPS functions

#####################


#####################
##
# evalutes the scaled CRPS for a vector y under normal assumption
# y     - (n x 1) observations
# mu    - (n x 1) predicted mean for each observations
# sigma - (n x 1) predicted standard deviation for each observations
##
SCRPS <- function(y, mu, sigma)
{
  return(-Exy(mu, sigma, y) / Exx(mu, sigma) - 0.5 * log(Exx(mu, sigma)))
}

SCRPS.sharpness <- function(mu, sigma)
{
  return(-1 - 0.5 * log(Exx(mu, sigma)))
}

E_SCRPS <- function(mu, sigma, mu_hat, sigma_hat)
{
  return(-EExy(mu, sigma, mu_hat, sigma_hat) / Exx(mu_hat, sigma_hat) - 0.5 *
           log(Exx(mu_hat, sigma_hat)))
}

#####################
#log-score functions
#####################
##
# evalutes the log score for a vector y under normal assumption
# y     - (n x 1) observations
# mu    - (n x 1) predicted mean for each observations
# sigma - (n x 1) predicted standard deviation for each observations
##
LS <- function(y, mu, sigma)
{
  return(dnorm(y, mean = mu, sd = sigma, log = TRUE))
}

LS.sharpness <- function(mu, sigma)
{
  return(-0.5 * log(2 * pi) - 0.5 - log(sigma))
}

E_LS <- function(mu, sigma, mu_hat, sigma_hat)
{
  return(-(sigma ^ 2 + (mu - mu_hat) ^ 2) / (2 * sigma_hat ^ 2) - log(sqrt(2 *
                                                                             pi) * sigma_hat))
}

#####################
#rCRPS functions
#####################
rCRPS <- function(y, c, mu, sigma)
{
  return(0.5 * Exxc(c, 0, sqrt(2) * sigma) - Exyc(y, c, mu, sigma))
}

rCRPS.sharpness <- function(c, mu, sigma)
{
  return(-0.5 * Exxc(c, 0, sqrt(2) * sigma))
}

E_rCRPS <- function(c, mu, sigma, mu_hat, sigma_hat)
{
  return(0.5 * Exxc(c, 0, sqrt(2) * sigma_hat) - EExyc(c, mu_hat, sigma_hat, mu, sigma))
}


#####################
#rSCRPS functions
#####################
rSCRPS <- function(y, c, mu, sigma)
{
  return(-Exyc(y, c, mu, sigma) / Exxc(c, 0, sqrt(2) * sigma) - 0.5 * log(Exxc(c, 0, sqrt(2) *
                                                                                 sigma)))
}

rSCRPS.sharpness <- function(c, mu, sigma)
{
  return(-1 - 0.5 * log(Exxc(c, 0, sqrt(2) * sigma)))
}

E_rSCRPS <- function(c, mu, sigma, mu_hat, sigma_hat)
{
  return(-EExyc(c, mu_hat, sigma_hat, mu, sigma) / Exxc(c, 0, sqrt(2) * sigma_hat) - 0.5 *
           log(Exxc(c, 0, sqrt(2) * sigma_hat)))
}


#####################
#utility functions
#####################

Exyc <- function(y, c, mu, sigma) {
  #compute E(c-g(|y-X|)) when X is N(mu,sigma) and g(x) = c - |x| for x<c.
  return(c - EcXnorm(c, y - mu, sigma))
}

Exxc <- function(c, mu, sigma) {
  #compute E(c-g(|X|)) when X is N(mu,sigma) and g(x) = c - |x| for x<c.
  return(c - EcXnorm(c, mu, sigma))
}

EExyc <- function(c, mu, sigma, mu.y, sigma.y) {
  #compute E(E(c-g(|y-X|))) when X is N(mu,sigma) Y is N(mu.y,sigma.y)
  #and g(x) = c - |x| for x<c.
  return(c - EEcXnorm(c, mu, sigma, mu.y, sigma.y))
}

EcXnorm <- function(c, mu, sigma) {
  #Compute E((c-|X|)*1(|X|<c)) when X is N(mu,sigma)
  #E((c-|X|)*1(|X|<c)) = E(1(|X|<c)) - E(|X|*1(|x|<c))
  Int1 = PcXnorm(-c, c, mu, sigma)
  Int2 = -EcIntnorm(-c, 0, mu, sigma) + EcIntnorm(0, c, mu, sigma)
  return(c * Int1 - Int2)
}

EEcXnorm <- function(c, mu, sigma, mu.y, sigma.y) {
  #Compute E(E((c-|Y-X|)*1(|Y-X|<c))) when X is N(y-mu,sigma)
  #and Y is N(mu.y,sigma.y)
  mu.xy = mu.y - mu
  sigma.xy = sqrt(sigma ^ 2 + sigma.y ^ 2)
  return(EcXnorm(c, mu.xy, sigma.xy))
}


PcXnorm <- function(a, b, mu, sigma) {
  #compute E(1(a < X <b)) for X = N(mu,sigma)
  return(pnorm(b, mean = mu, sd = sigma) - pnorm(a, mean = mu, sd = sigma))
}

EcIntnorm <- function(a, b, mu, sigma) {
  #compute int_a^b x pi(x) dx when pi is pdf for N(mu,sigma)
  A = (a - mu) / sigma
  B = (b - mu) / sigma
  return((sigma / sqrt(2 * pi)) * (exp(-0.5 * A ^ 2) - exp(-0.5 * B ^ 2)) + mu *
           (pnorm(B) - pnorm(A)))
}


EEcIntnorm <- function(a, b, mu, sigma, mu.y, sigma.y) {
  #compute E(int_a^b x pi(x) dx) when pi is pdf for N(y-mu,sigma)
  #and y is N(mu.y,sigma.y)
  mu.d  = mu.y - mu
  sigma.d = sqrt(sigma ^ 2 + sigma.y ^ 2)
  return(EcIntnorm(a, b, mu.d, sigma.d))
}

#computes E(phi(X)) when X is N(mu,sigma)
Ephi <- function(mu, sigma) {
  return(exp(-0.5 * mu ^ 2. / (sigma ^ 2 + 1)) / sqrt(2 * pi * (sigma ^ 2 +
                                                                  1)))
}

#computes E(Phi(X)) when X is N(mu,sigma)
EPhi <- function(mu, sigma) {
  return(pnorm(mu / sqrt(1 + sigma ^ 2)))
}
#computes E(X*Phi(X)) when X is N(mu,sigma)
EXPhi <- function(mu, sigma) {
  return(mu * EPhi(mu, sigma) + sigma ^ 2 * Ephi(mu, sigma))
}

#compute E[|X-X'|] when X is N(mu,sigma^2)
Exx <- function(mu, sigma) {
  #X-X' = N(0,2*sigma^2)
  return(Efnorm(0, sqrt(2) * sigma))
}
#compute E[|X-y|] when X is N(mu,sigma^2)
Exy <- function(mu, sigma, y) {
  #X-y = N(mu-y,sigma^2)
  return(Efnorm(mu - y, sigma))
}

#compute E[|X-Y|] when X is N(mu.x,sigma.x^2) and Y is N(mu.y,sigma.y^2)
EExy <-  function(mu.x, sigma.x, mu.y, sigma.y) {
  #E(|X-Y|) = E(E(|X-y| | Y))
  #E(sigma*sqrt(2/pi)*exp(-((mu-y)^2)/(2*sigma^2))+(mu-y)*(1-2*2*pnorm(-(mu-y)/sigma)))
  # X = (mu.x-y)/sigma.x = N((mu.x-mu.y)/sigma.x,sigma.y/sigma.x)
  #E(sigma*2*dnorm(X)+ sigma*X*(1-2*pnorm(-X))
  return(
    2 * sigma.x * Ephi((mu.x - mu.y) / sigma.x, sigma.y / sigma.x)
    + sigma.x * 2 * EXPhi((mu.x - mu.y) / sigma.x, sigma.y / sigma.x) - sigma.x *
      (mu.x - mu.y) / sigma.x
  )
}

Efnorm <- function(mu, sigma) {
  return(sigma * sqrt(2 / pi) * exp(-(mu ^ 2) / (2 * sigma ^ 2)) + mu * (1 -
                                                                           2 * pnorm(-mu / sigma)))
}

Vfnorm <- function(mu, sigma) {
  return(mu ^ 2 + sigma ^ 2 - Efnorm(mu, sigma) ^ 2)
}
