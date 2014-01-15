# Functions:
CalcMort <- function(th, ...) UseMethod("CalcMort")

CalcMort.matrix <- function(th, x) {
  exp(th[, 1] - th[, 2] * x) + th[, 3] + exp(th[, 4] + th[, 5] * x)
}

CalcMort.numeric <- function(th, x) {
  exp(th[1] - th[2] * x) + th[3] + exp(th[4] + th[5] * x)
}

CalcSurv <- function(th, ...) UseMethod("CalcSurv")

CalcSurv.matrix <- function(th, x) {
  exp(exp(th[, 1])/th[, 2] * (exp(-th[, 2] * x) - 1) - th[, 3] * x + 
        exp(th[, 4])/th[, 5] * (1 - exp(th[, 5] * x)))
}

CalcSurv.numeric <- function(th, x) {
  exp(exp(th[1])/th[2] * (exp(-th[2] * x) - 1) - th[3] * x + 
        exp(th[4])/th[5] * (1 - exp(th[5] * x)))
}

CalcCovTheta <- function(th) {
  covs %*% th
}

CalcMortPdf <- function(th, x) {
  log(CalcMort(th, x) * CalcSurv(th, x))
}

CalcMortPost <- function(th, thCov, xt, mortPdf) {
  #sum(mortPdf - log(CalcSurv(thCov, xt))) + 
  sum(mortPdf) + 
    sum(dtnorm(c(th), c(thetaPriorMean), c(thetaPriorSd),
               low = c(thetaLow), log = TRUE))
}

CalcEx <- function(th, minAge = 0, dx = 0.1, xMax = 1000) {
  sum(CalcSurv(th, seq(minAge, xMax, dx)) * dx)
}

CalcAgePrior <- function(x, thp) {
  CalcSurv(thp, x) / e0prior
}

CalcAgePost <- function(x, xl, th, mortPdf, dispState) {
  mortPdf - (1 - dispState) * (detectPar * (x - xl)) + 
    log(CalcAgePrior(x, th))
}

CalcDispLike <- function(xl, xar, lam, idIm, idEm, idPotEm) {
  dispLike <- xl * 0
  dispLike[idIm] <- plnorm(xar[idIm] - xdmin, lam[1], lam[2], log = TRUE)
  idNotEm <- idPotEm[!(idPotEm %in% idEm)]
  dispLike[idNotEm] <- log(1 - plnorm(xl[idNotEm] - xdmin, lam[1], lam[2]))
  dispLike[idEm] <- dlnorm(xl[idEm] - xdmin + 0.1, lam[1], lam[2], log = TRUE)
  return(dispLike)
}

CalcLamPost <- function(lam, dispLike) {
  sum(dispLike) + 
    dnorm(lam[1], lambdaPriorMean[1], lambdaPriorMean[2], log = TRUE) +
    1 / dgamma(lam[2]^2, lambdaPriorSd[1], lambdaPriorSd[2], log = TRUE)  
}


CalcDispPost <- function(xl, lam, dispLike, idEm) {
  dispPost <- log(1 - plnorm(xl - xdmin, lam[1], lam[2]))
#  dispPost[idEm] <- plnorm(xl[idEm] - xdmin + 0.1, lam[1], lam[2], log = TRUE)
  dispPost[idEm] <- dlnorm(xl[idEm] - xdmin + 0.1, lam[1], lam[2], log = TRUE)
  return(dispPost)
}


UpdateJumps <- function(jumps, updMat, iter, iterUpd, updTarg) {
  updRate <- apply(updMat[iter - ((iterUpd - 1):0), ], 2, sum) / iterUpd  
  updRate[updRate == 0] <- 1e-2
  jumps <- jumps * 
    matrix(updRate, nrow(jumps), ncol(jumps)) / updTarg
  return(jumps)
}


