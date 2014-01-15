source("code/functions.R")
source("code/functionsNew.R")
library(msm)
# Data simulation:
# Simulate ages at death and sexes:
n <- 2000
study <- c(1965, 2000)
minbirth <- 1965
model <- "go"
shape <- "bt"
nb <- 100
thf <- matrix(c(-1.4, 0.65, 0.07, -3.8, 0.2), 1, 5, 
              dimnames = list("f", c("a0", "a1", "c", "b0", "b1")))
thm <- matrix(c(-1.2, 0.7, 0.16, -3.5, 0.23), 1, 5, 
              dimnames = list("m", c("a0", "a1", "c", "b0", "b1")))
#thm <- matrix(c(-0.8, 0.8, 0.16, -3, 0.23), 1, 5, 
#              dimnames = list("m", c("a0", "a1", "c", "b0", "b1")))
class(thf) <- c(model, shape)
class(thm) <- c(model, shape)
thetaReal <- matrix(c(thf, thm), 2, 5, byrow = TRUE, 
                    dimnames = list(c("f", "m"), 
                                    c("a0", "a1", "c", "b0", "b1")))

xv <- seq(0, 25, 0.1)
class(thetaReal) <- c(model, shape)
Fxf <- 1 - CalcSurv(thf, xv)
Fxm <- 1 - CalcSurv(thm, xv)
sr <- 0.5
sexf <- rev(sort(rbinom(n, 1, sr)))
nf <- sum(sexf)
nm <- n - nf
bi <- sort(sample(minbirth:study[2], n, replace = TRUE))
xi <- bi*0
xi[sexf == 1] <- xv[findInterval(runif(nf), Fxf)]
xi[sexf == 0] <- xv[findInterval(runif(nm), Fxm)]
di <- bi + xi

# Build immigrant males:
bi2 <- sort(sample(minbirth:study[2], nm, replace = TRUE))
xi2 <- bi2*0
xi2 <- xv[findInterval(runif(nm), Fxm)]
di2 <- bi2 + xi2

# Simulate dispersal ages for males:
xdmin <- 1.75
lambda <- c(log(1.25), 1)

# Immigrants dispersal age:
xdi2 <- xdmin + rlnorm(nm, lambda[1], lambda[2])
iddisp <- which(xi2 > xdmin & xdi2 < xi2)
bi2 <- bi2[iddisp]
xi2 <- xi2[iddisp]
di2 <- di2[iddisp]
xdi2 <- xdi2[iddisp]
nm2 <- length(xi2)
xar2 <- apply(cbind(xdi2 + runif(nm2), xi2), 1, min)

# Residents dispersal age:
xdi <- xdmin + rlnorm(n, lambda[1], lambda[2])
idd <- which(xdi < xi & sexf == 0)
xdi[xdi >= xi | sexf == 1] <- 0

# Build data table:
Xi <- c(xi, xi2)
Bi <- c(bi, bi2)
Di <- c(di, di2)
Xdi <- c(xdi, xdi2)
Xari <- c(rep(0, n), xar2)
datRaw <- cbind(Xi, Bi, Di, Xdi, Xari)

# Build covariates:
immId <- c(rep(0, n), rep(1, nm2))
sex <- c(sexf, rep(0, nm2))
emigId <- rep(0, n + nm2)
emigId[idd] <- 1
covars <- cbind(sex, immId, emigId)
covs <- cbind(sex, 1 - sex)
colnames(covs) <- c("f", "m")

# Prep data for analysis:
Xli <- Xi - 0.1
Xli[Xi == 0] <- 0
Xli[covars[, 'emigId'] == 1] <- Xdi[covars[, 'emigId'] == 1]
obsDi <- Di * rbinom(length(Xi), 1, 0.3)
idNoDeath <- which(obsDi == 0)
dat <- data.frame(bi = Bi, di = obsDi, xli = Xli, xdi = Xdi, xai = Xari)

# Define priors:
thetaPriorMean <- matrix(c(-2, 0.01, 0, -4, 0.01), 2, 5, byrow = TRUE, 
                     dimnames = dimnames(thetaReal))
thetaPriorSd <- matrix(1, 2, 5, byrow = TRUE, 
                       dimnames = dimnames(thetaReal))
thetaLow <- thetaReal * 0
thetaLow[, c(1, 4)] <- -Inf
e0prior <- CalcEx(thetaPriorMean[1, ])

lambdaPriorMean <- c(0, 2)
lambdaPriorSd <- c(2, 2)
detectPar <- -log(0.00005) / 2
d0prior <- sum((seq(0, 100, 0.1) * dlnorm(seq(0, 100, 0.1), lambdaPriorMean[1], 
                           lambdaPriorMean[2])) * 0.1)

# Define intial values:
xNow <- dat$di - dat$bi
xNow[idNoDeath] <- dat$xli[idNoDeath] + runif(length(idNoDeath), 0.1, 10)
xTrun <- dat$xai
theNow <- thetaReal 
theNow[1, ] <- c(-3, 0.2, 0, -5, 0.01)
theNow[2, ] <- c(-3, 0.2, 0, -5, 0.01)
theCovNow <- CalcCovTheta(theNow)
mortPdfNow <- CalcMortPdf(theCovNow, xNow)
mortPostNow <- CalcMortPost(theNow, theCovNow, xTrun, mortPdfNow)
dispStateNow <- covars[, 'emigId']
agePostNow <- CalcAgePost(xNow, dat$xli, theCovNow, mortPdfNow, dispStateNow)
lamNow <- c(log(2.1), 1.5)
idIm <- which(covars[, 'immId'] == 1)
idEmNow <- which(dispStateNow == 1)
idPotEm <- which(covs[, 'm'] == 1 & dat$xli >= xdmin & covars[, 'immId'] == 0)
dispLikeNow <- CalcDispLike(dat$xli, dat$xai, lamNow, idIm, idEmNow, idPotEm)
lamPostNow <- CalcLamPost(lamNow, dispLikeNow)
dispPostNow <- CalcDispPost(dat$xli, lamNow, dispLikeNow, idEmNow)


# PREP MCMC:
niter <- 25000
nthe <- length(theNow)
namesThe <- paste(rep(colnames(theNow), each = 2), rep(rownames(theNow), 5), 
                  sep = '.')
theMat <- matrix(0, niter, nthe, dimnames = list(NULL, namesThe))
theMat[1, ] <- c(theNow)
lamMat <- matrix(0, niter, 2, dimnames = list(NULL, c("lam1", "lam2")))
lamMat[1, ] <- lamNow
theJumps <- matrix(rep(0.1, nthe), 2, 5, dimnames = dimnames(theNow))
lamJumps <- matrix(c(0.01, 0.01), 1, 2)
theUpdMat <- theMat
lamUpdMat <- lamMat
iterUpd <- 100
updTarg <- 0.25
UpdJumps <- TRUE

# MCMC:
for (iter in 2:niter) {
  # 1. Update mortality parameters:
  for (pp in 1:nthe) {
    theNew <- theNow
    theNew[pp] <- rtnorm(1, theNow[pp], theJumps[pp], low = thetaLow[pp])
    theCovNew <- CalcCovTheta(theNew)
    mortPdfNew <- CalcMortPdf(theCovNew, xNow)
    mortPostNew <- CalcMortPost(theNew, theCovNew, xTrun, mortPdfNew)
    postRatio <- exp(mortPostNew - mortPostNow)
    if (!is.na(postRatio) & postRatio > runif(1)) {
      theNow <- theNew
      theCovNow <- theCovNew
      mortPdfNow <- mortPdfNew
      mortPostNow <- mortPostNew
      if (UpdJumps) theUpdMat[iter, pp] <- 1
    }
  }
  if (any(theUpdMat[iter, ] == 1)) {
    agePostNow <- CalcAgePost(xNow, dat$xli, theCovNow, mortPdfNow, 
                              dispStateNow)
  }
  
  # 2. Update ages at death:
  xNew <- xNow
  xNew[idNoDeath] <- rtnorm(length(idNoDeath), xNow[idNoDeath], 0.1, 
                            low = dat$xli[idNoDeath])
  mortPdfNew <- CalcMortPdf(theCovNow, xNew)
  agePostNew <- CalcAgePost(xNew, dat$xli, theCovNow, mortPdfNew, 
                            dispStateNow)
  
  r <- exp(agePostNew - agePostNow)[idNoDeath]
  z <- runif(length(idNoDeath))
  idUpd <- idNoDeath[r > z]
  idUpd <-idUpd[!(is.na(idUpd))] 
  if (length(idUpd) > 0) {
    xNow[idUpd] <- xNew[idUpd]
    mortPdfNow[idUpd] <- mortPdfNew[idUpd]
    agePostNow[idUpd] <- agePostNew[idUpd]
  }
  mortPostNow <- CalcMortPost(theNow, theCovNow, xTrun, mortPdfNow)

  # 3. update dispersal parameters:
  for (pp in 1:2) {
    lamNew <- lamNow
    lamNew[pp] <- rtnorm(1, lamNow[pp], lamJumps[pp], 
                            lower = c(-Inf, 0)[pp])
    dispLikeNew <- CalcDispLike(dat$xli, dat$xai, lamNew, idIm, idEmNow, 
                                idPotEm)
    lamPostNew <- CalcLamPost(lamNew, dispLikeNew)
    r <- exp(lamPostNew - lamPostNow)
    if (!is.na(r)) {
      if (r > runif(1)) {
        lamNow <- lamNew
        dispLikeNow <- dispLikeNew
        lamPostNow <- lamPostNew
        if (UpdJumps) lamUpdMat[iter, pp] <- 1
      }
    }
  }
  if (any(lamUpdMat[iter, ] == 1)) {
    dispPostNow <- CalcDispPost(dat$xli, lamNow, dispLikeNow, idEmNow)
  }
  
  
  # 4. Update dispersal state:
  dispStateNew <- dispStateNow
  dispStateNew[idPotEm] <- rbinom(length(idPotEm), 1, 0.5)
  idEmNew <- which(dispStateNew == 1)
  dispLikeNew <- CalcDispLike(dat$xli, dat$xai, lamNow, idIm, idEmNew, idPotEm)
  dispPostNew <- CalcDispPost(dat$xli, lamNow, dispLikeNew, idEmNew)
  
  r <- exp(dispPostNew - dispPostNow)[idPotEm]
  z <- runif(length(idPotEm))
  idUpd <- idPotEm[r > z]
  idUpd <-idUpd[!(is.na(idUpd))] 
  if (length(idUpd) > 0) {
    dispStateNow[idUpd] <- dispStateNew[idUpd]
    dispLikeNow[idUpd] <- dispLikeNew[idUpd]
    dispPostNow[idUpd] <- dispPostNew[idUpd]
  }
  idEmNow <- which(dispStateNow == 1)
  lamPostNow <- CalcLamPost(lamNow, dispLikeNow)
  
  # 5. Dynamic Metropolis to update jumps:
  if (UpdJumps) {
    if (is.element(iter/iterUpd,c(1:100))) {
      theJumps <- UpdateJumps(theJumps, theUpdMat, iter, iterUpd, updTarg)
      lamJumps <- UpdateJumps(lamJumps, lamUpdMat, iter, iterUpd, updTarg)
    }
  }
  
  # 6. Fill in the output matrices:
  theMat[iter, ] <- c(theNow)
  lamMat[iter, ] <- lamNow
}

pdf("PlotsTest.pdf", width = 8, height = 12)
par(mfrow = c(6, 2), mar = c(2, 4, 1, 1))
for (i in 1:ncol(theMat)) {
  plot(theMat[, i], type = 'l', main = colnames(theMat)[i], 
       ylim = range(c(theMat[, i], thetaReal[i])))
  abline(h = thetaReal[i], col = 2)
}

for(i in 1:2) {
  plot(lamMat[, i], type = 'l', main = colnames(lamMat)[i], 
       ylim = range(c(lamMat[, i], lambda[i])))
  abline(h = lambda[i], col = 2)
}


# Construct survival and mortality curves:
burnin <- 1000
idkeep <- seq(burnin, niter, 10)
xv <- seq(0, 20, 0.1)
mortList <- lapply(c("f", "m"), function(ss) {
  ids <- which(substr(namesThe, nchar(namesThe), nchar(namesThe)) == ss)
  mort <- apply(theMat[idkeep, ids], 1, function(th) CalcMort(th, xv))
  mortave <- apply(mort, 1, mean)
  mortci <- apply(mort, 1, quantile, c(0.025, 0.975))
  mortfin <- rbind(mortave, mortci)
  return(mortfin)
})

survList <- lapply(c("f", "m"), function(ss) {
  ids <- which(substr(namesThe, nchar(namesThe), nchar(namesThe)) == ss)
  surv <- apply(theMat[idkeep, ids], 1, function(th) CalcSurv(th, xv))
  survave <- apply(surv, 1, mean)
  survci <- apply(surv, 1, quantile, c(0.025, 0.975))
  survfin <- rbind(survave, survci)
  return(survfin)
})

par(mfrow = c(1, 1), mar = c(5, 4, 1, 1))
plot(range(xv), c(0, 2), col = NA, xlab = "Age", ylab = "Mortality")
for (i in 1:2) {
  polygon(c(xv, rev(xv)), c(mortList[[i]][2, ], rev(mortList[[i]][3, ])),
          col = adjustcolor(i, alpha.f= 0.25), border = i)
  lines(xv, CalcMort(thetaReal[i, ], xv), col = i, lwd = 4)
}


alarm()


dev.off()
