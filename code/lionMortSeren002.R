# Serengeti lions:
rm(list = ls())

library(msm)
library(RColorBrewer)
if ("fernando" %in% list.files("/Users/")) {
  setwd("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLionsGithub/compLionMort/")
  seren <- read.csv("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/data/serengeti/seren.csv", header = TRUE)
} else {
  setwd("/Users/Viktualia/Documents/GitHub/compLionMort")
  seren <- read.csv("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/data/serengeti/seren.csv", header = TRUE)
}


# Source functions:
source("code/functions.R")

# want diagnostic plots?:
plotInd <- TRUE

### Extract variables:
# number of observations
n <- nrow(seren)
# time objects (birth, first, last, death)
study <- julian(as.Date(c("1966-04-01", "2013-08-01")), 
                origin = as.POSIXct("1950-01-01"))
birth <- julian(as.Date(seren[, "birthDate"]), 
                origin = as.POSIXct("1950-01-01")) # for age at birth
first <- julian(as.Date(seren[, "firstSeenDate"]), 
               origin = as.POSIXct("1950-01-01")) # for age at immigration
last <- julian(as.Date(seren[, "lsDate"]), # for age at emigration/death
               origin = as.POSIXct("1950-01-01"))
death <- last
death[is.na(seren$dead)] <- NA
# minimum dispersal age
# in Pusey & Packer 1987 all m disp by age 4.2, min 1.8 (1 of 12)
# Elliot et al (subm.) all m disp by 3.75, min 1.66 (no m survived younger than 2.6)
minDispAge <- 1.5
# sex
sex <- as.character(seren[, 'sex'])
# indicator and index for open fate after last seen ( 0: observed death | alive, 1: open fate)
unknownFate <- rep(0, n)
unknownFate[is.na(death)] <- 1
unknownFate[seren$alive == 1] <- 0
idNoDeath <- which(unknownFate == 1)
# left truncation
idLeftTr <- which(birth < study[1] | !is.na(seren$firstSeenDate))
ageTrunc <- apply(cbind(study[1] - birth, 0), 1, max) / 365.25  # born before study start
ageTrunc[!is.na(first)] <-
  c(first - birth)[!is.na(first)] / 365.25  # immigrants
# age at last seen
ageToLast <- (last - birth) / 365.25
# age at first seen (0 for natives, ageTrunc for immigrants)
ageToFirst <- rep(0, n)
ageToFirst[!is.na(first)] <- ageTrunc[!is.na(first)] 
# index for unknown sex
idNoSex <- which(sex == "u" | sex == "x")
# probability newborns female
probFem <- 0.45

### Propose (initial) parameter values:
# Mortality
model <- "go"; shape <- "bt"
ncovs <- 2
defPars <- SetDefaultTheta()
npars <- defPars$length * ncovs
names <- c("f", "m")
thetaNames <- paste(rep(defPars$name, ncovs), 
                    rep(names, each = defPars$len), sep = '.')
thetaStart <- matrix(defPars$start, nrow = ncovs, ncol = defPars$length, 
                     byrow = TRUE, dimnames = list(names, 
                                                   defPars$name))
class(thetaStart) <- c(model, shape)
# Dispersal pars and state indices
lambdaStart <- c(log(2), 1)
lambdaJump <- c(0.2, 0.2)
nDispPars <- length(lambdaStart)
idIMstart <- which(sex == "m" & !is.na(first) & ageToFirst >= minDispAge)
idEMstart <- which(sex == "m" & unknownFate == 1 & ageToLast >= minDispAge)
idSTstart <- (1:n)[!((1:n) %in% idEMstart)]  # not emigrator
# Ages
xStart <- c(last - birth) / 365.25
xStart[idNoDeath] <- xStart[idNoDeath] + sample(1:10, length(idNoDeath), 
                                                replace = TRUE)
# non-dectection probability (if alive and in study area)
detectPar <- -log(0.00005) / 2
# Initial sexes
sexFemStart <- rep(1, n)
sexFemStart[sex == "m"] <- 0
sexFemStart[idNoSex] <- rbinom(length(idNoSex), 1, probFem)
covarsStart  <- matrix(c(sexFemStart, 1 - sexFemStart), n, ncovs, dimnames = list(NULL, names))
# Dispersal state indicator
dispStateStart <- rep(0, n)
dispStateStart[idEMstart] <- 1

### For calculating priors:
xv <- seq(0, 100, 0.1)
thPrior <- matrix(defPars$priorMean, length(xv), defPars$length, byrow = TRUE)
class(thPrior) <- c(model, shape)
priorLamMean <- c(log(3), 2)
priorLamSd <- c(1, 1)
exPrior <- sum(CalcSurv(thPrior, xv) * 0.1)
dispStatePrior <- 0.7

### Build jumps matrix:
jumpMatStart <- matrix(defPars$jump, ncovs, defPars$length, byrow = TRUE,
                       dimnames = dimnames(thetaStart))

### Number of iterations and runs
niter <- 10000
niterRun <- 5

### Run dynamic Metropolis to find jumps
UpdJumps <- TRUE
niter <- 10000
burnin <- 5001
thin <- 20
keep <- seq(burnin, niter, thin)
outJump <- RunMCMC(1)

par(mfcol = c(5, 2), mar = c(3, 3, 1, 1))
for (i in 1:10) plot(outJump$pars[, i], type = 'l', 
                     main = colnames(outJump$pars)[i])

# Run MCMC:
UpdJumps <- FALSE
jumpMatStart <- outJump$jumpsMort
lambdaJump <- outJump$jumpsDisp
rm("outJump")
niter <- 40000
burnin <- 20001
thin <- 40
keep <- seq(burnin, niter, thin)

nsim <- 4
ncpus <- 4
require(snowfall)
sfInit(parallel = TRUE, cpus = ncpus)
sfExport(list = c(ls(), ".Random.seed"))
sfLibrary(msm, warn.conflicts = FALSE)
out <- sfClusterApplyLB(1:nsim, RunMCMC)
sfStop()
  
# plot the mean survival curves

if(plotInd){
  if ("fernando" %in% list.files("/Users/")) {
    pdf("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/plots/hwange04Nov.pdf")
  } else {
    pdf("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/plots/hwange04Nov.pdf")
  }
  
#  keep <- seq(15001, niter, 20)
  parMat <- out[[1]]$pars[keep, ]
  for (i in 2:nsim) {
    parMat <- rbind(parMat, out[[i]]$pars[keep, ])
  }
  
  parQuant <- cbind(apply(parMat, 2, mean), apply(parMat, 2, sd),
                    t(apply(parMat, 2, quantile, c(0.025, 0.975))))
  colnames(parQuant) <- c("Mean", "SE", "2.5%", "97.5%")
  
  thetaFemNew <-  matrix(parQuant[1:5, 1], 1, 5)
  colnames(thetaFemNew) <- thetaNames[1:5]
  thetaMalNew <-  matrix(parQuant[6:10, 1], 1, 5)
  colnames(thetaMalNew) <- thetaNames[6:10]
  xv <- seq(0, 25, 0.1)
  class(thetaMalNew) <- c(model, shape)
  class(thetaFemNew) <- class(thetaMalNew)
  
  ySurvFemNew <- CalcSurv(thetaFemNew, xv)
  ySurvMalNew <- CalcSurv(thetaMalNew, xv)
  
  
  plot(xv, ySurvFemNew, col = 2, lwd = 2, lty = 1, type = "l")
  lines(xv, ySurvMalNew, col = 4, lwd = 2, lty = 1)
  legend("topright", c("females", "males"),
         lwd = c(2,2), lty = c(1, 1), col = c(2,4))
  dev.off()
}

rm(list = setdiff(ls(), c("out", "nsim", "niter", "model", "shape", "ncovs", "names", "npars", 
                          "defPars", "thetaNames")))

if ("fernando" %in% list.files("/Users/")) {
  save.image("/Users/fernando/FERNANDO/PROJECTS/1.ACTIVE/JuliaLions/results/outputHwang04Nov2.Rdata")
} else {
  save.image("/Users/Viktualia/Dropbox/Projects/008_LionSexDiffMort/JuliaLions/results/outputHwang04Nov2.Rdata")
}

pdf("results/trace009.pdf", width = 15, height = 10)
par(mfrow = c(ncovs, defPars$length))
for (i in 1:npars) {
  for (sim in 1:nsim) {
    if (sim == 1) plot(out[[sim]]$pars[ ,i], type = 'l', main = thetaNames[i])
    if (sim != 1) lines(out[[sim]]$pars[ ,i], col = brewer.pal(npars-1, "Set1")[sim])
  }
}
dev.off()

# plot densities male and female separate plots:
#pdf("results/parDens005.pdf", width = 10, height = 10)
#par(mfrow = c(2, 5))
#for (i in 1:10) plot(density(parList[-c(1:1000), i]), type = 'l', 
#                    main = thetaNames[i], lwd = 3)
#dev.off()

pdf("results/parDens005.pdf", width = 12, height = 5)
par(mfrow = c(1, defPars$length))
for (i in 1:defPars$length) {
  plot(density(parList[-c(1:1000), i]), type = 'l', 
       main = c("a0", "a1", "c", "b0", "b1") [i], 
       lwd = 3, xlim = range(parList[-c(1:1000), c(i, i+2)]),
       ylim = c(0, max (c(max(density(parList[-c(1:1000), i])[[2]]), 
                          max(density(parList[-c(1:1000), 
                                              i+defPars$length])[[2]])))))
  lines(density(parList[-c(1:1000), i + 5]), col = 'dark green', lwd = 3)
  if (i == defPars$length) legend("topright", legend = c("female", "male"), 
                                  col = c('black', 'dark green'), lwd = c(3, 3))
}
dev.off()

# Calculate mean and quantiles for parameter estimates:
keep <- seq(1000, niter, 20)
parMat <- out[[1]]$pars[keep, ]
for (i in 2:nsim) {
  parMat <- rbind(parMat, out[[i]]$pars[keep, ])
}

parQuant <- cbind(apply(parMat, 2, mean), apply(parMat, 2, sd),
                  t(apply(parMat, 2, quantile, c(0.025, 0.975))))
colnames(parQuant) <- c("Mean", "SE", "2.5%", "97.5%")

mortList <- list()
survList <- list()
xv <- seq(0, 25, 0.1)
parnames <- colnames(parMat)
for (nc in 1:ncovs) {
  idpars <- which(substr(parnames, nchar(parnames), nchar(parnames)) == names[nc])
  thetamat <- parMat[, idpars]
  nPar <- ncol(thetamat)
  mort <- apply(thetamat, 1, function(tht) {
    tht2 <- matrix(tht, length(xv), nPar, byrow = TRUE)
    class(tht2) <- c(model, shape)
    CalcMort(tht2, xv)})
  mortList[[names[nc]]] <- cbind(apply(mort, 1, mean), 
                                 t(apply(mort, 1, quantile, c(0.025, 0.975))))
  surv <- apply(thetamat, 1, function(tht) {
    tht2 <- matrix(tht, length(xv), nPar, byrow = TRUE)
    class(tht2) <- c(model, shape)
    CalcSurv(tht2, xv)})
  survList[[names[nc]]] <- cbind(apply(surv, 1, mean), 
                                 t(apply(surv, 1, quantile, c(0.025, 0.975))))
  
}

rangesSurv <- sapply(names, function(nn) which(survList[[nn]][, 1] < 0.01)[1])

rangesMort <- sapply(names, function(nn) range(mortList[[nn]][1:rangesSurv[nn], ]))

par(mfrow = c(2, 1), mar = c(4, 5, 1, 1), family = 'serif')
plot(range(xv[1:max(rangesSurv)]), c(0, 1), col = NA, frame.plot = FALSE, 
     xlab = "", ylab = expression(paste("Survival, ", italic(S(x)))))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(survList[[i]][1:rangesSurv[i], 2], 
            rev(survList[[i]][1:rangesSurv[i], 3])), 
          col = c('dark red', 'dark green')[i], 
          border =  c('dark red', 'dark green')[i])
  lines(xv[1:rangesSurv[i]], survList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}
legend('topright', names, pch = 15, cex = 2, col = c('dark red', 'dark green'),
       bty = 'n')

plot(range(xv[1:max(rangesSurv)]), c(0, max(rangesMort)), col = NA, 
     frame.plot = FALSE, xlab = expression(paste("Age, ", italic(x))), 
     ylab = expression(paste("Mortality, ", mu(italic(x)))))
for (i in 1:2) {
  polygon(c(xv[1:rangesSurv[i]], rev(xv[1:rangesSurv[i]])), 
          c(mortList[[i]][1:rangesSurv[i], 2], 
            rev(mortList[[i]][1:rangesSurv[i], 3])), 
          col = c('dark red', 'dark green')[i], 
          border =  c('dark red', 'dark green')[i])
  lines(xv[1:rangesSurv[i]], mortList[[i]][1:rangesSurv[i], 1], col = 'white',
        lwd = 2)
}



