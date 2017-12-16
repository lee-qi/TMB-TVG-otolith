setwd("C:/Users/lee qi/Google Drive/TMB_TVG")

dir.check <- getwd()
dir.new <- "Real params"

require(foreach)
require(doParallel)

Ncores <- detectCores() - 2

cl <- makeCluster(Ncores)
registerDoParallel(cl)

#---------------------
# Initial param values
#---------------------
#meanK <- c(0.227977, 0.249566)
#meanW <- c(26.5548, 24.8187)
#
#betaK <- c(0.005, 0.003)
#betaW <- c(-0.004, -0.002)

#Ep <- rep(0, Nyears)
#logsigEp <- 0.00001

meanK <- 0.02
meanW <- 2.3

betaK <- 1.872263e-27
betaW <- 8.734038e-02

sigK <- 5.026118e-02
sigW <- 3.381055e-02
sigInc <- 3.313147e-03

sigEp <- 1
rho <- 0.8

maxage <- 100
a <- 1
Nfishvec <- c(20, 50, 100, 200, 500)
minagevec <- seq(from = 0.1, to = 0.7, by = 0.2)

Nyears <- 150 # Extra 1 year of Epsilons because starting increment is also affected by time

#----------------
# Simulating Data
#----------------
intyr <- 1

#nyrs <- round(0.9*maxage)

Nsex <- 1

Nmaxa <- 1
Nmina <- length(minagevec)
Nsamp <- length(Nfishvec)
Nsims <- 100

sim_inc <- function(K, W, ep, betK, betW, intep, baseAge, Incsig = sigInc, Ksig = sigK, Wsig = sigW) {
	We <- rnorm(1, 0, Wsig)
	Ke <- rnorm(1, 0, Ksig)

	baseWinf <- W * exp(We)
	baseKappa <- K * exp(Ke)

	Winf <- baseWinf * exp(betW * ep)
	Kappa <- baseKappa * exp(betK * ep)

	tempNyrs <- length(Winf)
	incs <- vector(length = tempNyrs)

	intwid <- (baseWinf * exp(betW * intep)) * (1 - exp(-(baseKappa * exp(betK * intep)) * (baseAge - 1)))
	intwid <- rnorm(1, intwid, sigInc)

	temp <- (exp(-Kappa[1]) - 1) * (intwid - Winf[1])
	incs[1] <- rnorm(1, temp, sigInc)
	while(incs[1]<=0) {
		incs[1] <- rnorm(1, temp, sigInc)
	}
	tempwid <- intwid + incs[1]

	for(yr in 2:tempNyrs) {
		temp <- (exp(-Kappa[yr]) - 1) * (tempwid - Winf[yr])
		incs[yr] <- rnorm(1, temp, sigInc)
		while(incs[yr]<=0) {
			incs[yr] <- rnorm(1, temp, sigInc)
		}
		tempwid <- tempwid + incs[yr]
	}

	res <- NULL

	res$We <- We
	res$Ke <- Ke
	res$intwid <- intwid
	res$incs <- incs
	return(res)
}

#-------------------------
# Start of simulation loop
#-------------------------
#loopres <- list()
#for(f in 1:length(Nfishvec)){
#	for(a in 1:length(maxagevec)){
#		for(aa in 1:length(minagevec)){

loopres <- 	foreach(f = 1:Nsamp) %:%
		foreach(aa = 1:Nmina) %:%

			foreach(s = 1:Nsims) %dopar% {
			Nfish <- Nfishvec[f]
			minage <- round(minagevec[aa] * maxage)

			set.seed(8888*as.numeric(paste(s,aa,f, sep = "")))
		
#			Nyrsperfish <- sample(minage:maxage, Nfish, replace = TRUE)
#			startyrsperfish <- sample(intyr:(intyr+maxage), Nfish, replace = TRUE)
#			endyr <- max(Nyrsperfish + startyrsperfish - 1)
#			styr <- min(startyrsperfish)
#			Nyears <- endyr - styr + 1 + 1 # Extra 1 year of Epsilons because starting increment is also affected by time
			
			#yrs <- (styr:endyr) - styr + 1

#			startyrsperfish <- sample(intyr:(intyr+maxage), Nfish, replace = TRUE)			
			sexperfish <- sample(1:Nsex, Nfish, replace = TRUE)
			
			sex <- sexperfish
			baseage <- NULL
			finage <- NULL
			Nyrsperfish <- NULL
			startyrsperfish <- NULL
			for(i in 1:Nfish) {
				baseage[i] <- sample(1:round(0.1*maxage), 1)
				finage[i] <- sample(minage:maxage, 1)
				Nyrsperfish[i] <- finage[i] - baseage[i]

				while(Nyrsperfish[i] < 3) {
					baseage[i] <- sample(1:round(0.1*maxage), 1)
					finage[i] <- sample(minage:maxage, 1)
					Nyrsperfish[i] <- finage[i] - baseage[i]
				}

				startyrsperfish[i] <- sample(intyr:((intyr + Nyears - Nyrsperfish[i] - 1)), 1)

				if(i==1) {
			#		sex <- rep(sexperfish[i], Nyrsperfish[i])

					indiv <- rep(i, Nyrsperfish[i])
					yrs <- seq(from = startyrsperfish[i], to = (startyrsperfish[i] + Nyrsperfish[i] - 1))

					age <- seq(from = baseage[i], to = (finage[i] - 1))
				}
				else{
			#		sex <- c(sex, rep(sexperfish[i], Nyrsperfish[i]))
					indiv <- c(indiv, rep(i, Nyrsperfish[i]))
					yrs <- c(yrs, seq(from = startyrsperfish[i], to = (startyrsperfish[i] + Nyrsperfish[i] - 1)))

					age <- c(age, seq(from = baseage[i], to = (baseage[i] + Nyrsperfish[i] - 1)))
				}
			}

			
			# AR-1
			simOm <- rnorm(Nyears, 0, sigEp)
			simEps <- NULL
			simEps[1] <- simOm[1]
			
			for(i in 2:Nyears) {
				simEps[i] <- (rho * simEps[i-1]) + (sqrt(1-rho^2)*simOm[i])
			}
			
			simEps <- simEps - mean(simEps)
			
			intwid <- vector(length = Nfish)
			widths <- vector(length = length(age))
			We <- vector(length = Nfish)
			Ke <- vector(length = Nfish)
			for(i in 1:Nfish) {
					vecpos <- which(indiv == i)
					temp <- sim_inc(K = meanK[sex[i]], W = meanW[sex[i]],
											  ep = simEps[yrs[vecpos]+1],
											  intep = simEps[yrs[vecpos][1]],
											  betK = betaK[sex[i]], betW = betaW[sex[i]],
											  baseAge = baseage[i])
					widths[vecpos] <- temp$incs
					intwid[i] <- temp$intwid
					We[i] <- temp$We
					Ke[i] <- temp$Ke
			}
			
#			dev.new()
#			plot(0, type = "n", xlim = range(age), ylim = range(widths),
#				xlab = "Age", ylab = "Width (mm)")
#			cols <- rainbow(Nfish)
#			for(i in 1:Nfish) {
#				temp <- which(indiv == i)
#				lines(x = age[temp], y = widths[temp], lwd = 2, col = cols[i])
#			}
#
#			dev.new()
#			plot(0, type = "n", xlim = range(yrs), ylim = range(widths),
#				xlab = "Year", ylab = "Width (mm)")
#			cols <- rainbow(Nfish)
#			for(i in 1:Nfish) {
#				temp <- which(indiv == i)
#				lines(x = yrs[temp], y = widths[temp], lwd = 2, col = cols[i])
#			}

			#-----------------------
			# Creating lists for TMB
			#-----------------------
			dat <- list()
			dat$indiv <- indiv - 1
			dat$yrs <- yrs
			dat$sex <- sex - 1
			dat$widths <- widths
			dat$age <- age
			dat$Nfish <- Nfish
			dat$Nyears <- Nyears
			dat$Nsex <- Nsex
			
			pars <- list()
			pars$logmeanK <- log(meanK)
			pars$logmeanWinf <- log(meanW)
			pars$logbetaK <- rep(log(0.1), Nsex)
			pars$logbetaW <- rep(log(0.1), Nsex)
			pars$Omega <- rep(0.001, Nyears)
			pars$Ke <- rep(.00001, Nfish)
			pars$We <- rep(.00001, Nfish)
			pars$logsigKe <- log(0.001)
			pars$logsigWe <- log(0.001)
			pars$logsigInc <- log(0.0001)
			
			#pars$betaK <- rep(betaK, Nsex)
			#pars$betaW <- rep(betaW, Nsex)
			#pars$Omega <- simEps
			#pars$Ke <- Ke
			#pars$We <- We
			#pars$logsigKe <- log(sigK)
			#pars$logsigWe <- log(sigW)
			#pars$logsigInc <- log(sigInc)
			
			
			require(TMB)
			setwd(dir.check)
			compile("TVG.cpp")
			dyn.load(dynlib("TVG"))
			
			objfun <- MakeADFun(dat, pars, DLL = "TVG", hessian = TRUE, silent = TRUE,
								random = c("Omega", "Ke", "We")
			#					random = c("Ke", "We")
			#					map = list(
			#								Ke = factor(rep(NA, Nfish)),
			#								We = factor(rep(NA, Nfish)),
			#								logsigWe = factor(NA),
			#								logsigKe = factor(NA),
			#								Omega = factor(rep(NA, Nyears)),
			#								betaK = factor(rep(NA, Nsex)), betaW = factor(rep(NA, Nsex))
			#								logmeanK = factor(rep(5, Nsex)), logmeanWinf = factor(rep(6, Nsex))
			#								)
								)
			
			fit1 <- nlminb(objfun$par, objfun$fn, objfun$gr,
							control = list(iter.max = 1e5, eval.max = 1e5))
			
			# Generating report output
			if(fit1$convergence != 0) {
				return("Not converged")
			}
			
			if(fit1$convergence == 0) {
				objrep <- objfun$report()
				truepars <- c(meanK, meanW, betaK, betaW, sigK, sigW, sigInc)
				estpars <- c(objrep$meanK, objrep$meanWinf,
							 objrep$betaK, objrep$betaW,
							 objrep$sigK, objrep$sigW,
							 objrep$sigInc)
			
				parsRE <- (estpars - truepars) / truepars
				names(parsRE) <- c("meanK", "meanW", "betaK", "betaW", "sigK", "sigW", "sigInc")
			
				trueKe <- Ke
				estKe <- objrep$Ke
			
				KeRE <- (estKe - trueKe) / trueKe
			
				trueWe <- We
				estWe <- objrep$We
			
				WeRE <- (estWe - trueWe) / trueWe
			
				trueEps <- simEps
				estEps <- objrep$Ep
				EpsRE <- (estEps - trueEps) / trueEps
		
				lmEps <- summary(lm(trueEps~0+estEps))
			
				replist <- NULL
				replist$trueEps <- trueEps
				replist$estEps <- estEps
				replist$lmEps <- lmEps
				replist$adjRsqEps <- lmEps$adj.r.squared
				replist$parsRE <- parsRE
				replist$KeRE <- KeRE
				replist$WeRE <- WeRE
				replist$EpsRE <- EpsRE

				replist$minage <- minage
				replist$maxage <- maxage
				replist$Nfish <- Nfish
		
				save(replist, file = file.path(getwd(), dir.new, paste0("f", f, "a", a, "aa", aa, "s", s,".RData")))
				return(replist)
			}
		}
#	}
#}
#}

loopres <- foreach(f = 1:Nsamp) %:%
	foreach(a = 1:Nmaxa) %:%
		foreach(aa = 1:Nmina) %:%
			foreach(s = 1:Nsims) %dopar% {
				repfile <- file.path(getwd(), dir.new, paste0("f", f, "a", 3, "aa", aa, "s", s,".RData"))
				if(file.exists(repfile)) {
					load(file = repfile)
					return(replist)
				}
				if(!file.exists(repfile)) {
					return("Not converged")
				}
			}

stopCluster(cl)


save.image(file = file.path(dir.check, dir.new, "all.RData"))

source(file.path(dir.check, "Results.R"))

objrep <- sdreport(objfun)
objrep$cov.fixed
#summary(objrep)
#objfun$report()
#
#dev.new()
#plot(x = widths, y = objfun$report()$predInc)
#abline(0, 1, col = "red")
#
#dev.new()
#plot(x = exp(Ke), y = exp(objfun$report()$Ke))
#abline(0, 1, col = "red")
#
#dev.new()
#plot(x = We, y = objfun$report()$We)
#abline(0, 1, col = "red")
#
#
#dev.new()
#plot(x = simEps, y = objfun$report()$Ep, xlab = "True Epsilons", ylab = "Estimated Epsilons")
#abline(0, 1, col = "red")
#
##Sex-specific map param_vector
#
##Indiv-specific
#
##Different relationships?
#
##Eps Line
#dev.new()
#plot(objfun$report()$Ep, type = "l", col = "red", lwd = 2,
#	ylim =  range(objfun$report()$Ep, simEps),
#	xlab = "Year", ylab = "Year Effects (Epsilon)")
#lines(simEps, col = "black", lwd = 2, type = "l")
#legend("topright", legend = c("Truth", "Predicted"), lwd = 2, col = c("black", "red"))
#
#lm(simEps~0+objfun$report()$Ep)
#
## Get var-covar matrix
#covar.mat <- solve(objfun$he())
#rownames(covar.mat) <- colnames(covar.mat) <- names(objfun$par)
#write.csv(covar.mat, "VarCovar_Matrix.csv")
#
## Get correlation matrix
#correl.mat <- cov2cor(solve(objfun$he()))
#rownames(correl.mat) <- colnames(correl.mat) <- names(objfun$par)
#write.csv(correl.mat, "Correlations_Matrix.csv")
#
## Get standard errors
#sqrt(diag(solve(objfun$he())))

