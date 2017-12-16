setwd("C:/Users/lee qi/Google Drive/TMB_TVG")
setwd("/Users/LeeQi/Google Drive/TMB_TVG")

dir.check <- getwd()
dir.new <- "Change betaK and betaW"

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

meanK <- c(0.1)
meanW <- c(1.5)

betaKvec <- c(-0.001, -0.01, -0.1, 0.13)
betaWvec <- c(0.001, 0.01, 0.1, -0.12)

#betaK <- -0.001
#betaW <- 0.09

sigK <- 0.05
sigW <- 0.03
sigInc <- 0.003

sigEp <- 1
rho <- 0.5

maxage <- c(25)
Nfishvec <- c(20,50,100)
minage <- 0.5

Nyears <- 150 # Extra 1 year of Epsilons because starting increment is also affected by time

#----------------
# Simulating Data
#----------------
intyr <- 1

#nyrs <- round(0.9*maxage)

Nsex <- 1

NbK <- length(betaKvec)
NbW <- length(betaWvec)
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
	whilecount <- 0
	while(incs[1]<=0 & whilecount < 50) {
		incs[1] <- rnorm(1, temp, sigInc)
		whilecount <- whilecount+1
	}
	if(whilecount == 50) {incs[1] <- 0.0000001}
	tempwid <- intwid + incs[1]

	for(yr in 2:tempNyrs) {
		temp <- (exp(-Kappa[yr]) - 1) * (tempwid - Winf[yr])
		incs[yr] <- rnorm(1, temp, sigInc)
		whilecount <- 0
		while(incs[yr]<=0 & whilecount < 50) {
			incs[yr] <- rnorm(1, temp, sigInc)
			whilecount <- whilecount + 1
		}
		if(whilecount == 50) {incs[yr] <- 0.0000001}
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

loopres <- foreach(bk = 1:NbK) %:%
			foreach(bw = 1:NbW) %:%
			foreach(f = 1:Nsamp) %:%
				foreach(s = 1:Nsims, .packages = c("TMB", "TMBhelper")) %dopar% {

				Nfish <- Nfishvec[f]
				betaK <- betaKvec[bk]
				betaW <- betaWvec[bw]
	
				set.seed(8888*as.numeric(paste(s,f,bk,bw, sep = "")))
				
				sexperfish <- sample(1:Nsex, Nfish, replace = TRUE)
				
				sex <- sexperfish
				baseage <- sample(1:round(0.1*maxage), Nfish, replace = TRUE)
				finage <- sample(minage:maxage, Nfish, replace = TRUE)
				Nyrsperfish <- finage - baseage
				startyrsperfish <- NULL
				for(i in 1:Nfish) {
				
					while(Nyrsperfish[i] < 3) {
						baseage[i] <- sample(1:round(0.1*maxage), 1)
						finage[i] <- sample(minage:maxage, 1)
						Nyrsperfish[i] <- finage[i] - baseage[i]
								}
			
				startyrsperfish[i] <- sample(intyr:((intyr + Nyears - Nyrsperfish[i] - 1)), 1)
			
				if(i==1) {
			#		sex <- rep(sexperfish[i], Nyrsperfish[i])
					startyrsperfish[i] <- intyr
					indiv <- rep(i, Nyrsperfish[i])
					yrs <- seq(from = startyrsperfish[i], to = (startyrsperfish[i] + Nyrsperfish[i] - 1))
			
					age <- seq(from = baseage[i], to = (finage[i] - 1))
				}
				else{
					if(i == Nfish) {
						startyrsperfish[Nfish] <- Nyears - Nyrsperfish[Nfish]
					}
			#		sex <- c(sex, rep(sexperfish[i], Nyrsperfish[i]))
					indiv <- c(indiv, rep(i, Nyrsperfish[i]))
					yrs <- c(yrs, seq(from = startyrsperfish[i], to = (startyrsperfish[i] + Nyrsperfish[i] - 1)))
			
					age <- c(age, seq(from = baseage[i], to = (baseage[i] + Nyrsperfish[i] - 1)))
				}
			}
			
			
			# AR-1
			#simOm <- rnorm(Nyears, 0, sigEp)
			simOm <- rnorm(Nyears, ((-sigEp/2)*((1-rho)/sqrt(1-rho^2))), sigEp)
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
			
			#dev.new()
			#plot(0, type = "n", xlim = range(age), ylim = range(widths),
			#	xlab = "Age", ylab = "Width (mm)")
			#cols <- rainbow(Nfish)
			#for(i in 1:Nfish) {
			#	temp <- which(indiv == i)
			#	lines(x = age[temp], y = widths[temp], lwd = 2, col = cols[i])
			#}
			#
			#dev.new()
			#plot(0, type = "n", xlim = range(yrs), ylim = range(widths),
			#	xlab = "Year", ylab = "Width (mm)")
			#cols <- rainbow(Nfish)
			#for(i in 1:Nfish) {
			#	temp <- which(indiv == i)
			#	lines(x = yrs[temp], y = widths[temp], lwd = 2, col = cols[i])
			#}
			
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
			pars$betaK <- rep(0.1, Nsex)
			pars$betaW <- rep(0.1, Nsex)
			pars$Omega <- rep(0.001, Nyears)
			pars$Ke <- rep(.00001, Nfish)
			pars$We <- rep(.00001, Nfish)
			pars$logsigKe <- log(0.001)
			pars$logsigWe <- log(0.001)
			pars$logsigInc <- log(0.0001)
					
			setwd(dir.check)
			compile("TVG.cpp")
			dyn.load(dynlib("TVG"))
			
			objfun <- MakeADFun(dat, pars, DLL = "TVG", hessian = TRUE, silent = TRUE,
								random = c("Omega", "Ke", "We")
								)
			
			fit1 <- nlminb(objfun$par, objfun$fn, objfun$gr,
							control = list(iter.max = 1e5, eval.max = 1e5))
			
			
			optgrad <- abs(Optimize(objfun, newtonsteps = 2, getsd = FALSE)$diagnostics$final_gradient)
			
			if(sum(which(optgrad>=10e-6))>0) {
				return("Not converged")
			}
			
			if(sum(which(optgrad>=10e-6))==0) {
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

				RMSE <- min(c(sqrt((sum((estEps - trueEps)^2))/Nyears),
							  sqrt((sum((estEps + trueEps)^2))/Nyears)))
			
				lmEps <- summary(lm(trueEps~0+estEps))
			
				replist <- NULL
				replist$objrep <- objrep
				replist$datfile <- dat
				replist$pars <- pars
				replist$objfun <- objfun
				replist$fit <- fit1
				replist$trueEps <- trueEps
				replist$estEps <- estEps
				replist$lmEps <- lmEps
				replist$adjRsqEps <- lmEps$adj.r.squared
				replist$RMSE <- RMSE
				replist$parsRE <- parsRE
				replist$KeRE <- KeRE
				replist$WeRE <- WeRE
				replist$EpsRE <- EpsRE
			
				replist$minage <- minage
				replist$maxage <- maxage
				replist$Nfish <- Nfish
			
				save(replist, file = file.path(getwd(), dir.new, paste0("bk", bk, "bw", bw, "f", f, "s", s,".RData")))
				return(replist)
			}
		}
#	}
#}
#}

loopres <- foreach(bk = 1:NbK) %:%
		foreach(bw = 1:NbW) %:%
		foreach(f = 1:Nsamp) %:%
			foreach(s = 1:Nsims) %dopar% {
				repfile <- file.path(dir.check, dir.new, paste0("bk", bk, "bw", bw, "f", f, "s", s,".RData"))
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

source(file.path(dir.check, "Results_ChangebetKbetW.R"))

