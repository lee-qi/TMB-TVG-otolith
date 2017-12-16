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
pars$logbetaK <- rep(0.1, Nsex)
pars$logbetaK[1] <- log(pars$logbetaK[1])
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


testgrad <- function(obfun = objfun) {
Optimize(obfun, newtonsteps = 2)
return(1)
}

optres <- testgrad()

# Generating report output
if(!exists("optres")) {
	return("Not converged")
}

if(exists("optres")) {
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