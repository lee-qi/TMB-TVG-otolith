setwd("C:/Users/lee qi/Google Drive/TMB_TVG")

dir.check <- getwd()

real.chron <- read.csv("splitnose_master_chronology.csv")

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

meanKvec <- c(0.2, 0.1, .03)
meanWvec <- c(1, 1.5, 2.2)

betaK <- -0.001
betaW <- 0.09

sigK <- 0.05
sigW <- 0.03
sigInc <- 0.003

sigEp <- 1
rho <- 0.5

maxagevec <- c(15, 25, 100)
Nfishvec <- c(20, 50, 100, 200, 500)
minagevec <- seq(from = 0.1, to = 0.7, by = 0.2)

Nyears <- 150 # Extra 1 year of Epsilons because starting increment is also affected by time

Nfish <- 3

#----------------
# Simulating Data
#----------------
intyr <- 1

#nyrs <- round(0.9*maxage)

Nsex <- 1

Nmaxa <- length(maxagevec)
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


rawdat <- read.csv("splitnose_increments.csv", stringsAsFactors = F)
newdat <- rawdat

for(fish in 1:Nfish) {
	newdat$ID[rawdat$ID==unique(rawdat$ID)[fish]] <- fish
}

indiv <- as.numeric(newdat$ID)

#-------------------------
# Start of simulation loop
#-------------------------
f <- 1
s <- 1
aa <- Nmina
simEpsmat <- matrix(ncol = Nyears, nrow = Nmaxa)

#	dev.new()
#	plot(0, type = "n", xlim = range(age), ylim = range(widths),
#		xlab = "Age", ylab = "Width (mm)")
#	cols <- rainbow(Nfish)
#	for(i in 1:Nfish) {
#		temp <- which(indiv == i)
#		lines(x = age[temp], y = widths[temp], lwd = 2, col = cols[i])
#	}
#
dev.new()
par(mfrow = c(2,1), mar = c(1,3,0.2,0.2), oma = c(1,1,1,1))
plot(0, type = "n", xlim = c(1, Nyears), ylim = c(0, 0.2),
	xaxs="i", xaxt = 'n', ylab = "Width (mm)", xlab = "", xpd = NA)
cols <- rainbow(Nmaxa)

for(i in 1:Nfish) {
	tempdat <- subset(rawdat, ID == unique(rawdat$ID)[i])
	lines(x = (tempdat$Calyr - min(rawdat$Calyr) + 1), y = tempdat$width.mm, col = "black", lwd = 2)
}
#	for(i in 1:Nfish) {
#		temp <- which(indiv == i)
#		lines(x = yrs[temp], y = widths[temp], lwd = 2, col = cols[i])
#	}


for(a in 1:Nmaxa) {
	maxage <- maxagevec[a]
	minage <- round(minagevec[aa] * maxage)

	meanK <- meanKvec[a]
	meanW <- meanWvec[a]

	set.seed(8888*as.numeric(paste(s,aa,a,f, sep = "")))
	
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

	
	simEpsmat[a,] <- simEps

	for(i in 1:Nfish) {
		temp <- which(indiv == i)
		lines(x = yrs[temp], y = widths[temp], lwd = 2, col = cols[a])
	}

}
legend("top", lwd = 2, legend = c("Actual Rockfish", "Simulated SL", "Simulated ML", "Simulated LL"), col = c("black", cols))

plot(x = min(real.chron$year):max(real.chron$year), y = real.chron$splitnose_master_chronology_normalized, type = "l", col = "black", lwd = 2,
	xlab = "Year", ylab = "Year Effects (Epsilon)", xpd = NA, ylim = c(-3.5, 2.5))
for(i in 1:Nmaxa) {
	lines(x = min(real.chron$year):max(real.chron$year), y = simEpsmat[i,(1+Nyears - length(unique(real.chron$year))):Nyears], lwd = 2, col = cols[i])
}


##############################
# Plotting against age
##############################
	
dev.new()
#par(mfrow = c(2,1))
plot(0, type = "n", xlim = c(0, maxagevec[Nmaxa]), ylim = c(0, 0.2),
	xlab = "Age", ylab = "Width (mm)")
cols <- rainbow(Nmaxa)


for(a in 1:Nmaxa) {
	Nfish <- 3
	maxage <- maxagevec[a]
	minage <- round(minagevec[aa] * maxage)

	meanK <- meanKvec[a]
	meanW <- meanWvec[a]

	set.seed(8888*as.numeric(paste(s,aa,a,f, sep = "")))
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
	

#	dev.new()
#	plot(0, type = "n", xlim = range(age), ylim = range(widths),
#		xlab = "Age", ylab = "Width (mm)")
#	cols <- rainbow(Nfish)
	for(i in 1:Nfish) {
		temp <- which(indiv == i)
		lines(x = age[temp], y = widths[temp], lwd = 2, col = cols[a])

	}

}

for(i in 1:66) {
	tempdat <- subset(rawdat, ID == unique(rawdat$ID)[i])
	lines(x = tempdat$Age, y = tempdat$width.mm, col = "black", lwd = 2)
}

#legend("topright", )