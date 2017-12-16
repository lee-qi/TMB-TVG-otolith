cnt <- 1
while(loopres[[1]][[cnt]]=="Not converged") {
cnt <- cnt + 1
}

dir.check <- getwd()
setwd(file.path(dir.check, "Plots"))
Npars <- length(loopres[[1]][[cnt]]$parsRE)
parnames <- names(loopres[[1]][[cnt]]$parsRE)


every.scenario.med.EpsRE <- array(dim = c(Nmaxa, Nyears))

every.scenario.SI50.EpsRE <-
every.scenario.SI95.EpsRE <- array(dim = c(Nmaxa, (Nyears * 2)))

every.scenario.med.parsRE <- array(dim = c(Nmaxa, Npars))
every.scenario.parsRE <- array(dim = c(Nmaxa, Nsims, Npars))

every.scenario.adjRsqEps <- array(dim = c(Nmaxa, Nsims))

every.scenario.med.KeRE <- 
every.scenario.med.WeRE <- 
every.scenario.med.adjRsqEps <- 
every.scenario.converge <- array(dim = c(Nmaxa))


for(maxa in 1:Nmaxa) {
	convcount <- 0

	Nfish <- Nfishvec[maxa]
	temp.parsRE <- matrix(ncol = Npars, nrow = Nsims)
	temp.KeRE <- temp.WeRE <- matrix(nrow = Nsims, ncol = Nfish)
	temp.EpsRE <- matrix(ncol = Nyears, nrow = Nsims)
	temp.adjRsqEps <- vector(length = Nsims)


	for(sim in 1:Nsims) {
		if(loopres[[maxa]][[sim]] == "Not converged") {
			convcount <- convcount + 1

			temp.adjRsqEps[sim] <- 
			temp.parsRE[sim,]   <- 
			temp.KeRE[sim,]     <- 
			temp.WeRE[sim,]     <- 
			temp.EpsRE[sim,]    <- NA

			
		}

		else {

			temp.adjRsqEps[sim] <- loopres[[maxa]][[sim]]$adjRsqEps
			temp.parsRE[sim,] <- loopres[[maxa]][[sim]]$parsRE
			temp.KeRE[sim,] <- loopres[[maxa]][[sim]]$KeRE
			temp.WeRE[sim,] <- loopres[[maxa]][[sim]]$WeRE

			estEps <- loopres[[maxa]][[sim]]$estEps
			truEps <- loopres[[maxa]][[sim]]$trueEps

			pos <- estEps - truEps
			neg <- -estEps - truEps

			if(sqrt(mean(pos^2) <= sqrt(mean(neg^2))))  {
				temp.EpsRE[sim,] <- pos
			}
			if(sqrt(mean(pos^2) > sqrt(mean(neg^2)))) {
				temp.EpsRE[sim,] <- neg
			}
						

		}
	}

		every.scenario.med.EpsRE[maxa,] <- apply(temp.EpsRE, MARGIN = 2, FUN = median, na.rm = T) #MARGIN = 2 applies to columns, 1 to rows
		every.scenario.SI50.EpsRE[maxa,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = T)
		every.scenario.SI50.EpsRE[maxa,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = T)
		every.scenario.SI95.EpsRE[maxa,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.05, na.rm = T)
		every.scenario.SI95.EpsRE[maxa,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.95, na.rm = T)


		every.scenario.med.parsRE[maxa,] <- apply(temp.parsRE, MARGIN = 2, FUN = median, na.rm = T)
		every.scenario.parsRE[maxa,,] <- temp.parsRE

		every.scenario.med.adjRsqEps[maxa] <- median(temp.adjRsqEps, na.rm = T)
		every.scenario.adjRsqEps[maxa,] <- temp.adjRsqEps
		every.scenario.med.KeRE[maxa] <- median(temp.KeRE, na.rm = T)
		every.scenario.med.WeRE[maxa] <- median(temp.WeRE, na.rm = T)

		every.scenario.converge[maxa] <- 1 - signif(convcount / Nsims, digits = 2)

}

png(paste(dir.new, "EpsRE.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(Nmaxa,1), mar = c(.3,.3,.3,.3), oma = c(2,2,1,1))

for(maxa in 1:Nmaxa) {

		  plot(0, xlim=c(0, Nyears), ylim=c(-.8,.8),xaxs="i", ylab = "", xlab = "", ann = F, xaxt = 'n', yaxt = 'n', type = "n")
		  polygon(y=every.scenario.SI95.EpsRE[maxa,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
		  polygon(y=every.scenario.SI50.EpsRE[maxa,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
		  lines(y=every.scenario.med.EpsRE[maxa,], x = 1:Nyears, type = "l", col="black", lwd=2)
		  abline(h = 0, col = "black", lwd = 0.5, lty = 2)
		  if(maxa == Nmaxa) {axis(side = 1, cex.axis = 1.5)}
		  axis(side = 2, cex.axis = 1.5)

}

par(mfrow = c(1,1))
dev.off()


png(paste(dir.new, "parsRE.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(Nmaxa,1), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))

for(maxa in 1:Nmaxa) {
		
		boxplot(every.scenario.parsRE[maxa,,], horizontal = TRUE,
				names = parnames, outline = F, ylim = c(-.5,.5), xaxt = 'n', yaxt = 'n')
		abline(v = 0, col = "black", lwd = 0.5, lty = 2)

		if(maxa == Nmaxa) {axis(side = 1, at = c(-.4, 0, .4), cex.axis = 1.5)}
		axis(side = 2, las = 1, label = c(expression(kappa), expression("w"[infinity]), expression(beta[kappa]), expression(beta["w"[infinity]]), expression(sigma[kappa]), expression(sigma["w"]), expression(sigma["inc"])),
			 at = seq(from = 1, to = length(parnames)), cex.axis = 1.5)
}

dev.off()

