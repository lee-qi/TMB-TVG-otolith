cnt <- 1
while(loopres[[1]][[1]][[1]][[cnt]]=="Not converged") {
cnt <- cnt + 1
}

Npars <- length(loopres[[1]][[1]][[1]][[cnt]]$parsRE)
parnames <- names(loopres[[1]][[1]][[1]][[cnt]]$parsRE)

dir.check <- getwd()
setwd(file.path(dir.check,"Plots"))

every.scenario.med.EpsRE <- array(dim = c(Nsamp, Nmaxa, Nmina, Nyears))

every.scenario.SI50.EpsRE <-
every.scenario.SI95.EpsRE <- array(dim = c(Nsamp, Nmaxa, Nmina, (Nyears * 2)))

every.scenario.med.parsRE <- array(dim = c(Nsamp, Nmaxa, Nmina, Npars))
every.scenario.parsRE <- array(dim = c(Nsamp, Nmaxa, Nmina, Nsims, Npars))

every.scenario.RMSE <-
every.scenario.adjRsqEps <- 
every.scenario.lmsigma <- array(dim = c(Nsamp, Nmaxa, Nmina, Nsims))

every.scenario.med.KeRE <- 
every.scenario.med.WeRE <- 
every.scenario.med.adjRsqEps <- 
every.scenario.converge <- array(dim = c(Nsamp, Nmaxa, Nmina))


for(fish in 1:Nsamp) {
	for(maxa in 1:Nmaxa) {
		for(mina in 1:Nmina) {
			convcount <- 0

			temp.parsRE <- matrix(ncol = Npars, nrow = Nsims)
			temp.KeRE <- temp.WeRE <- matrix(nrow = Nsims, ncol = Nfishvec[fish])
			temp.EpsRE <- matrix(ncol = Nyears, nrow = Nsims)
			temp.adjRsqEps <- vector(length = Nsims)
			temp.RMSE <- vector(length = Nsims)
			temp.lmsigma <- vector(length = Nsims)


			for(sim in 1:Nsims) {
				if(loopres[[fish]][[maxa]][[mina]][[sim]] == "Not converged") {
					convcount <- convcount + 1

					temp.adjRsqEps[sim] <- 
					temp.RMSE[sim] <- 
					temp.lmsigma[sim] <-
					temp.parsRE[sim,]   <- 
					temp.KeRE[sim,]     <- 
					temp.WeRE[sim,]     <- 
					temp.EpsRE[sim,]    <- NA

#					print(paste0("f", fish, "a", maxa, "aa", mina, "s", sim, "/n"))
					
				}

				else {
#					loopres[[fish]][[maxa]][[mina]][[sim]]$trueEps
#					loopres[[fish]][[maxa]][[mina]][[sim]]$estEps
#					loopres[[fish]][[maxa]][[mina]][[sim]]$lmEps
					temp.adjRsqEps[sim] <- loopres[[fish]][[maxa]][[mina]][[sim]]$adjRsqEps
					temp.RMSE[sim] <- loopres[[fish]][[maxa]][[mina]][[sim]]$RMSE
					temp.lmsigma[sim] <- loopres[[fish]][[maxa]][[mina]][[sim]]$lmEps$sigma
					temp.parsRE[sim,] <- loopres[[fish]][[maxa]][[mina]][[sim]]$parsRE
					temp.KeRE[sim,] <- loopres[[fish]][[maxa]][[mina]][[sim]]$KeRE
					temp.WeRE[sim,] <- loopres[[fish]][[maxa]][[mina]][[sim]]$WeRE

					estEps <- loopres[[fish]][[maxa]][[mina]][[sim]]$estEps
					truEps <- loopres[[fish]][[maxa]][[mina]][[sim]]$trueEps

					pos <- truEps - estEps
					neg <- truEps + estEps

					if(sqrt(mean(pos^2) <= sqrt(mean(neg^2))))  {
						temp.EpsRE[sim,] <- pos
					}
					if(sqrt(mean(pos^2) > sqrt(mean(neg^2)))) {
						temp.EpsRE[sim,] <- neg
					}
				
#					loopres[[fish]][[maxa]][[mina]][[sim]]$minage
#					loopres[[fish]][[maxa]][[mina]][[sim]]$maxage
#					loopres[[fish]][[maxa]][[mina]][[sim]]$Nfish
				}
			}

				every.scenario.med.EpsRE[fish,maxa,mina,] <- apply(temp.EpsRE, MARGIN = 2, FUN = median, na.rm = T) #MARGIN = 2 applies to columns, 1 to rows
				every.scenario.SI50.EpsRE[fish,maxa,mina,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = T)
				every.scenario.SI50.EpsRE[fish,maxa,mina,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = T)
				every.scenario.SI95.EpsRE[fish,maxa,mina,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.05, na.rm = T)
				every.scenario.SI95.EpsRE[fish,maxa,mina,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.95, na.rm = T)


				every.scenario.med.parsRE[fish,maxa,mina,] <- apply(temp.parsRE, MARGIN = 2, FUN = median, na.rm = T)
				every.scenario.parsRE[fish,maxa,mina,,] <- temp.parsRE

				every.scenario.med.adjRsqEps[fish,maxa,mina] <- median(temp.adjRsqEps, na.rm = T)
				every.scenario.RMSE[fish,maxa,mina,] <- temp.RMSE
				every.scenario.adjRsqEps[fish,maxa,mina,] <- temp.adjRsqEps
				every.scenario.lmsigma[fish,maxa,mina,] <- temp.lmsigma
				every.scenario.med.KeRE[fish,maxa,mina] <- median(temp.KeRE, na.rm = T)
				every.scenario.med.WeRE[fish,maxa,mina] <- median(temp.WeRE, na.rm = T)

				every.scenario.converge[fish,maxa,mina] <- 1 - signif(convcount / Nsims, digits = 2)

		}
	}
}

for(fish in 1:Nsamp) {
	png(paste(dir.new, "EpsRE", "Nsamp", Nfishvec[fish],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
	par(mfrow=c(Nmaxa,Nmina), mar = c(.3,.3,.3,.3), oma = c(2,2,1,1), cex.axis = 1.3)
	
	for(maxa in 1:Nmaxa) {
		for(mina in 1:Nmina) {
			  plot(0, xlim=c(0, Nyears), ylim=c(-1.5,1.5), xaxs="i", ylab = "", xlab = "", xaxt = 'n', yaxt = 'n', type = "n")
			  polygon(y=every.scenario.SI95.EpsRE[fish,maxa,mina,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
			  polygon(y=every.scenario.SI50.EpsRE[fish,maxa,mina,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
			  lines(y=every.scenario.med.EpsRE[fish,maxa,mina,], x = 1:Nyears, type = "l", col="black", lwd=2)
			  abline(h = 0, col = "black", lwd = 0.5, lty = 2)
			  if(maxa == Nmaxa) {
			  	if(mina == Nmina) {
			  		axis(side = 1, at = c(0, Nyears/2, Nyears), cex.axis = 1.5)
			  	}
			  	else{
			  		axis(side = 1, at = c(0, Nyears/2), cex.axis = 1.5)}
			  	}
			  if(mina == 1) {axis(side = 2, at = c(-1,0,1), cex.axis = 1.5)}
	
		}

	}

	par(mfrow = c(1,1))
	dev.off()
}


for(fish in 1:Nsamp) {
	png(paste(dir.new, "parsRE", "Nsamp", Nfishvec[fish],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
	par(mfrow=c(Nmaxa,Nmina), mar = c(.3,.3,.3,.3), oma = c(2,3,1,1), cex.axis = 1.6)

	for(maxa in 1:Nmaxa) {
		for(mina in 1:Nmina) {
			
			boxplot(every.scenario.parsRE[fish,maxa,mina,,], horizontal = TRUE, 
					names = parnames, outline = F, ylim = c(-0.5,0.5), xaxt = 'n', yaxt = 'n')
			abline(v = 0, col = "black", lwd = 0.5, lty = 2)

			if(maxa == Nmaxa) {axis(side = 1, at = c(-.4, 0, .4), cex.axis = 1.5)}
			if(mina == 1) {axis(side = 2, las = 1, label = c(expression(kappa), expression("w"[infinity]), expression(beta[kappa]), expression(beta["w"[infinity]]), expression(sigma[kappa]), expression(sigma["w"]), expression(sigma["inc"])), 
				at = seq(from = 1, to = length(parnames)))}

		}
	}

	par(mfrow = c(1,1))

	dev.off()
}

for(fish in 1:Nsamp) {
	png(paste(dir.new, "RsqEps", "Nsamp", Nfishvec[fish],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
	par(mfrow=c(Nmaxa,Nmina), mar = c(0.3,0.3,0.3,0.3), oma = c(2,3,1,1), cex.axis = 1.3)

	for(maxa in 1:Nmaxa) {
		for(mina in 1:Nmina) {
			
			boxplot(every.scenario.adjRsqEps[fish,maxa,mina,], horizontal = FALSE,
					outline = F, ylim = c(0,1), xaxt = 'n', yaxt = 'n')

			abline(h = 0.77, col = "black", lty = 2)
			mtext(side = 1, at = 0.75, line = -1, 
					text = bquote(sigma == .(
							signif(median(every.scenario.lmsigma[fish,maxa,mina,], na.rm = T),2))))
			if(mina == 1) {axis(side = 2, at = c(0, 0.5, 0.77, 1), cex.axis = 1.5,
								label = c(0, 0.5, 0.77, 1))}
#			if(mina == 1) {axis(side = 2, las = 1, label = parnames, at = seq(from = 1, to = length(parnames)))}

		}
	}
	dev.off()
}



mina <- 2
png(paste(dir.new, "MinAge", mina, "EpsRE", ".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfcol=c(Nsamp,Nmaxa), mar = c(.3,.3,.3,.3), oma = c(2,3,1,1), cex.axis = 1.3)

for(fish in 1:Nsamp) {	
	for(maxa in 1:Nmaxa) {
		plot(0, xlim=c(0, Nyears), ylim=c(-1.5,1.5),xaxs="i", ylab = "", xlab = "", ann = , xaxt = 'n', yaxt = 'n', type = "n")
		polygon(y=every.scenario.SI95.EpsRE[fish,maxa,mina,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
		polygon(y=every.scenario.SI50.EpsRE[fish,maxa,mina,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
		lines(y=every.scenario.med.EpsRE[fish,maxa,mina,], x = 1:Nyears, type = "l", col="black", lwd=2)
		abline(h = 0, col = "black", lwd = 0.5, lty = 2)
		if(fish == 1) {axis(side = 2, at = c(-1, 0, 1))}
		if(maxa == Nmaxa) {
			if(fish == Nsamp) {
				axis(side = 1, at = c(0, Nyears/2, Nyears), cex.axis = 1.5)
			}
			else{
				axis(side = 1, at = c(0, Nyears/2), cex.axis = 1.5)}
			}

	}
}
par(mfcol = c(1,1))
dev.off()

png(paste(dir.new, "MinAge", mina, "RsqEps", ".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfcol=c(Nsamp,Nmaxa), mar = c(.3,.3,.3,.3), oma = c(2,2,1,1), cex.axis = 1.5)

for(fish in 1:Nsamp) {	
	for(maxa in 1:Nmaxa) {
			
		boxplot(every.scenario.adjRsqEps[fish,maxa,mina,], horizontal = FALSE,
				outline = F, ylim = c(0,1), xaxt = 'n', yaxt = 'n')

		abline(h = 0.77, col = "black", lty = 2)
		mtext(side = 1, at = 0.75, line = -1, 
				text = bquote(sigma == .(
						signif(median(every.scenario.lmsigma[fish,maxa,mina,], na.rm = T),2))))
		if(fish == 1) {axis(side = 2, at = c(0, 0.5, 0.77, 1), cex.axis = 1.5,
							label = c(0, 0.5, 0.77, 1))}
#			if(mina == 1) {axis(side = 2, las = 1, label = parnames, at = seq(from = 1, to = length(parnames)))}

	}
}

par(mfrow = c(1,1))
dev.off()
