cnt <- 1
while(loopres[[1]][[1]][[cnt]]=="Not converged") {
cnt <- cnt + 1
}
Npars <- length(loopres[[1]][[1]][[cnt]]$parsRE)
parnames <- names(loopres[[1]][[1]][[cnt]]$parsRE)

dir.check <- getwd()
setwd(file.path(dir.check,"Plots"))

every.scenario.med.EpsRE <- array(dim = c(Nmaxa, Nvar, Nyears))

every.scenario.SI50.EpsRE <-
every.scenario.SI95.EpsRE <- array(dim = c(Nmaxa, Nvar, (Nyears * 2)))

every.scenario.med.parsRE <- array(dim = c(Nmaxa, Nvar, Npars))
every.scenario.parsRE <- array(dim = c(Nmaxa, Nvar, Nsims, Npars))

every.scenario.adjRsqEps <- 
every.scenario.lmsigma <- array(dim = c(Nmaxa, Nvar, Nsims))

every.scenario.med.KeRE <- 
every.scenario.med.WeRE <- 
every.scenario.med.adjRsqEps <- 
every.scenario.converge <- array(dim = c(Nmaxa, Nvar))


for(maxa in 1:Nmaxa) {
	for(varpar in 1:Nvar) {
		convcount <- 0

		temp.parsRE <- matrix(ncol = Npars, nrow = Nsims)
		temp.KeRE <- temp.WeRE <- matrix(nrow = Nsims, ncol = Nfishvec[maxa])
		temp.EpsRE <- matrix(ncol = Nyears, nrow = Nsims)
		temp.adjRsqEps <- vector(length = Nsims)
		temp.lmsigma <- vector(length = Nsims)


		for(sim in 1:Nsims) {
			if(loopres[[maxa]][[varpar]][[sim]] == "Not converged") {
				convcount <- convcount + 1

				temp.adjRsqEps[sim] <- 
				temp.lmsigma[sim] <-
				temp.parsRE[sim,]   <- 
				temp.KeRE[sim,]     <- 
				temp.WeRE[sim,]     <- 
				temp.EpsRE[sim,]    <- NA

#					print(paste0("f",  "a", maxa, "aa", varpar, "s", sim, "/n"))
					
				}

				else {
#					loopres[[maxa]][[varpar]][[sim]]$trueEps
#					loopres[[maxa]][[varpar]][[sim]]$estEps
#					loopres[[maxa]][[varpar]][[sim]]$lmEps
					temp.adjRsqEps[sim] <- loopres[[maxa]][[varpar]][[sim]]$adjRsqEps
					temp.lmsigma[sim] <- loopres[[maxa]][[varpar]][[sim]]$lmEps$sigma
					temp.parsRE[sim,] <- loopres[[maxa]][[varpar]][[sim]]$parsRE
					temp.KeRE[sim,] <- loopres[[maxa]][[varpar]][[sim]]$KeRE
					temp.WeRE[sim,] <- loopres[[maxa]][[varpar]][[sim]]$WeRE

					estEps <- loopres[[maxa]][[varpar]][[sim]]$estEps
					truEps <- loopres[[maxa]][[varpar]][[sim]]$trueEps
		
					pos <- estEps - truEps
					neg <- -estEps - truEps
		
					if(sqrt(mean(pos^2) <= sqrt(mean(neg^2))))  {
						temp.EpsRE[sim,] <- pos
					}
					if(sqrt(mean(pos^2) > sqrt(mean(neg^2)))) {
						temp.EpsRE[sim,] <- neg
					}
									
#					loopres[[maxa]][[varpar]][[sim]]$varparge
#					loopres[[maxa]][[varpar]][[sim]]$maxage
#					loopres[[maxa]][[varpar]][[sim]]$Nfish
				}
			}

				every.scenario.med.EpsRE[maxa,varpar,] <- apply(temp.EpsRE, MARGIN = 2, FUN = median, na.rm = T) #MARGIN = 2 applies to columns, 1 to rows
				every.scenario.SI50.EpsRE[maxa,varpar,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = T)
				every.scenario.SI50.EpsRE[maxa,varpar,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = T)
				every.scenario.SI95.EpsRE[maxa,varpar,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.05, na.rm = T)
				every.scenario.SI95.EpsRE[maxa,varpar,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.95, na.rm = T)


				every.scenario.med.parsRE[maxa,varpar,] <- apply(temp.parsRE, MARGIN = 2, FUN = median, na.rm = T)
				every.scenario.parsRE[maxa,varpar,,] <- temp.parsRE

				every.scenario.med.adjRsqEps[maxa,varpar] <- median(temp.adjRsqEps, na.rm = T)
				every.scenario.adjRsqEps[maxa,varpar,] <- temp.adjRsqEps
				every.scenario.lmsigma[maxa,varpar,] <- temp.lmsigma
				every.scenario.med.KeRE[maxa,varpar] <- median(temp.KeRE, na.rm = T)
				every.scenario.med.WeRE[maxa,varpar] <- median(temp.WeRE, na.rm = T)

				every.scenario.converge[maxa,varpar] <- 1 - signif(convcount / Nsims, digits = 2)

	}
}

png(paste(dir.new, "EpsRE.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(Nmaxa,Nvar), mar = c(.3,.3,.3,.3), oma = c(2,3,1,1))

for(maxa in 1:Nmaxa) {
	for(varpar in 1:Nvar) {
		  plot(0, xlim=c(0, Nyears), ylim=c(-1.5,1.5),xaxs="i", ylab = "", xlab = "", ann = , xaxt = 'n', yaxt = 'n', type = "n")
		  polygon(y=every.scenario.SI95.EpsRE[maxa,varpar,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
		  polygon(y=every.scenario.SI50.EpsRE[maxa,varpar,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
		  lines(y=every.scenario.med.EpsRE[maxa,varpar,], x = 1:Nyears, type = "l", col="black", lwd=2)
		  abline(h = 0, col = "black", lwd = 0.5, lty = 2)
		  if(maxa == Nmaxa) {
		  	if(varpar == Nvar) {
		  		axis(side = 1, at = c(0, Nyears/2, Nyears), cex.axis = 1.5)
			  	}
			  	else{
			  		axis(side = 1, at = c(0, Nyears/2), cex.axis = 1.5)}
			  	}
		  if(varpar == 1) {axis(side = 2, at = c(-1,0,1), cex.axis = 1.5)}

	}

}

par(mfrow = c(1,1))
dev.off()

png(paste(dir.new, "parsRE.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(Nmaxa,Nvar), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))

for(maxa in 1:Nmaxa) {
	for(varpar in 1:Nvar) {
		
		boxplot(every.scenario.parsRE[maxa,varpar,,], horizontal = TRUE,
				names = parnames, outline = F, ylim = c(-.5,.5), xaxt = 'n', yaxt = 'n')
		abline(v = 0, col = "black", lwd = 0.5, lty = 2)

		if(maxa == Nmaxa) {axis(side = 1, at = c(-.4, 0, .4), cex.axis = 1.5)}
		if(varpar == 1) {axis(side = 2, las = 1, label = c(expression(kappa), expression("w"[infinity]), expression(beta[kappa]), expression(beta["w"[infinity]]), expression(sigma[kappa]), expression(sigma["w"]), expression(sigma["inc"])),
						at = seq(from = 1, to = length(parnames)), cex.axis = 1.5)}


	}
}
dev.off()

png(paste(dir.new, "RsqEps.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(Nmaxa,Nvar), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))

for(maxa in 1:Nmaxa) {
	for(varpar in 1:Nvar) {
		
		boxplot(every.scenario.adjRsqEps[maxa,varpar,], horizontal = FALSE,
				outline = F, ylim = c(0,1), xaxt = 'n', yaxt = 'n')

		abline(h = 0.77, col = "black", lty = 2)
		mtext(side = 1, at = 0.75, line = -1, 
				text = bquote(sigma == .(
						signif(median(every.scenario.lmsigma[maxa,varpar,], na.rm = T),2))))

		if(varpar == 1) {axis(side = 2, at = c(0, 0.5, 0.77, 1), cex.axis = 1.5,
							  label = c(0, 0.5, 0.77, 1))}
#		if(maxa == Nmaxa) {mtext(side = 1, text = sigIncvec[varpar]*1000, padj = 0.5)}

	}
}
dev.off()
