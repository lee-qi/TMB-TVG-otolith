cnt <- 1
while(loopres[[1]][[1]][[1]][[cnt]]=="Not converged") {
cnt <- cnt + 1
}
Npars <- length(loopres[[1]][[1]][[1]][[cnt]]$parsRE)
parnames <- names(loopres[[1]][[1]][[1]][[cnt]]$parsRE)

dir.check <- getwd()
setwd(file.path(dir.check,"Plots"))

every.scenario.med.EpsRE <- array(dim = c(NbK, NbW, Nsamp, Nyears))

every.scenario.SI50.EpsRE <-
every.scenario.SI95.EpsRE <- array(dim = c(NbK, NbW, Nsamp, (Nyears * 2)))

every.scenario.med.parsRE <- array(dim = c(NbK, NbW, Nsamp, Npars))
every.scenario.parsRE <- array(dim = c(NbK, NbW, Nsamp, Nsims, Npars))

every.scenario.adjRsqEps <- 
every.scenario.lmsigma <- array(dim = c(NbK, NbW, Nsamp, Nsims))

every.scenario.med.KeRE <- 
every.scenario.med.WeRE <- 
every.scenario.med.adjRsqEps <- 
every.scenario.converge <- array(dim = c(NbK, NbW, Nsamp))


	for(bk in 1:NbK) {
		for(bw in 1:NbW) {
			for(f in 1:Nsamp) {
			convcount <- 0

			temp.parsRE <- matrix(ncol = Npars, nrow = Nsims)
			temp.KeRE <- temp.WeRE <- matrix(nrow = Nsims, ncol = Nfishvec[f])
			temp.EpsRE <- matrix(ncol = Nyears, nrow = Nsims)
			temp.adjRsqEps <- vector(length = Nsims)
			temp.lmsigma <- vector(length = Nsims)


			for(sim in 1:Nsims) {
				if(loopres[[bk]][[bw]][[f]][[sim]] == "Not converged") {
					convcount <- convcount + 1

					temp.adjRsqEps[sim] <- 
					temp.lmsigma[sim] <- 
					temp.parsRE[sim,]   <- 
					temp.KeRE[sim,]     <- 
					temp.WeRE[sim,]     <- 
					temp.EpsRE[sim,]    <- NA

#					print(paste0("f",  "a", bk, "aa", bw,f "s", sim, "/n"))
					
				}

				else {
#					loopres[[bk]][[bw]][[f]][[sim]]$trueEps
#					loopres[[bk]][[bw]][[f]][[sim]]$estEps
#					loopres[[bk]][[bw]][[f]][[sim]]$lmEps
					temp.adjRsqEps[sim] <- loopres[[bk]][[bw]][[f]][[sim]]$adjRsqEps
					temp.lmsigma[sim] <- loopres[[bk]][[bw]][[f]][[sim]]$lmEps$sigma
					temp.parsRE[sim,] <- loopres[[bk]][[bw]][[f]][[sim]]$parsRE
					temp.KeRE[sim,] <- loopres[[bk]][[bw]][[f]][[sim]]$KeRE
					temp.WeRE[sim,] <- loopres[[bk]][[bw]][[f]][[sim]]$WeRE

					estEps <- loopres[[bk]][[bw]][[f]][[sim]]$estEps
					truEps <- loopres[[bk]][[bw]][[f]][[sim]]$trueEps
		
					pos <- estEps - truEps
					neg <- -estEps - truEps
		
					if(sqrt(mean(pos^2) <= sqrt(mean(neg^2))))  {
						temp.EpsRE[sim,] <- pos
					}
					if(sqrt(mean(pos^2) > sqrt(mean(neg^2)))) {
						temp.EpsRE[sim,] <- neg
					}

#					loopres[[bk]][[bw]][[f]][[sim]]$wege
#					loopres[[bk]][[bw]][[f]][[sim]]$kege
#					loopres[[bk]][[bw]][[f]][[sim]]$Nfish
				}
			}

				every.scenario.med.EpsRE[bk,bw,f,] <- apply(temp.EpsRE, MARGIN = 2, FUN = median, na.rm = T) #MARGIN = 2 applies to columns, 1 to rows
				every.scenario.SI50.EpsRE[bk,bw,f,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.25, na.rm = T)
				every.scenario.SI50.EpsRE[bk,bw,f,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.75, na.rm = T)
				every.scenario.SI95.EpsRE[bk,bw,f,] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.05, na.rm = T)
				every.scenario.SI95.EpsRE[bk,bw,f,(Nyears*2):(Nyears+1)] <- apply(temp.EpsRE, MARGIN = 2, FUN = quantile, probs = 0.95, na.rm = T)


				every.scenario.med.parsRE[bk,bw,f,] <- apply(temp.parsRE, MARGIN = 2, FUN = median, na.rm = T)
				every.scenario.parsRE[bk,bw,f,,] <- temp.parsRE

				every.scenario.med.adjRsqEps[bk,bw,f] <- median(temp.adjRsqEps, na.rm = T)
				every.scenario.adjRsqEps[bk,bw,f,] <- temp.adjRsqEps
				every.scenario.lmsigma[bk,bw,f,] <- temp.lmsigma
				every.scenario.med.KeRE[bk,bw,f] <- median(temp.KeRE, na.rm = T)
				every.scenario.med.WeRE[bk,bw,f] <- median(temp.WeRE, na.rm = T)

				every.scenario.converge[bk,bw,f] <- 1 - signif(convcount / Nsims, digits = 2)
			}

	}
}

finNbK <- NbK
finNbW <- NbW
NbK <- finNbK - 1
NbW <- finNbW - 1
for(f in 1:Nsamp) {
	png(paste(dir.new, "EpsRE", "Nsamp", Nfishvec[f],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
	par(mfrow=c(NbK,NbW), mar = c(.3,.3,.3,.3), oma = c(2,2,1,1))
	
	for(bk in 1:NbK) {
		for(bw in 1:NbW) {
			  plot(0, xlim=c(0, Nyears), ylim=c(-1.5,1.5),xaxs="i", ylab = "", xlab = "", ann = , xaxt = 'n', yaxt = 'n', type = "n")
			  polygon(y=every.scenario.SI95.EpsRE[bk,bw,f,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
			  polygon(y=every.scenario.SI50.EpsRE[bk,bw,f,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
			  lines(y=every.scenario.med.EpsRE[bk,bw,f,], x = 1:Nyears, type = "l", col="black", lwd=2)
			  abline(h = 0, col = "black", lwd = 0.5, lty = 2)
			  if(bk == NbK) {
			  	if(bw==NbW) {
			  		axis(side = 1, at = c(0, Nyears/2, Nyears), cex.axis = 1.5)}
			  	else{axis(side = 1, at = c(0, Nyears/2), cex.axis = 1.5)}
			  }
			  if(bw == 1) {axis(side = 2, at = c(-1,0,1), cex.axis = 1.5)}
	
		}

	}

	par(mfrow = c(1,1))
	dev.off()
}

for(f in 1:Nsamp) {
	png(paste(dir.new, "parsRE", "Nsamp", Nfishvec[f],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
	par(mfrow=c(NbK,NbW), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))

	for(bk in 1:NbK) {
		for(bw in 1:NbW) {
			
			boxplot(every.scenario.parsRE[bk,bw,f,,], horizontal = TRUE,
					names = parnames, outline = F, ylim = c(-0.5,0.5), xaxt = 'n', yaxt = 'n')
			abline(v = 0, col = "black", lwd = 0.5, lty = 2)

			if(bk == NbK) {axis(side = 1, at = c(-0.4, 0, 0.4), cex.axis = 1.5)}
			if(bw == 1) {axis(side = 2, las = 1, label = c(expression(kappa), expression("w"[infinity]), expression(beta[kappa]), expression(beta["w"]), expression(sigma[kappa]), expression(sigma["w"[infinity]]), expression(sigma["inc"])), 
				at = seq(from = 1, to = length(parnames)), cex.axis = 1.5)}
		}


	}
	dev.off()
}


for(f in 1:Nsamp) {
png(paste(dir.new, "RsqEps", "Nsamp", Nfishvec[f],".png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(NbK,NbW), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))

for(bk in 1:NbK) {
	for(bw in 1:NbW) {
		
		boxplot(every.scenario.adjRsqEps[bk,bw,f,], horizontal = FALSE,
				outline = F, ylim = c(0,1), xaxt = 'n', yaxt = 'n')

		abline(h = 0.77, col = "black", lty = 2)
		mtext(side = 1, at = 0.75, line = -1, 
				text = bquote(sigma == .(
						signif(median(every.scenario.lmsigma[bk,bw,f,], na.rm = T),2))))

		if(bw == 1) {axis(side = 2, at = c(0, 0.5, 0.77, 1), cex.axis = 1.5,
						  label = c(0, 0.5, 0.77, 1))}
#		if(bw == 1) {axis(side = 2, las = 1, label = parnames, at = seq(from = 1, to = length(parnames)))}
		}

	}
	dev.off()
}


bk <- finNbK
bw <- finNbW
#dev.new()
#par(mfrow=c(1,Nsamp), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))
#
#for(f in 1:Nsamp) {
#		boxplot(every.scenario.parsRE[bk,bw,f,,], horizontal = TRUE,
#				names = parnames, outline = F, ylim = c(-0.5,0.5), xaxt = 'n', yaxt = 'n')
#		abline(v = 0, col = "black", lwd = 0.5, lty = 2)
#
#		axis(side = 1, at = c(-0.5, 0, 0.5))
#		if(f == 1) {axis(side = 2, las = 1, label = c(expression(kappa), expression("w"[infinity]), expression(beta[kappa]), expression(beta["w"]), expression(sigma[kappa]), expression(sigma["w"]), expression(sigma["inc"])), 
#			at = seq(from = 1, to = length(parnames)))}
#}

png(paste(dir.new, "truebets", "RsqEps.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(1,Nsamp), mar = c(0.2,0.2,0.2,0.2), oma = c(2,3,1,1))
for(f in 1:Nsamp) {
		boxplot(every.scenario.adjRsqEps[bk,bw,f,], horizontal = FALSE,
				outline = F, ylim = c(0,1), xaxt = 'n', yaxt = 'n')

		abline(h = 0.77, col = "black", lty = 2)
		mtext(side = 1, at = 0.75, line = -1, 
				text = bquote(sigma == .(
						signif(median(every.scenario.lmsigma[bk,bw,f,], na.rm = T),2))))
		if(f == 1) {axis(side = 2, at = c(0, 0.5, 0.77, 1), cex.axis = 1.5,
						 label = c(0, 0.5, 0.77, 1))}
}
dev.off()

png(paste(dir.new, "truebets", "EpsRE.png", sep = "_"), width = 7, height = 7, res = 600, units = 'in')
par(mfrow=c(1,Nsamp), mar = c(.3,.3,.3,.3), oma = c(2,3,1,1))
for(f in 1:Nsamp) {
	plot(0, xlim=c(0, Nyears), ylim=c(-1.5,1.5),xaxs="i", ylab = "", xlab = "", ann = , xaxt = 'n', yaxt = 'n', type = "n")
	polygon(y=every.scenario.SI95.EpsRE[bk,bw,f,],x=c(1:Nyears,Nyears:1), col = "gray60", border = NA)
	polygon(y=every.scenario.SI50.EpsRE[bk,bw,f,],x=c(1:Nyears,Nyears:1),col="gray75", border = NA)
	lines(y=every.scenario.med.EpsRE[bk,bw,f,], x = 1:Nyears, type = "l", col="black", lwd=2)
	abline(h = 0, col = "black", lwd = 0.5, lty = 2)
	if(f!=Nsamp) {
		axis(side = 1, at = c(0, Nyears/2), cex.axis = 1.5)
	}
	if(f==Nsamp) {
		axis(side = 1, at = c(0, Nyears/2, Nyears), cex.axis = 1.5)
	}
	  if(f == 1) {axis(side = 2, cex.axis = 1.5)
	}
}
dev.off()

