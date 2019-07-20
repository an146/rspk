#install.package('RcppFaddeeva')
library('RcppFaddeeva')
voigt <- function(xval, xpeak, fG, fL, height) {
	sigma <- fG / (2 * (2*log(2))^0.5)
	gamma <- fL / 2
	z <- complex(real = xval - xpeak, imaginary = gamma) / (sqrt(2) * sigma)
	peakZ <- complex(imaginary = gamma) / (sqrt(2) * sigma)
	normFactor <- height / (Re(Faddeeva_w(peakZ)) / (sigma * sqrt(2 * pi)))
	return(Re(Faddeeva_w(z)) / (sigma * sqrt(2 * pi)) * normFactor)
}

voigt2 <- function(xval, peak) {
	return(voigt(xval, peak[1], peak[2], peak[3], peak[4]))
}

voigt3 <- function(xval, peaks) {
	#TODO: optimize, iterate functionally
	ret <- numeric(length(xval))
	if (is.null(peaks)) {
	} else if (is.null(dim(peaks))) {
		ret <- voigt2(xval, peaks)
	} else {
		for (i in 1:dim(peaks)[1]) {
			ret <- ret + voigt2(xval, peaks[i,])
		}
	}
	return(ret)
}

accumerror <- function(errors) {
	ret <- 0
	for (e in errors)
		ret <- ret + e * e
	return(ret)
}

voigterror <- function(data, peaks) {
	v <- voigt3(data$x, peaks)
	return(accumerror(data$y - v))
}

voigtlearn <- function(data, peaks, modx, modh, clearn) {
	deriv_coeff = 1e-9
	#learn_coeff = c(1e-6, 1e-6, 1e-6, 1e-6)
	learn_coeff = c(clearn, clearn, clearn, clearn)
	if (is.null(dim(peaks)))
		peaks <- matrix(peaks, ncol=4, nrow=1)
	nextpeaks <- peaks
	xval <- data$x
	yval <- data$y
	eref <- voigterror(data, peaks)
	if (modx)
		startj <- 1
	else
		startj <- 2
	if (modh)
		endj <- 4
	else
		endj <- 3
	for (i in 1:dim(peaks)[1]) {
		for (j in startj:endj) {
			peaks1 <- peaks
			peaks1[i, j] <- peaks1[i, j] + deriv_coeff
			deriv <- (voigterror(data, peaks1) - eref) / deriv_coeff
			nextpeaks[i, j] <- nextpeaks[i, j] - deriv * learn_coeff[j]
		}
	}
	return(nextpeaks)
}

voigtlearn2 <- function(data, peaks, count=1, modx=FALSE, modh=FALSE, clearn=1e-7) {
	for (i in 1:count) {
		peaks <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
	}
	return(peaks)
}

extractpeak <- function(data, peaks) {
	v <- voigt3(data$x, peaks)
	yval <- data$y - v
	maxrow <- 0
	maxvalue <- 0
	for (i in 1:dim(data)[1]) {
		if (data[i, 2] > maxvalue) {
			maxvalue <- data[i, 2]
			maxrow <- i
		}
	}
	newx <- data[maxrow, 1]
	newy <- data[maxrow, 2]
	newpeak <- c(newx, 5, 5, newy)
	df <- data.frame(x=data$x, y=yval)
	#df1 <- subset(df, x > newx - 30 & x < newx + 30)
	df1 <- df
	newpeak <- voigtlearn2(df, newpeak, 10000, modx=FALSE, modh=FALSE)
	if (is.null(peaks))
		return(matrix(newpeak, nrow=1, ncol=4))
	return(rbind(peaks, newpeak))
}

table = read.table('table.txt', header = TRUE)
#table <- subset(table, x >= 1300 & x <= 1400)

#approx <- NULL
#approx <- extractpeak(table, approx)
#table$y <- table$y - voigt3(table$x, approx)
#print(approx)
#print(voigterror(table, approx))

#for (i in 1:100) {
	#approx <- voigtlearn2(table, approx, 10000, modx=TRUE, modh=TRUE)
	#print(approx)
	#print(voigterror(table, approx))
#}

#for (i in 1:20) {
	#approx <- extractpeak(table, approx)
	#table$y <- table$y - voigt3(table$x, approx)
	#print(approx)
	#print(voigterror(table, approx))
#}

approx <- rbind(c(1362.13,6.263586,6.588041,71.761264),
		c(1510.40,6.370530,6.662593,68.944266),
		c( 612.20,5.567818,5.563286,66.666902),
		c(1311.62,6.045249,6.125649,62.198678),
		c( 774.54,5.354478,5.401860,37.188178),
		c(1183.85,5.276569,5.349737,27.016442),
		c( 401.84,5.089704,5.111561,25.262224),
		c( 570.80,5.036813,4.947490,21.730405),
		c(1648.68,5.081232,5.079127,17.814100),
		c(1199.29,4.945842,4.806231,15.922607),
		c( 453.54,5.042014,5.040382,12.180987),
		c(1125.49,5.018237,5.003095, 9.441611),
		c( 659.73,4.976499,4.921704, 6.473517),
		c(1088.16,5.003444,4.992668, 6.198420),
		c( 526.86,4.998828,4.974003, 5.809724),
		c(1604.83,4.997256,4.968168, 5.296209),
		c(1012.62,5.000064,4.999955, 5.191784),
		c(1574.76,4.985285,4.942382, 3.899192),
		c(1044.57,4.999984,4.999897, 2.876482),
		c( 855.87,4.999994,4.999964, 2.465833),
		c( 828.90,4.997510,4.990354, 1.812561))


for (i in 1:1000) {
	approx <- voigtlearn2(table, approx, 1000, modx=FALSE, modh=FALSE, clearn=1e-8)
	print(approx)
	print(voigterror(table, approx))
}
