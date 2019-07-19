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

voigtlearn <- function(data, peaks) {
	deriv_coeff = 1e-7
	learn_coeff = 1e-5
	if (is.null(dim(peaks)))
		peaks <- matrix(peaks, ncol=4, nrow=1)
	nextpeaks <- peaks
	xval <- data$x
	yval <- data$y
	eref <- voigterror(data, peaks)
	for (i in 1:dim(peaks)[1]) {
		for (j in 1:dim(peaks)[2]) {
			peaks1 <- peaks
			peaks1[i, j] <- peaks1[i, j] + deriv_coeff
			deriv <- (voigterror(data, peaks1) - eref) / deriv_coeff
			nextpeaks[i, j] <- nextpeaks[i, j] - deriv * learn_coeff
		}
	}
	return(nextpeaks)
}

voigtlearn2 <- function(data, peaks, count=1) {
	for (i in 1:count) {
		peaks <- voigtlearn(data, peaks)
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
	newpeak <- c(data[i, 1], 15, 15, data[i, 2])
	newpeak <- voigtlearn2(data.frame(x=data$x, y=yval), newpeak, 50000)
	if (is.null(peaks))
		return(matrix(newpeak, nrow=1, ncol=4))
	return(rbind(peaks, newpeak))
}

table = read.table('table.txt', header = TRUE)
approx <- NULL
approx <- extractpeak(table, approx)
print(approx)
print(voigterror(table, approx))
for (i in 1:20) {
	approx <- voigtlearn2(table, approx, 20000)
	print(approx)
	print(voigterror(table, approx))
}
for (i in 1:20) {
	approx <- extractpeak(table, approx)
	print(voigterror(table, approx))
	print(voigt3(table$x, approx))
	print(approx)
}

for (i in 1:1000) {
	approx <- voigtlearn2(table, approx, 10000)
	print(voigterror(table, approx))
}
