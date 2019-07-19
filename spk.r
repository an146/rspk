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
	ret <- 0
	for (i in 1:dim(peaks)[1]) {
		ret = ret + voigt2(xval, peaks[i,])
	}
	return(ret)
}

accumerror <- function(errors) {
	ret <- 0
	for (e in errors)
		ret <- ret + e * e
	return(ret)
}

voigtlearn <- function(xval, peaks) {
	deriv_coeff = 1e-5
	learn_coeff = 1e-3
	nextpeaks <- peaks
	ref <- voigt3(xval, peaks)
	for (i in 1:dim(peaks)[1]) {
		for (j in 1:dim(peaks)[2]) {
			peaks1 <- peaks
			peaks1[i, j] <- peaks1[i, j] + deriv_coeff
			nextpeaks[i, j] <- nextpeaks[i, j] - accumerror((voigt3(xval, peaks1) - ref) / deriv_coeff) * learn_coeff
		}
	}
	return(nextpeaks)
}

voigtlearn2 <- function(xval, peaks, count=1) {
	for (i in 1:count) {
		peaks <- voigtlearn(xval, peaks)
	}
	return(peaks)
}

voigterror <- function(xval, peaks, yval) {
	v <- voigt3(xval, peaks)
	return(accumerror(v - yval))
}

approx = rbind(c(400, 30, 60, 40), c(900, 60, 40, 60), c(1080, 60, 60, 40))
table = read.table('table.txt', header = TRUE)

print(voigt3(table$x, approx))
print(voigterror(table$x, approx, table$y))
print(voigterror(table$x, voigtlearn2(table$x, approx, 100), table$y))
