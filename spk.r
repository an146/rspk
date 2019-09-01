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
		if (dim(peaks)[1] > 0) {
			for (i in 1:dim(peaks)[1]) {
				ret <- ret + voigt2(xval, peaks[i,])
			}
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
	deriv_coeff = 1e-12
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

voigtlearn2 <- function(data, peaks, count=1, modx=FALSE, modh=FALSE, clearn=1e-5) {
	for (i in 1:count) {
		peaks <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
	}
	return(peaks)
}

voigtlearn3 <- function(data, peaks, modx=FALSE, modh=FALSE, clearn=1e-5) {
	ve <- voigterror(data, peaks)
	ve1 <- ve
	peaks1 <- peaks
	while (ve <= ve1) {
		peaks <- peaks1
		ve <- ve1
		peaks1 <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
		ve1 <- voigterror(data, peaks)
	}
	return(peaks)
}

extractpeak <- function(data, peaks, left) {
	maxrow <- 0
	maxvalue <- -Inf
	for (i in 1:length(left)) {
		if (maxvalue < left[i]) {
			maxvalue <- left[i]
			maxrow <- i
		}
	}
	newx <- data$x[maxrow]
	newy <- maxvalue
	newpeak <- c(newx, 5, 5, newy)
	df <- data.frame(x=data$x, y=left)
	df <- subset(df, x > newx - 3 & x < newx + 3)
	newpeak <- voigtlearn3(df, newpeak, modx=TRUE, modh=TRUE, clearn=1e-6)
	return(rbind(peaks, newpeak))
}

init <- function() {
	#Reading data to fit
	#
	table <- read.table('table.txt', header = TRUE)
	#table <- subset(table, x >= 1300 & x <= 1400)
	assign("table", table, envir = .GlobalEnv)

	#Reading saved fit
	#
	if (file.exists('fit.txt')) {
		fit <- read.table('fit.txt', header=TRUE, col.names=c('xpeak','fG','fL','height'))
	} else {
		fit <- data.frame(xpeak=numeric(),fG=numeric(),fL=numeric(),height=numeric())
	}
	approx <- cbind(fit$xpeak, fit$fG, fit$fL, fit$height)
	print(approx)
	left <- table$y - voigt3(table$x, approx)
	ve <- voigterror(table, approx)
	#assign("fit", fit, envir = .GlobalEnv)
	assign("approx", approx, envir = .GlobalEnv)
	assign("left", left, envir = .GlobalEnv)
	assign("ve", ve, envir = .GlobalEnv)
}

save_approx <- function() {
	t = data.frame(xpeak=approx[,1],fG=approx[,2],fL=approx[,3],height=approx[,4])
	write.table(t,'fit.txt')
}

do_extract_peak <- function() {
	approx <- extractpeak(table, approx, left)
	ve1 <- ve
	ve <- voigterror(table, approx)
	if (ve < ve1) {
		v <- voigt2(table$x, approx[dim(approx)[1],])
		left <- left - v
		print(approx)
		print(ve)
		assign("approx", approx, envir = .GlobalEnv)
		assign("left", left, envir = .GlobalEnv)
		assign("ve", ve, envir = .GlobalEnv)
		save_approx()
		return(TRUE)
	} else {
		return(FALSE)
	}
}

do_extract_peaks <- function(wanted_peaks) {
	while (do_extract_peak())
		print("peak extracted");
}

do_learn <- function(count=1, clearn=1e-5) {
	approx <- voigtlearn2(table, approx, count, modx=TRUE, modh=TRUE, clearn=clearn)
	left <- table$y - voigt3(table$x, approx)
	ve <- voigterror(table, approx)
	print(approx)
	print(ve)
	assign("approx", approx, envir = .GlobalEnv)
	assign("left", left, envir = .GlobalEnv)
	assign("ve", ve, envir = .GlobalEnv)
	save_approx()
}

do_learn_till_error <- function(err, clearn) {
	while (ve > err) {
		do_learn(1, clearn=clearn)
		do_extract_peak()
	}
}

do_plot <- function() {
	#Plot
	#
	v <- voigt3(table$x, approx)
	vt <- data.frame(x=table$x, y=v)
	plot(table$x, table$y, type="l", col="red")
	lines(table$x, vt$y, col="blue")
}

init()
do_extract_peaks()
do_learn_till_error(1400, 1e-5)
#do_extract_peaks(wanted_peaks=120)
#do_learn(10, 1e-4)
do_plot()
