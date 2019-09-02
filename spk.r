#install.packages('RcppFaddeeva')
library('RcppFaddeeva')

Faddeeva_w_ <- function(z) {
	if (is.atomic(z) && length(z) == 1 && is.complex(z)) {
		return((-2*z)*Faddeeva_w(z) + exp(-z*z)*2i/sqrt(pi)*exp(z*z))
	} else {
		stopifnot(length(z) > 0)
		for (i in 1:length(z)) {
			z[i] <- (-2*z[i])*Faddeeva_w(z[i]) + exp(-z[i]*z[i])*2*complex(imaginary=1)/sqrt(pi)*exp(z[i]*z[i])
		}
		return(z)
	}
}

voigt <- function(xval, xpeak, sigma, ygamma, height) {
	ygamma <- Re(ygamma)
	z <- complex(real = xval - xpeak, imaginary = rep(ygamma, length(xval))) / (sqrt(2) * sigma)
	peakZ <- complex(imaginary = ygamma) / (sqrt(2) * sigma)
	return(height * Re(Faddeeva_w(z)) / Re(Faddeeva_w(peakZ)))
}

voigts <- function(xval, peaks) {
	ret <- numeric(length(xval))
	if (nrow(peaks) > 0) {
		for (i in 1:nrow(peaks)) {
			v <- voigt(xval, peaks[i, 'xpeak'], peaks[i, 'sigma'], peaks[i, 'ygamma'], peaks[i, 'height'])
			ret <- ret + v
		}
	}
	return(ret)
}

voigt_xpeak <- function(ref, peak) {
	z <- complex(real = ref$x - peak$xpeak, imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	z_ <- -1
	peakZ <- complex(imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	return(peak$height * Re(z_) * Re(Faddeeva_w(z)) / Re(Faddeeva_w(peakZ)))
}

voigt_sigma <- function(ref, peak) {
	z <- complex(real = ref$x - peak$xpeak, imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	z_ <- -z / peak$sigma
	peakZ <- complex(imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	peakZ_ <- -peakZ / peak$sigma
	return(Re(peak$height * (z_ * Re(Faddeeva_w_(z)) * Re(Faddeeva_w(peakZ)) - z * Re(Faddeeva_w_(z)) * peakZ_ * Re(Faddeeva_w_(peakZ))) / Re(Faddeeva_w(peakZ)) / Re(Faddeeva_w(peakZ))))
}

voigt_ygamma <- function(ref, peak) {
	z <- complex(real = ref$x - peak$xpeak, imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	z_ <- complex(imaginary = 1) / (sqrt(2) * peak$sigma)
	peakZ <- complex(imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	peakZ_ <- z_
	return(Re(peak$height * (z_ * Re(Faddeeva_w_(z)) * Re(Faddeeva_w(peakZ)) - Re(Faddeeva_w(z)) * peakZ_ * Re(Faddeeva_w_(peakZ)) / Re(Faddeeva_w(peakZ)) / Re(Faddeeva_w(peakZ)))))
}

voigt_height <- function(ref, peak) {
	z <- complex(real = ref$x - peak$xpeak, imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	peakZ <- complex(imaginary = peak$ygamma) / (sqrt(2) * peak$sigma)
	return(Re(Faddeeva_w(z)) / Re(Faddeeva_w(peakZ)))
}

accumerror <- function(errors) {
	ret <- 0
	for (e in errors)
		ret <- ret + e * e
	return(ret)
}

voigterror <- function(data, peaks) {
	v <- voigts(data$x, peaks)
	return(accumerror(data$y - v))
}

voigtlearn <- function(data, peaks, modx, modh, clearn) {
	learn_coeff = clearn
	nextpeaks <- peaks
	if (nrow(peaks) > 0) {
		for (i in 1:nrow(peaks)) {
			v <- voigts(data$x, peaks[i,])
			ref <- data.frame(x=data$x, y=data$y - voigts(data$x, peaks) + v)
			if (modx) {
				nextpeaks$xpeak[i] <- nextpeaks$xpeak[i] + sum(2 * (ref$y - v) * voigt_xpeak(ref, peaks[i,]) * learn_coeff)
			}
			v <- voigts(data$x, peaks[i,])
			nextpeaks$sigma[i] <- nextpeaks$sigma[i] + sum(2 * (ref$y - v) * voigt_sigma(ref, peaks[i,]) * learn_coeff)
			nextpeaks$ygamma[i] <- nextpeaks$ygamma[i] + sum(2 * (ref$y - v) * voigt_ygamma(ref, peaks[i,]) * learn_coeff)
			if (modh) {
				nextpeaks$height[i] <- nextpeaks$height[i] + sum(2 * (ref$y - v) * voigt_height(ref, peaks[i,]) * learn_coeff)
			}
			v <- voigts(data$x, peaks[i,])
		}
	}
	ve <- voigterror(data, peaks)
	ve1 <- voigterror(data, nextpeaks)
	stopifnot(!is.nan(ve))
	if (is.nan(ve1) || ve <= ve1)
		return(peaks)
	else
		return(nextpeaks)
}

voigtlearn2 <- function(data, peaks, count=1, modx=FALSE, modh=FALSE, clearn=1e-5) {
	for (i in 1:count) {
		peaks <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
	}
	return(peaks)
}

voigtlearn3 <- function(data, peaks, maxcount=1000000, modx=FALSE, modh=FALSE, clearn=1e-5) {
	ve <- voigterror(data, peaks)
	peaks1 <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
	ve1 <- voigterror(data, peaks1)
	n <- 0
	while (n < maxcount && ve1 < ve) {
		peaks <- peaks1
		ve <- ve1
		peaks1 <- voigtlearn(data, peaks, modx=modx, modh=modh, clearn=clearn)
		ve1 <- voigterror(data, peaks1)
		n <- n + 1
	}
	return(peaks)
}

extractpeak <- function(data, peaks, left) {
	maxrow <- 0
	maxvalue <- -Inf
	for (i in 1:nrow(left)) {
		if (maxvalue < left$y[i]) {
			maxvalue <- left$y[i]
			maxrow <- i
		}
	}
	newx <- data$x[maxrow]
	newy <- maxvalue
	fG <- 5
	fL <- 5
	sigma <- fG / (2 * (2*log(2))^0.5)
	ygamma <- fL / 2
	newpeak <- data.frame(xpeak=newx, sigma=sigma, ygamma=ygamma, height=newy)
	df <- subset(left, x > newx - 5 & x < newx + 5)
	newpeak <- voigtlearn3(df, newpeak, maxcount=1000, modx=TRUE, modh=TRUE, clearn=1e-5)
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
		fit <- read.table('fit.txt', header=TRUE, col.names=c('xpeak','sigma','ygamma','height'))
	} else {
		fit <- data.frame(xpeak=numeric(),sigma=numeric(),ygamma=numeric(),height=numeric())
	}
	print(fit)
	left <- data.frame(x=table$x, y=table$y - voigts(table$x, fit))
	ve <- voigterror(table, fit)
	assign("fit", fit, envir = .GlobalEnv)
	assign("left", left, envir = .GlobalEnv)
	assign("ve", ve, envir = .GlobalEnv)
}

save_fit <- function() {
	write.table(fit, 'fit.txt')
	print("saved")
}

do_extract_peak <- function() {
	fit <- extractpeak(table, fit, left)
	ve_prev <- ve
	ve <- voigterror(table, fit)
	if (ve < ve_prev) {
		v <- voigts(table$x, fit[nrow(fit)[1],])
		left <- left - v
		print(fit)
		print(ve)
		assign("fit", fit, envir = .GlobalEnv)
		assign("left", left, envir = .GlobalEnv)
		assign("ve", ve, envir = .GlobalEnv)
		save_fit()
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
	fit <- voigtlearn2(table, fit, count, modx=TRUE, modh=TRUE, clearn=clearn)
	left <- data.frame(x=table$x, y=table$y - voigts(table$x, fit))
	ve <- voigterror(table, fit)
	print(fit)
	print(ve)
	assign("fit", fit, envir = .GlobalEnv)
	assign("left", left, envir = .GlobalEnv)
	assign("ve", ve, envir = .GlobalEnv)
	save_fit()
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
#do_extract_peaks()
do_learn_till_error(1400, 1e-5)
#do_extract_peaks(wanted_peaks=120)
#do_learn(10, 1e-4)
do_plot()
