# plotDiversityEstimate.R

# Clear workspace:
rm(list=ls())

# Load data from file:
df <- read.table('div.txt', header=T)

# Functions:

getDivAtTime <- function(df, tidx) {

	res <- list()

	L <- (length(df)-2)/2 - 1 # trust me

	NY <- 0
	nY2 <- 0
	NV <- 0
	nV2 <- 0

	for (h in seq(0,L)) {

		NY <- NY + df[[3 + h*2]][tidx]
		nY2 <- nY2 + df[[3 + h*2]][tidx]^2

		NV <- NV + df[[4 + h*2]][tidx]
		nV2 <- nV2 + df[[4 + h*2]][tidx]^2

	}

	res$ydiv <- NY^2/nY2
	res$vdiv <- NV^2/nV2

	return(res)
}

getDiv <- function (df) {

	res <- list()
	res$ydiv <- vector(mode='numeric', length=length(df$t))
	res$vdiv <- vector(mode='numeric', length=length(df$t))

	for (tidx in seq(1,length(df$t))) {

		thisres <- getDivAtTime(df, tidx)
		res$ydiv[tidx] <- thisres$ydiv
		res$vdiv[tidx] <- thisres$vdiv

	}

	return(res)
}
