# plotDiversityEstimate.R

# Clear workspace:
rm(list=ls())

########################
### Define functions ###
########################

g <- function(h, L) {
	return (3^h*choose(L,h))
}

getDivAtTime <- function(df, tidx) {

	div <- list()

	L <- (length(df)-2)/2 - 1 # trust me

	NY <- 0
	nY2 <- 0
	NV <- 0
	nV2 <- 0

	for (h in seq(0,L)) {

		NY <- NY + df[[3 + h*2]][tidx]
		nY2 <- nY2 + df[[3 + h*2]][tidx]^2/g(h, L)

		NV <- NV + df[[4 + h*2]][tidx]
		nV2 <- nV2 + df[[4 + h*2]][tidx]^2/g(h, L)

	}

	div$y <- NY^2/nY2
	div$v <- NV^2/nV2

	return(div)
}

getDiv <- function (df) {

	div <- list()
	div$y <- vector(mode='numeric', length=length(df$t))
	div$v <- vector(mode='numeric', length=length(df$t))

	for (tidx in seq(1,length(df$t))) {

		thisdiv <- getDivAtTime(df, tidx)
		div$y[tidx] <- thisdiv$y
		div$v[tidx] <- thisdiv$v

	}

	div$t <- df$t

	return(div)
}

########################


# Load data from file:
df <- read.table('div.txt', header=T)

# Calculate diversity:
div <- getDiv(df)

########################
###   Plot figure    ###
########################

#pdf('div_deterministic.pdf', onefile=F, width=7, height=6)

plot(div$t, div$y, 'l', col='blue', xlab='Time (days)', ylab='Inverse Simpson Index', main='Deterministic diversity dynamics')
lines(div$t, div$v, col='red')
legend('bottomright', inset=.05, c('Infected cells', 'Virions'), lty=1, col=c('blue','red'))

#dev.off()

########################
