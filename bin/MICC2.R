##### Rscript
require(VGAM)
args<-commandArgs(TRUE)

Default.Max <- 500
#===================================================
# 0. Truncated Zeta function
ZetaTruncated <- function(theta, max=Default.Max, Truncated=FALSE){
	if( Truncated ) {
		x <- c(1:max)
		theta.len <- length(theta)
		if(theta.len>1) {
			i <- 1
			v <- rep( 0, theta.len )
			while( i<=max ) {
				v <- v + x[i]^-theta
				i <- i + 1
			}
		}
		else {
			v <- sum( 1/x^theta )
		}
	}
	else {
		v <- zeta(theta)
	}
	v
}
D1.ZetaTruncated <- function(theta, max=Default.Max, Truncated=FALSE){
	if( Truncated ) {
		x <- c(1:max)
		theta.len <- length(theta)
		if(theta.len>1) {
			i <- 1
			v <- rep( 0, theta.len )
			while( i<=max ) {
				v <- v - log(x[i]) * x[i]^-theta
				i <- i + 1
			}
		}
		else {
			v <- sum( -log(x)/x^theta )
		}
	}
	else {
		v <- zeta(theta, deriv=1)
	}
	v
}

#===================================================
# 1. true interaction
F1 <- function(cAB, theta=theta){
	v <- 1/ZetaTruncated(theta) * cAB^-theta
	v
}
logF1 <- function(cAB, theta=theta) {
	v <- -log(ZetaTruncated(theta)) - theta*log(cAB)
	v
}
Diff.logF1.theta <- function(theta, cAB=cAB) {
	-D1.ZetaTruncated(theta) / ZetaTruncated(theta) - log(cAB)
}

Theta1 <- function(distance, params=params) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	distance.finite <- distance < Inf
	distance.len <- length(distance)
	v <- rep( par1, distance.len )
	#v[ distance.finite ] <- (par1*par2*distance[distance.finite] + 1) / (par2*distance[distance.finite] + 1)
	par.fdconstance <- 1000
	pd <- -par2*distance[distance.finite]
	v[distance.finite] <- (par1 * distance[distance.finite] + par3*par2*par.fdconstance) / ( par2*par.fdconstance + distance[distance.finite] ) + par4 / (10*distance[distance.finite])
	v
}
D1.Theta1.params <- function(params, distance=distance) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	distance.finite <- distance < Inf
	distance.len <- length(distance)
	par.fdconstance <- 1000
	v <- rep( c(1, 0, 0, 0), distance.len )
	v <- matrix( v, 4, distance.len )
	v <- t(v)
	pd <- ( par2*par.fdconstance + distance[distance.finite] )
	v[distance.finite, ] <- c( distance[distance.finite] / pd,
		(par3-par1)*distance[distance.finite]*par.fdconstance / pd^2,
		par2*par.fdconstance / pd,
		1 / (10*distance[distance.finite])
		)
	v
}

rF1 <- function(n, theta=theta, max=max, Truncated=FALSE){
	u0 <- runif(n, 0, 1)
	x <- 1
	cumF <- rep(0, n)
	r0 <- rep(0, n)
	Finished <- rep(TRUE, n)
	UnFinished <- rep(TRUE, n)
	UnFinished.label <- 1
	zeta.theta <- ZetaTruncated(theta, max=max, Truncated=Truncated)
	while( UnFinished.label ) {
		pF <- 1/zeta.theta[UnFinished] * x^-theta[UnFinished]
		cumF.unfinished <- cumF[UnFinished]
		u0.unfinished <- u0[UnFinished]
		r0.unfinished <- r0[UnFinished]
		UnFinished.unfinished <- UnFinished[UnFinished]
		cumF.unfinished <- cumF.unfinished + pF
		cumF[UnFinished] <- cumF.unfinished
		CompareF <- u0.unfinished < cumF.unfinished
		r0.unfinished[CompareF] <- x
		r0[UnFinished] <- r0.unfinished
		UnFinished.unfinished[CompareF] <- FALSE
		UnFinished[UnFinished] <- UnFinished.unfinished
		x <- x+1
		UnFinished.label <- sum(!CompareF)
	}
	r0
}

#===================================================
# 2. random collision
F2 <- function(cAB, theta=theta) {
	v <- 1/ZetaTruncated(theta) * cAB^-theta
	v
}
logF2 <- function(cAB, theta=theta) {
	v <- logF1(cAB, theta)
	v
}
Diff.logF2.theta <- function(theta, cAB=cAB) {
	-D1.ZetaTruncated(theta) / ZetaTruncated(theta) - log(cAB)
}

Theta2 <- function(distance, params=params) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	theta0 <- params[5]
	pars <- c(par1,par2,par3,par4)
	v <- Theta1(distance, pars)
	v <- v + theta0
	v
}
D1.Theta2.params <- function(params, distance=distance) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	theta0 <- params[5]
	pars <- c(par1,par2,par3,par4)
	v <- D1.Theta1.params( pars, distance )
	v0 <- 1
	v <- cbind(v, v0)
	v
}

rF2 <- function(n, theta=theta, max=max, Truncated=FALSE){
	u0 <- runif(n, 0, 1)
	x <- 1
	cumF <- rep(0, n)
	r0 <- rep(0, n)
	Finished <- rep(TRUE, n)
	UnFinished <- rep(TRUE, n)
	UnFinished.label <- 1
	zeta.theta <- ZetaTruncated(theta, max=max, Truncated=Truncated)
	while( UnFinished.label ) {
		pF <- 1/zeta.theta[UnFinished] * x^-theta[UnFinished]
		cumF.unfinished <- cumF[UnFinished]
		u0.unfinished <- u0[UnFinished]
		r0.unfinished <- r0[UnFinished]
		UnFinished.unfinished <- UnFinished[UnFinished]
		cumF.unfinished <- cumF.unfinished + pF
		cumF[UnFinished] <- cumF.unfinished
		CompareF <- u0.unfinished < cumF.unfinished
		r0.unfinished[CompareF] <- x
		r0[UnFinished] <- r0.unfinished
		UnFinished.unfinished[CompareF] <- FALSE
		UnFinished[UnFinished] <- UnFinished.unfinished
		x <- x+1
		UnFinished.label <- sum(!CompareF)
	}
	r0
}

#===================================================
# 3. random ligation
F3 <- function(cAB, cA, cB, N) {
	p1 <- dhyper( cAB, cA, 2*N-cA, cB )
	p0 <- dhyper( 0, cA, 2*N-cA, cB )
	p1/(1-p0)
}
logF3 <- function(cAB, cA, cB, N) {
	logp1 <- dhyper( cAB, cA, 2*N-cA, cB, log=TRUE )
	p0 <- dhyper( 0, cA, 2*N-cA, cB )
	v <- logp1 - log(1-p0)
	v
}

rF3 <- function(n, cA=cA, cB=cB, N=N){
	u0 <- runif(n, 0, 1)
	p0 <- dhyper( 0, cA, 2*N-cA, cB )
	x <- 1
	cumF <- rep(0, n)
	r0 <- rep(0, n)
	Finished <- rep(TRUE, n)
	UnFinished <- rep(TRUE, n)
	UnFinished.label <- 1
	while( UnFinished.label ) {
		pF <- dhyper( x, cA, 2*N-cA, cB, log=TRUE )
		pF <- exp(pF-log(1-p0))
		cumF.unfinished <- cumF[UnFinished]
		u0.unfinished <- u0[UnFinished]
		r0.unfinished <- r0[UnFinished]
		UnFinished.unfinished <- UnFinished[UnFinished]
		cumF.unfinished <- cumF.unfinished + pF[UnFinished]
		cumF[UnFinished] <- cumF.unfinished
		CompareF <- u0.unfinished < cumF.unfinished
		r0.unfinished[CompareF] <- x
		r0[UnFinished] <- r0.unfinished
		UnFinished.unfinished[CompareF] <- FALSE
		UnFinished[UnFinished] <- UnFinished.unfinished
		x <- x+1
		UnFinished.label <- sum(!CompareF)
	}
	r0
}

#===================================================
# 4. prior
PriorDistance <- function(distance, params=params) {
	ld <- log(distance)
	par1 <- params[1]
	par2 <- params[2]
	lambda0 <- params[3]
	distance.finite <- distance < Inf
	distance.len <- length(distance)
	lambda <- rep( lambda0, distance.len )
	v <- exp( par1 * ld[distance.finite] + par2 )
	lambda[ distance.finite ] <- lambda0 * v / ( 1 + v )
	lambda
}
LogPriorDistance <- function(distance, params=params) {
	ld <- log(distance)
	par1 <- params[1]
	par2 <- params[2]
	lambda0 <- params[3]
	distance.finite <- distance < Inf
	distance.len <- length(distance)
	lambda <- rep( log(lambda0), distance.len )
	v <- par1 * ld[distance.finite] + par2
	lambda[ distance.finite ] <- log(lambda0) + v - log( 1 + exp(v) )
	lambda
}
D1.PriorDistance.params <- function(params, distance=distance) {
	ld <- log(distance)
	par1 <- params[1]
	par2 <- params[2]
	lambda0 <- params[3]
	distance.finite <- distance < Inf
	distance.len <- length(distance)
	lambda <- rep( c(0, 0, 1), distance.len )
	lambda <- matrix( lambda, 3, distance.len )
	lambda <- t(lambda)
	v <- exp( par1 * ld[distance.finite] + par2 )
	lambda[ distance.finite, ] <- c( lambda0 * ld[distance.finite] * v / ( 1 + v )^2,
					lambda0 * v / ( 1 + v )^2,
					v / ( 1 + v ) )
	lambda
}

PriorcAcB <- function(cA, cB, params=params, Minus=FALSE) {
	lcA <- log(cA)
	lcB <- log(cB)
	par1 <- params[1]
	par2 <- params[2]
	v1 <- exp( par1 * lcA + par2 )
	v2 <- exp( par1 * lcB + par2 )
	if(!Minus) {
		mu <- 1 - v1 * v2 / ((1+v1)*(1+v2))
	}
	else {
		mu <- v1 * v2 / ((1+v1)*(1+v2))
	}
	mu[ mu<1e-300] <- 1e-300
	mu
}
LogPriorcAcB <- function(cA, cB, params=params) {
	lcA <- log(cA)
	lcB <- log(cB)
	par1 <- params[1]
	par2 <- params[2]
	v1 <- exp( par1 * lcA + par2 )
	v2 <- exp( par1 * lcB + par2 )
	mu <- log( 1 + v1 + v2 ) - log( 1 + v1 ) - log( 1 + v2 )
	#mu[mu==0] <- -1e-16
}

D1.PriorcAcB.params <- function(params, cA=cA, cB=cB) {
	lcA <- log(cA)
	lcB <- log(cB)
	par1 <- params[1]
	par2 <- params[2]
	v1 <- exp( par1 * lcA + par2 )
	v2 <- exp( par1 * lcB + par2 )
	D1.v1.par1 <- lcA * v1
	D1.v1.par2 <- v1
	D1.v2.par1 <- lcB * v2
	D1.v2.par2 <- v2
	logis.v1 <- v1/(1+v1)
	logis.v2 <- v2/(1+v2)
	D1.logis.v1.par1 <- D1.v1.par1 / (1+v1)^2
	D1.logis.v1.par2 <- D1.v1.par2 / (1+v1)^2
	D1.logis.v2.par1 <- D1.v2.par1 / (1+v2)^2
	D1.logis.v2.par2 <- D1.v2.par2 / (1+v2)^2
	cA.len <- length(cA)
	mu <- rep( c(0, 0), cA.len )
	mu <- matrix( mu, 2, cA.len )
	mu <- t(mu)
	mu[,] <- c( - (D1.logis.v1.par1 * logis.v2 + D1.logis.v2.par1 * logis.v1),
				- (D1.logis.v1.par2 * logis.v2 + D1.logis.v2.par2 * logis.v1)
				)
	mu
}

#===================================================
# 5. first order differentiation
Loglik.theta12.params <- function(params, cAB=cAB, label=label, distance=distance) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	theta0 <- params[5]
	pars <- c(par1,par2,par3,par4)
	theta1 <- Theta1( distance, pars )
	theta2 <- theta1 + theta0
	#
	D <- logF1( cAB, theta1) * label[,1] + logF2( cAB, theta2 ) * label[,2]
	-sum(D)
}
D1.Loglik.theta12.params <- function(params, cAB=cAB, label=label, distance=distance) {
	par1 <- params[1]
	par2 <- params[2]
	par3 <- params[3]
	par4 <- params[4]
	theta0 <- params[5]
	pars <- c(par1,par2,par3,par4)
	#
	Dif1.Theta1.params <- D1.Theta1.params( pars, distance )
	Dif1.Theta2.params <- D1.Theta2.params( params, distance )
	theta1 <- Theta1( distance, pars )
	theta2 <- theta1 + theta0
	Diff.logF1.theta.v <- Diff.logF1.theta( theta1, cAB )
	Diff.logF2.theta.v <- Diff.logF2.theta( theta2, cAB )
	#
	D <- Diff.logF1.theta.v * Dif1.Theta1.params[,1] * label[,1] + Diff.logF2.theta.v * Dif1.Theta2.params[,1] * label[,2]
	Dif1 <- sum(D)
	#
	D <- Diff.logF1.theta.v * Dif1.Theta1.params[,2] * label[,1] + Diff.logF2.theta.v * Dif1.Theta2.params[,2] * label[,2]
	Dif2 <- sum(D)
	#
	D <- Diff.logF1.theta.v * Dif1.Theta1.params[,3] * label[,1] + Diff.logF2.theta.v * Dif1.Theta2.params[,3] * label[,2]
	Dif3 <- sum(D)
	#
	D <- Diff.logF1.theta.v * Dif1.Theta1.params[,4] * label[,1] + Diff.logF2.theta.v * Dif1.Theta2.params[,4] * label[,2]
	Dif4 <- sum(D)
	#
	D <- Diff.logF2.theta.v * Dif1.Theta2.params[,5] * label[,2]
	Dif5 <- sum(D)
	-c(Dif1, Dif2, Dif3, Dif4, Dif5)
}


Loglik.PriorDistance.params <- function(params, label=label, distance=distance) {
	lambda <- PriorDistance( distance, params )
	D <- label[,1] * log(1-lambda)
	D <- D + label[,2] * log(1-lambda)
	D <- D + label[,3] * log(lambda)
	-sum(D)
}
D1.Loglik.PriorDistance.params <- function(params, label=label, distance=distance) {
	Dif1.PriorDistance.params <- D1.PriorDistance.params( params, distance )
	lambda <- PriorDistance( distance, params )
	#
	D <- -label[,1] / (1-lambda)
	D <- D - label[,2] / (1-lambda)
	D <- D + label[,3] / lambda
	Dif1 <- sum(D * Dif1.PriorDistance.params[,1])
	#
	D <- -label[,1] / (1-lambda)
	D <- D - label[,2] / (1-lambda)
	D <- D + label[,3] / lambda
	Dif2 <- sum(D * Dif1.PriorDistance.params[,2])
	#
	D <- -label[,1] / (1-lambda)
	D <- D - label[,2] / (1-lambda)
	D <- D + label[,3] / lambda
	Dif3 <- sum(D * Dif1.PriorDistance.params[,3])
	-c( Dif1, Dif2, Dif3 )
}

Loglik.PriorcAcB.params <- function(params, label=label, cA=cA, cB=cB) {
	lambda <- PriorcAcB( cA, cB, params)
	mlambda <- PriorcAcB( cA, cB, params, Minus=TRUE)
	D <- label[,1] * log(mlambda)
	D <- D + label[,2] * log(lambda)
	-sum(D)
}
D1.Loglik.PriorcAcB.params <- function(params, label=label, cA=cA, cB=cB) {
	Dif1.PriorcAcB.params <- D1.PriorcAcB.params( params, cA, cB )
	lambda <- PriorcAcB( cA, cB, params )
	mlambda <- PriorcAcB( cA, cB, params, Minus=TRUE)
	#
	D <- -label[,1] / mlambda
	D <- D + label[,2] / lambda
	Dif1 <- sum(D * Dif1.PriorcAcB.params[,1])
	#
	D <- -label[,1] / mlambda
	D <- D + label[,2] / lambda
	Dif2 <- sum(D * Dif1.PriorcAcB.params[,2])
	-c( Dif1, Dif2 )
}

#===================================================
# 6. Maxmize logliklihood and get best parameters

Solve.Theta12.params <- function(cAB=cAB, label=label, distance=distance, InitialValue=InitialValue, lower=c(2,0,1,0,0), upper=c(20,20,20,20,20)) {
	par <- optim(InitialValue, fn=Loglik.theta12.params, gr=D1.Loglik.theta12.params, cAB=cAB, label=label, distance=distance, lower=lower, upper=upper, method="L-BFGS-B")
	par
}

Solve.PriorDistance.params <- function(label=label, distance=distance, InitialValue=InitialValue, lower=c(0,-50,0.01), upper=c(20,0,0.9999)) {
	par <- optim(InitialValue, fn=Loglik.PriorDistance.params, gr=D1.Loglik.PriorDistance.params, label=label, distance=distance, lower=lower, upper=upper, method="L-BFGS-B")
	par
}

Solve.PriorcAcB.params <- function(label=label, cA=cA, cB=cB, InitialValue=InitialValue, lower=c(0,-40), upper=c(20,0)) {
	par <- optim(InitialValue, fn=Loglik.PriorcAcB.params, gr=D1.Loglik.PriorcAcB.params, label=label, cA=cA, cB=cB, lower=lower, upper=upper, method="L-BFGS-B")
	par
}

#===================================================
# 7. Estimate Initial values

Estimate.Theta12.params <- function(cAB=cAB, distance=distance, InitialValue=c(2,1,1,1,1), lower=c(1.001,0,1,0,0), upper=c(10,10,10,10,10), MaxConfident=5) {
	dataConfident <- (distance<1E6)
	data.cAB <- cAB[dataConfident]
	data.distance <- distance[dataConfident]
	data.label <- matrix(0, length(data.cAB), 3)
	confidentN = sum(data.cAB>MaxConfident)  #by Guipeng Li
	if ( confidentN<100 ){
		stop("Don't have enough confident interactions to learn the model.", call.=FALSE)
	}
	data.label[data.cAB > MaxConfident,1] <- 1
	data.label[data.cAB <= MaxConfident,2] <- 1
	par <- Solve.Theta12.params(cAB=data.cAB, label=data.label, distance=data.distance, InitialValue=InitialValue, lower=lower, upper=upper)
}

Estimate.PriorDistance.params <- function(cAB=cAB, distance=distance, InitialValue=c(3,3,0.5), lower=c(0,-20,0.01), upper=c(20,0,0.999)) {
	label <- matrix(0, length(cAB), 3)
	label[cAB>5, 1] <- 1
	label[cAB<5 & distance<1E3, 2] <- 1
	label[cAB<5 & distance>1E3, 3] <- 1
	par1 <- runif(1,0.9,1.1)
	par2 <- -par1 * log(1000)
	lambda0 <- sum(label[,3])/length(cAB)
	par <- c( par1, par2, lambda0 )
	par
}

Estimate.PriorcAcB.params <- function(cAB=cAB, distance=distance, cA=cA, cB=cB, InitialValue=c(1,-4,0.9), lower=c(0,-20), upper=c(20,0)) {
	label <- matrix(0, length(cAB), 3)
	label[cAB>5, 1] <- 1
	label[cAB<5 & distance<1E3, 2] <- 1
	label[cAB<5 & distance>1E3, 3] <- 1
	par1 <- runif(1,3,4)
	tmp <- quantile( c( log(cA[label[,1]==1]), log(cB[label[,1]==1]) ), 0.1)
	par2 <- -par1 * tmp
	par <- c( par1, par2 )
	par
}

#===================================================
# 7. EM Algorithm
# log likelihoood
CompAllLogProb <- function( cAB, distance, cA, cB, params ) {
	N <- sum(cAB)
	pars <- params$theta12.par[1:4]
	theta0 <- params$theta12.par[5]
	theta1 <- Theta1(distance, pars)
	theta2 <- theta1 + theta0
	lambda <- params$lambda.par
	mu <- params$mu.par
	LogF1.v <- logF1(cAB, theta1)
	LogF2.v <- logF2(cAB, theta2)
	LogF3.v <- logF3(cAB, cA, cB, N)
	Loglambda.v <- LogPriorDistance(distance, lambda)
	Logmu.v <- LogPriorcAcB(cA, cB, mu)
	LogProb <- list( LogF1.v=LogF1.v, LogF2.v=LogF2.v, LogF3.v=LogF3.v, Loglambda.v=Loglambda.v, Logmu.v=Logmu.v )
}
CompLoglik <- function( LogProb ) {
	LogF1.v <- LogProb$LogF1.v
	LogF2.v <- LogProb$LogF2.v
	LogF3.v <- LogProb$LogF3.v
	Loglambda.v <- LogProb$Loglambda.v
	Logmu.v <- LogProb$Logmu.v
	#
	F1.v <- exp(LogProb$LogF1.v)
	F2.v <- exp(LogProb$LogF2.v)
	F3.v <- exp(LogProb$LogF3.v)
	lambda.v <- exp(LogProb$Loglambda.v)
	mu.v <- exp(LogProb$Logmu.v)
	#
	Logmmu.v <- log( 1 - mu.v )
	verysmall <- which(Logmmu.v== -Inf)
	Logmmu.v[verysmall] <- log( -Logmu.v[verysmall] ) + log( 1 + Logmu.v[verysmall]/2 )

	LogProb1 <- Logmmu.v + log( 1 - lambda.v ) + LogProb$LogF1.v
	LogProb2 <- Logmu.v + log( 1 - lambda.v ) + LogProb$LogF2.v
	LogProb3 <- LogProb$Loglambda.v + LogProb$LogF3.v
	#
	LogLik <- LogProb1
	LogLik <- LogLik + log( 1 + exp( LogProb2-LogProb1 ) + exp( LogProb3-LogProb1 ) )
	#
	PostProb1 <- exp(LogProb1 - LogLik)
	PostProb2 <- exp(LogProb2 - LogLik)
	PostProb3 <- exp(LogProb3 - LogLik)
	v <- list( LogLik=sum(LogLik), PostProb=cbind( PostProb1, PostProb2, PostProb3 ) )
	v
}

EMIter <- function( data, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5 ) {
	cAB <- data$cAB
	distance <- data$distance
	cA <- data$cA
	cB <- data$cB

	# parameter initialization
	if(is.null(params.init)) {
		theta12.par <- Estimate.Theta12.params(cAB, distance)
		lambda.par <- Estimate.PriorDistance.params(cAB, distance)
		mu.par <- Estimate.PriorcAcB.params(cAB, distance, cA, cB)
		params.init <- list( theta12.par=theta12.par$par, lambda.par=lambda.par, mu.par=mu.par )
	}

	CONVERGE <- 0
	RELCONVERGE <- 0
	num.restart <- 1
	while( num.restart < restart && CONVERGE==0 && RELCONVERGE==0 ) {
		num.step <- 1
		randshift <- runif(4, 0.8, 1)
		oldparams <- params.init
		oldparams$theta12.par <- params.init$theta12.par
		oldparams$lambda.par <- params.init$lambda.par
		oldparams$mu.par <- params.init$mu.par
		# E step
		LogProb <- CompAllLogProb( cAB, distance, cA, cB, oldparams )
		v <- CompLoglik( LogProb )
		oldLogLik <- v$LogLik
		label <- v$PostProb
		params <- oldparams
		while( num.step < step && CONVERGE==0 && RELCONVERGE==0 ) {
			# M step
			params$theta12.par <- Solve.Theta12.params(cAB, label, distance, oldparams$theta12.par)$par
			params$lambda.par <- Solve.PriorDistance.params(label, distance, oldparams$lambda.par)$par
			params$mu.par <- Solve.PriorcAcB.params(label, cA, cB, oldparams$mu.par)$par
			# E step
			LogProb <- CompAllLogProb( cAB, distance, cA, cB, params )
			v <- CompLoglik( LogProb )
			LogLik <- v$LogLik
			label <- v$PostProb
			absdiff <- LogLik - oldLogLik
			reldiff <- 1 - LogLik / oldLogLik
			if( absdiff < abstol ) {
				CONVERGE <- 1
			}
			else if ( reldiff < reltol ) {
				RELCONVERGE <- 1
			}
			else {
				oldparams <- params
				oldLogLik <- LogLik
			}
			## for bug fix
			cat( date(), "step=", num.step, " || absdiff=", absdiff, " || reldiff=", reldiff, "\n" )
			#cat( params$lambda.par, " || ", params$mu.par, "\n" )
			#cat( params$theta12.par, "\n" )
			num.step <- num.step + 1
		}
		num.restart <- num.restart + 1
	}
	if( CONVERGE==1 ) {
		cat( "MICC CONVERGED!\n" );
	}
	else if( RELCONVERGE==1 ) {
		cat( "MICC RELATIVELY CONVERGED!\n" );
	}
	else {
		cat( "MICC NOT CONVERGED!\n" );
	}
	results <- list( params=params, LogLik=LogLik, PostProb=label )
}

#===================================================
# estimate FDR
rMyMultinom <- function( n, p ) {
	u0 <- runif(n, 0, 1)
	x <- 1
	cumF <- rep(0, n)
	r0 <- rep(0, n)
	Finished <- rep(TRUE, n)
	UnFinished <- rep(TRUE, n)
	UnFinished.label <- 1
	while( UnFinished.label ) {
		pF <- p[UnFinished,x]
		cumF.unfinished <- cumF[UnFinished]
		u0.unfinished <- u0[UnFinished]
		r0.unfinished <- r0[UnFinished]
		UnFinished.unfinished <- UnFinished[UnFinished]
		cumF.unfinished <- cumF.unfinished + pF
		cumF[UnFinished] <- cumF.unfinished
		CompareF <- u0.unfinished < cumF.unfinished
		r0.unfinished[CompareF] <- x
		r0[UnFinished] <- r0.unfinished
		UnFinished.unfinished[CompareF] <- FALSE
		UnFinished[UnFinished] <- UnFinished.unfinished
		x <- x+1
		UnFinished.label <- sum(!CompareF)
	}
	r0
}

rMICC <- function( data, params ) {
	cAB <- data$cAB
	distance <- data$distance
	cA <- data$cA
	cB <- data$cB
	N <- sum(cAB)
	assign("Zipf.Max", max(cAB))
	num.sample <- length(cAB)
	LogProb <- CompAllLogProb( cAB, distance, cA, cB, params )
	LogF1.v <- LogProb$LogF1.v
	LogF2.v <- LogProb$LogF2.v
	LogF3.v <- LogProb$LogF3.v
	Loglambda.v <- LogProb$Loglambda.v
	Logmu.v <- LogProb$Logmu.v
	F1.v <- exp(LogProb$LogF1.v)
	F2.v <- exp(LogProb$LogF2.v)
	F3.v <- exp(LogProb$LogF3.v)
	lambda.v <- exp(LogProb$Loglambda.v)
	mu.v <- exp(LogProb$Logmu.v)
	# rand generating labels
	Prior <- cbind( ( 1 - mu.v )*( 1 - lambda.v ), mu.v *( 1 - lambda.v ), lambda.v )
	label <- rMyMultinom( num.sample, Prior )
	# rand generating cAB
	theta1.par <- params$theta12.par[1:4]
	theta1 <- Theta1( distance, theta1.par )
	num.label1 <- sum(label==1)
	r1 <- rF1(num.label1, theta=theta1[label==1], max=Zipf.Max, Truncated=TRUE)
	#
	theta2.par <- params$theta12.par
	theta2 <- Theta2(distance, theta2.par)
	num.label2 <- sum(label==2)
	r2 <- rF2(num.label2, theta=theta2[label==2], max=Zipf.Max, Truncated=TRUE)
	#
	num.label3 <- sum(label==3)
	r3 <- rF3(num.label3, cA[label==3], cB=cB[label==3], N=N)
	#
	r <- rep(NA, num.sample)
	r[label==1] <- r1
	r[label==2] <- r2
	r[label==3] <- r3
	list( r=r, label=label )
}

FDRestimate <- function( data, params, PostProb ) {
	rnd <- rMICC( data, params )
	cA <- data$cA
	cB <- data$cB
	cAB <- rnd$r

	label <- rnd$label
	distance <- data$distance
	LogProb <- CompAllLogProb( cAB, distance, cA, cB, params )
	v <- CompLoglik( LogProb )
	rndPostProb <- v$PostProb[,1]
	#
	PostProb <- 1-PostProb
	rndPostProb <- 1-rndPostProb
	PostProb.sort <- sort(PostProb)
	breaks <- unique(PostProb.sort)
	t1 <- min(rndPostProb)-1e-5
	if( min(breaks)>min(rndPostProb) ) {
		breaks <- c(t1, breaks)
	}
	t1 <- max(PostProb,rndPostProb)+0.01
	breaks <- c(breaks,t1)
	PostProb.breaks <- hist(PostProb, breaks, plot=FALSE)$counts
	rndPostProb.true.breaks <- hist(rndPostProb[label==1], breaks, plot=FALSE)$counts
	rndPostProb.false.breaks <- hist(rndPostProb[label!=1], breaks, plot=FALSE)$counts
	#
	rndPostProb.true.cum <- cumsum(rndPostProb.true.breaks)
	rndPostProb.false.cum <- cumsum(rndPostProb.false.breaks)
	FDR <- rndPostProb.false.cum / (rndPostProb.true.cum+0.01)
	names(FDR) <- 1-breaks[1:(length(breaks)-1)]
	PostProbtoname <- as.character(1-PostProb)
	PostProb.fdr <- FDR[PostProbtoname]
	PostProb.fdr[PostProb.fdr>1] <- 1
	PostProb.fdr
}


###########################################################################################
########### MICCMain.R
########### functions:
###########            InputMatrixFormatted: process input PET clusters into MICC main model input
###########            MICCMainLearn: fit model parameters
###########            FDRcompute: compute FDR
###########			   MICCoutput: implement MICC model and get output
###########################################################################################
####=============================================================================

InputMatrixFormatted <- function(data) {
	x <- data
	cAB <- x[,7]
	cA <- x[,8]
	cB <- x[,9]
	intra <- as.character(x[,1])==as.character(x[,4])
	inter <- as.character(x[,1])!=as.character(x[,4])
	distance <- (x[,5]+x[,6]-x[,3]-x[,2]) / 2 / 1000
	distance[inter] <- Inf
	data_formatted <- list( cAB=cAB, cA=cA, cB=cB, distance=distance )
}


MICCMainLearn <- function(data_formatted, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5) {
	Par <- EMIter( data_formatted, params.init=params.init, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident )
	Par
}

FDRcompute <- function( data_formatted, params, PostProb ) {
	fdr <- FDRestimate( data_formatted, params, PostProb )
	fdr
}

MICCoutput <- function( data, outfilename, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5 ) {
	data_formatted <- InputMatrixFormatted(data)
	Par <- MICCMainLearn( data_formatted, params.init=params.init, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident )
	params <- Par$params
	PostProb <- Par$PostProb
        cat('Calculating FDR... \n')
	fdr <- FDRcompute( data_formatted, params, PostProb[,1] )
	output.colnames <- c("chr.", "start", "end", "chr.", "start", "end", "cAB", "cA", "cB", "-log10(1-PostProb)", "fdr")
	y <- cbind(data, -log10(1-PostProb[,1]), fdr)
	colnames(y) <- output.colnames
        cat('Writing to file... \n')
	write.table( y, file=outfilename, sep="\t", row.names=F, quote=F )
}

InputMatrixFormatted2 <- function(data, cutoff=2) { ## by Guipeng Li, for ChIA-PET2
	x <- subset(data,V11>=cutoff)
	cAB <- x[,11]
	cA <- x[,9]
	cB <- x[,10]
	intra <- as.character(x[,1])==as.character(x[,4])
	inter <- as.character(x[,1])!=as.character(x[,4])
        cat('Intra: ', sum(intra), "\tInter:", sum(inter),'\n' )
	distance <- (x[,5]+x[,6]-x[,3]-x[,2]) / 2 / 1000
	#distance[inter] <- 1E10
	distance[inter] <- Inf
	data_formatted <- list( cAB=cAB, cA=cA, cB=cB, distance=distance )
}

## by Guipeng Li, for ChIA-PET2
MICCoutput2 <- function( data, outfilename, params.init=NULL, reltol=1e-5, abstol=1e-3, step=200, restart=5, MinConfident=5, cutoff=2) {
	data_formatted <- InputMatrixFormatted2(data,cutoff)
	Par <- MICCMainLearn( data_formatted, params.init=params.init, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident )
	params <- Par$params
	PostProb <- Par$PostProb
        cat('Calculating FDR... \n')
	fdr <- FDRcompute( data_formatted, params, PostProb[,1] )
	output.colnames <- c("chr", "start", "end", "chr", "start", "end", "peakA","peakB","cA", "cB","cAB", "-log10(1-PostProb)", "fdr")
	y <- cbind(subset(data,V11>=cutoff), -log10(1-PostProb[,1]), fdr)
	colnames(y) <- output.colnames
        cat('Writing to file... \n')
	write.table( y, file=outfilename, sep="\t", row.names=F, quote=F )
}


## added by Guipeng Li
argslen=length(args)
inputbedpe1="input1"
inputbedpe2="input2"
outputbedpe="output.MICC"
reltol=1e-8
abstol=1e-4
step=100
restart=10
cutoff=2
MinConfident=5
if (argslen<3) {
  stop("At least 3 arguments must be supplied (input, output).", call.=FALSE)
}else if (argslen>=3){
	inputbedpe1<-args[1]
	inputbedpe2<-args[2]
	outputbedpe<-args[3]
}
if (argslen>=4){
	cutoff<-as.numeric(args[4])
}
if (argslen>=5){
	MinConfident<-as.numeric(args[5])
}
if (argslen>=6){
	reltol<-as.numeric(args[6])
}
if (argslen>=7){
	step<-as.numeric(args[7])
}
if (argslen>=8){
	restart<-as.numeric(args[8])
}
if (argslen==9){
	abstol<-as.numeric(args[9])
}else if (argslen>9){
  stop("Too many arguments supplied.", call.=FALSE)
}
cat("Running MICC...\n")
cat('Intra data: ', inputbedpe1,'\n')
cat('Inter data: ', inputbedpe2,'\n')
cat('Output file: ', outputbedpe,'\n')
cat('PET count cutoff: ', cutoff,'\n')
cat('Minimun Confident PET count: ', MinConfident,'\n')
cat('reltol: ', reltol,'\n')

cat("Loading intra data...\n")
dat1=read.table(inputbedpe1,head=F)
cat("Loading inter data...\n")
dat2=read.table(inputbedpe2,head=F)
cat("Cacluating...\n")
MICCoutput2( rbind(dat1,dat2), outputbedpe, params.init=NULL, reltol=reltol, abstol=abstol, step=step, restart=restart, MinConfident=MinConfident, cutoff=cutoff )
cat("Job finished\n")
## end
