## Useage:
## res_vca <- vca(x_n, p, snr_input=1, SNR=SNR)
## Each column is a data point. 

library("MASS")

dirichlet_rnd<-function(a, dir_dim) {

#DIRICHLET_RND Random matrices from dirichlet distribution.
#   R = DIRICHLET_RND(A,DIM) returns a matrix of random numbers chosen   
#   from the dirichlet distribution with parameters vector A.
#   Size of R is (N x N) where N is the size of A or (N x DIM) if DIM is given.

rows <- dim(a)[1]
columns <- dim(a)[2]
if (missing(dir_dim)) {
	dir_dim <-rows * columns
} else {
	if (prod(dim(dir_dim)) != 1) {
		stop('The second parameter must be a scalar.')
	}
}
if (  (rows != 1 ) & (columns != 1)  ) {
	stop('Requires a vector as an argument.')
}
if (   any(a<0 | !array( sapply(a, is.double), dim=dim(a) ) )  ){
	stop('Parameters of Dirichlet function must be real positive values.')
}

# fastest method that generalize method used to sample from
# the BETA distribuition: Draw x1,...,xk from independent gamma 
# distribuitions with common scale and shape parameters a1,...,ak, 
# and for each j let rj=xj/(x1+...+xk).

N <- rows * columns
x <- array(NA, dim=c(dir_dim, N))
for (i in 1:N) {
	x[,i]<-array(rgamma(dir_dim, shape=a[i], scale=a[i]), dim=c(dir_dim,1))
}
r <- x/replicate(N, rowSums(x))
# For more details see "Bayesian Data Analysis" Appendix A
} #end of dirichlet_rnd

# Internal functions
estimate_snr <-function(R, r_m, x) {
	L <- dim(R)[1]
	N <- dim(R)[2]
	p <- dim(x)[1]
	N <- dim(x)[2]
	P_y <- sum(c(R)^2)/N
	P_x <- sum(c(x)^2)/N + t(r_m)%*%r_m
	snr_est <- 10*log10(  (P_x - p/L*P_y)/(P_y-P_x)  )
	return(snr_est)
} #end

vca<-function(R, p, verbose = 'on', snr_input = 0, SNR = 0) {

# Initializations
L <- dim(R)[1]
N <- dim(R)[2]

if ( (p<0) | (p>L) | ((p%%1)!=0) ) {
	stop('ENDMEMBER parameter must be integer between 1 and L')
}

# SNR Estimates

if (snr_input == 0) {
	r_m <-rowMeans(R)
	R_m <-replicate(N, r_m)
	R_o <-R - R_m
	res_svd <- svd(R_o%*%t(R_o)/N, nu=p, nv=p)
	Ud <- res_svd$u
	Sd <- res_svd$d
	Vd <- res_svd$v
	x_p <-t(Ud)%*%(R_o)

	SNR <- estimate_snr(R, r_m, x_p)

	if (verbose == 'on') {
		message(sprintf("SNR estimated = %g[dB]\n", SNR))
	}
} else {
	if (verbose == 'on') {
		message(sprintf("input    SNR = %g[dB]\t", SNR))
	}
}

	SNR_th <- 15+10*log10(p)

# Choosing Projective Projection or 
#          projection to p-1 subspace

	if (SNR < SNR_th) {
		if (verbose == 'on') {
			message(sprintf('... Select the projective proj.\n'))
		}
	
		d <- p - 1
		if (snr_input == 0) { # it means that the projection is already computed
			Ud <- Ud[,1:d]
		} else {
			r_m <- rowMeans(R)
			R_m <- replicate(N, r_m)	# mean of each band
			R_o <- R - R_m			# data with zero-mean 

			res_svd <- svd(R_o%*%t(R_o)/N, nu=d, nv=d)	# computes the p-projection matrix 
			Ud <- res_svd$u
			Sd <- res_svd$d
			Vd <- res_svd$v

			x_p <-t(Ud)%*%(R_o)				# project thezeros mean data onto p-subspace
		}

		Rp <- Ud%*%x_p[1:d,] + replicate(N, r_m)	# again in dimension L

		x <- x_p[1:d,]	#  x_p =  Ud' * R_o; is on a p-dim subspace
		c <- max(colSums(x^2))^.5
		y <- rbind(x, c*array(1, dim=c(1,N)))
	} else {
		if (verbose == 'on') {
			message(1,'... Select proj. to p-1\n')
		}

		d <- p

		res_svd <- svd(R%*%t(R)/N, nu=d, nv=d)	# computes the p-projection matrix 
		Ud <- res_svd$u
		Sd <- res_svd$d
		Vd <- res_svd$v

		x_p <- t(Ud)%*%R
		Rp <-Ud%*%x_p[1:d,]	# again in dimension L (note that x_p has no null mean)

		x <- t(Ud) %*% R
		u <- rowMeans(x)	#equivalent to  u = Ud' * r_m
		y <- x / t(replicate(d, colSums( x*replicate(N, u) )))
	}

	# VCA algorithm
	indice <- array(0, dim=c(1,p))
	A <- array(0, dim=c(p,p))
	A[p,1] <- 1

	for (i in 1:p) {
		w <- array(runif(p*1), dim=c(p,1))
		f <- w - A%*%ginv(A)%*%w
		f <- f / sqrt(sum(f^2))

		v <- t(f)%*%y
		v_max <- max(abs(v))
		indice[i] <- which.max(abs(v))
		A[,i] <- y[,indice[i]]	# same as x(:,indice(i))
	}


	Ae <- Rp[,indice]
	res <- NULL
	res$Ae <- Ae
	res$indice<-indice
	res$Rp <- Rp
	return(res)

} #end of vca func


# demo_vca
if (F) {
demo_number <- 2
verbose <- 'on'

if (demo_number == 1) {
	# abundances on strips
	# 3 endmembers (p=3) mixed in the image
	a <- rbind(c(1.0, 0.5, 0.2, 0.0, 0.2, 0.0) , 
              c(0.0, 0.3, 0.7, 1.0, 0.5, 0.0) ,
              c(0.0, 0.2, 0.1, 0.0, 0.3, 1.0))
	p <- dim(a)[1]
	nt<- dim(a)[2]
	Lines <-36
	Columns<-36
	N <- Lines * Columns # number of pixels
	dim_t <- Columns/nt
	s_o <- array(t(replicate(dim_t*Lines , c(t(a)))),dim=c(N, p))
} else if (demo_number == 2 ) {
	p <- 3 # number of endmembers

	# abundance with Dirichlet distributions
	Lines <- 36
	Columns <- 36
	N <- Lines * Columns # number of pixels
	s_o <- dirichlet_rnd(array(rep(1, p)/p ,dim=c(1,p)), N)
} else if (demo_number == 3 ) {
	N <- Lines * Columns
	p <- 6
	A <- array(0, dim=c(L, p))
	s <- array(0, dim=c(p, N))
	x_n <- x
	rm(x)
} else {
	stop('unknown demo number')
}


if (  (demo_number == 1 ) | (demo_number == 2) ) {
	# signatures from USGS Library
	A<-as.matrix(read.table("csv/A.csv",sep=","))
	names<-as.matrix(read.table("csv/names.csv",sep=","))
	wavlen<-as.matrix(read.table("csv/wavlen.csv",sep=","))

	if (verbose == "on") {
		message(sprintf('Selected signatures:\n%s', intToUtf8(t(names))))
	}
	L <- dim(A)[1]
	p_max <- dim(A)[2]
	A = A[,1:p];
	# illumination fluctuation
	q <- 10
	s <- s_o * array(rbeta(N*p, q, 1), dim=c(N,p))
	x <- A %*% t(s)

	# adding noise
	SNR <- 30 # dB
	varianc <-sum(c(x)^2)/10^(SNR/10)/L/N
	n <- sqrt(varianc) * array(rnorm(L*N), dim=c(L,N))
	x_n <- x + n
}

# Unmixing Procedure

if (verbose == 'on') {
	message('Unmixing with VCA algorithm\n')
}

#Useage:
res_vca <- vca(x_n, p, snr_input=1, SNR=SNR)

}


