###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### cpgs function
#### Version v1.1.5
#### 07.05.19
#### Author : Benjamin Planque
#### License CC-BY-NC
###############################################################################################

##### cpgs function
# convex polytope gibbs sampler

cpgs <- function(N,A,b,x0={}){
  p <- dim(A)[2]                         # dimension
  m <- dim(A)[1]                         # num constraint ineqs
  runup <- 50;                          # runup necessary to method
  discard <- runup;
  if (isempty(x0)==TRUE){
    x0 = chebycenter(A,b);              # Chebyshev center of the polytope
  }
  
  # Check input arguments
  if (m < (p+1)){
    stop(cat['At least ',num2str(p+1),' inequalities required',sep=''])
  }
  # Initialisation
  X <- matrix(data=0,nrow = N+runup+discard,ncol = p)           # initialise the output
  n <- 0                                  # num generated so far
  x <- x0 
  
  # Initialize variables for keeping track of sample mean, covariance
  # and isotropic transform.
  M <- matrix(data=0,nrow = p,ncol = 1);                         # Incremental mean.
  S2 <- matrix(data=0,nrow = p,ncol = p);                        # Incremental sum of
  # outer prodcts.
  S <- diag(p);
  T1 <- diag(p);
  W = A;
  
  while (n < (N+runup+discard)){               # sampling loop
    y <- x;
    # compute approximate stochastic transformation
    if (n == runup){
      T0 <- chol(t(S))
      T1 <- t(T0)
      W <- A%*%T1;
    }
    y <- mldivide(T1,y);
    
    # choose p new components
    for (i in 1:p){
      # Find points where the line with the (p-1) components x_i
      # fixed intersects the bounding polytope.
      e <- rep(0,p); e[i] <- 1
      z <- W%*%e
      d <- (b - W%*%y)/z
      tmin <- max(d[z<0]); tmax = min(d[z>0])
      # choose a point on that line segment
      y <- y + (tmin+(tmax-tmin)*runif(1))*e
    }
    x <- T1%*%y
    X[n+1,] <- t(x)
    n <- n + 1
    
    # Incremental mean and covariance updates
    delta0 <- x - M       # delta new point wrt old mean
    M <- M + delta0/n     # sample mean
    delta1 <- x - M;      # delta new point wrt new mean
    if (n > 1){
      S2 <- S2 + (n-1)/(n*n)*(delta0%*%t(delta0))+(delta1%*%t(delta1))
      S0 <- S
      S <- S2/(n-1);           # sample covariance
    } else {
      S <- diag(p)
    }
  }
  
  X <- X[(discard+runup+1):(N+discard+runup),]
  return(X)
}