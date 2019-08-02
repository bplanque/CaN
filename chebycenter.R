###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### Chebycenter function
#### Version v1.0
#### 24.06.19
#### Author : Benjamin Planque, Hilaire Drouineau
###############################################################################################

chebycenter <- function(A,b){
  n <- dim(A)[1]
  p <- dim(A)[2];
  an <- rowSums(A^2)^0.5
  A1 <- matrix(data=0,nrow = n,ncol = p+1)
  A1[,1:p] <- A;
  A1[,p+1] <- an;
  f <- matrix(data=0,nrow=p+1,ncol=1)
  f[p+1] <- -1;
  
  #d <- linprog(cc=f,A=A1,b=b);
  #x <- d$x[1:p];
  
  d <- lp(direction="min",objective.in = f,const.mat = A1,const.rhs = as.numeric(b),const.dir = rep("<=",n))
  x <- d$solution
  return(x[-p-1])
}