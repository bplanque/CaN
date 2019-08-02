###############################################################################################
#### The Non-Deterministic Network Model (NDND) Mullon et al., 2009 , Planque et al., 2014
#### cpgs2 function
#### Version v1.0
#### 24.06.19
#### Author : hilaire Drouineau
#### License ???
###############################################################################################

## Author : Hilaire Drouineau and Benjamin Planque

# library(RcppEigen)
# library(Rcpp)
# library(inline)
# library(linprog)

#creation de la fonction CPP    
cpgsCPP='#include <Eigen/Cholesky>
#include <Rcpp.h>
using namespace Rcpp;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

const Map<MatrixXd>A(as<Map<MatrixXd> >(AA));
const Map<MatrixXd>b(as<Map<MatrixXd> >(bb));
const Map<MatrixXd>x0(as<Map<MatrixXd> >(xx0));

int N=Rcpp::as<int>(NN); 
int p=A.cols();
int m=A.rows();
int runup=100;
int discard=100;

// Check input arguments
if (m < (p+1)){
throw std::range_error("dimensions mismatch");
}
// Initialisation

MatrixXd X(N+runup+discard,p);
int n=0;
MatrixXd x(p,1);
MatrixXd y(p,1);
x=x0;

// Initialize variables for keeping track of sample mean, covariance
// and isotropic transform.
MatrixXd M(p,1);
M.setZero();
MatrixXd S2(p,p);
S2.setZero();

// outer products.
MatrixXd S(p,p);
MatrixXd S0(p,p);
S.setIdentity();

MatrixXd T1(p,p);
T1.setIdentity();
MatrixXd W(m,p);

W = A;
MatrixXd d(m,1);
MatrixXd delta0(p,1);
MatrixXd delta1(p,1);
MatrixXd z(m,1);
while (n < (N+runup+discard)){               //sampling loop
y=x;
NumericVector alea2=runif(p);
// compute approximate stochastic transformation
if (n == runup){
T1=S.transpose().llt().matrixU().transpose();
W = A*T1;
}
y=T1.inverse()*y;

// choose p new components
VectorXd e(p);
for (int i=0;i<p;++i){
//Find points where the line with the (p-1) components x_i
//fixed intersects the bounding polytope.
e.setZero();
e(i)= 1;
z = (W*e); //prevent any divisions by 0

d=(b - W*y);
d=d.cwiseQuotient(z); 
double tmin=-9e9;
double tmax=9e9;
for (int j=0;j<m;++j){
if (z(j)<0 && tmin<d(j)) tmin=d(j);
if (z(j)>0 && tmax>d(j)) tmax=d(j);
}
y(i)+=(tmin+(tmax-tmin)*alea2(i));
}
x=T1*y;
X.row(n)= x.col(0);
++n;
// Incremental mean and covariance updates
delta0 = x - M; // delta new point wrt old mean
M+= delta0/(double)n;     // sample mean
delta1= x - M;      // delta new point wrt new mean


if (n > 1){
S2 +=(n-1)/(double)(n*n)*(delta0*delta0.transpose())+(delta1*delta1.transpose());
S0 = S;
S = S2/(double)(n-1);           // sample covariance
} else {
S.setIdentity();
}
}

return wrap(X.bottomRows(N));'

#compilation (à ne faire qu'une seule fois)
settings=getPlugin("RcppEigen")
settings$env$PKG_CXXFLAGS <- paste('-O3 -fopenmp', settings$env$PKG_CXXFLAGS)
settings$env$PKG_LIBS <- paste('-fopenmp -lgomp', settings$env$PKG_LIBS)

cpgs2<-cxxfunction(signature(NN="integer",AA="matrix",bb="vector",xx0="vector"),cpgsCPP,plugin="RcppEigen",verbose=TRUE,
                   settings=settings)