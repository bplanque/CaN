# modelling Chance and Necessity: CaN
# An R attempt to build a foodweb assessment model based on CaN principles
#
# Benjamin Planque / Christian Mullon - August 2019
# License CC-BY-NC 4.0


# initialisation ----------------------------------------------------------
graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

require(pracma)                             # Loading the pracma package -- Mathematical tools and functions
require(expm)                               # Loading the expm library -- require to compute matrix powers
require(linprog)
source('./chebycenter.R')
source('./cpgs.R')                          # Convex Polytope Gibbs Sampler (in R)
#require(Rcpp)
#require(RcppEigen)
#require(inline)
#source('./cpgs2.R')                         # Convex Polytope Gibbs Sampler (in cpp)

# loading input data ------------------------------------------------------

# Biomass data - this must include estimates + min and max values (constraints)
# Catches/landings - this must include estimates + min and max values (constraints)

Biomasses <- read.delim(file = 'Biomasses.mini.csv',header = TRUE, sep=';')
Years <- Biomasses[,1]
ny <- length(Years)
Species.names <- colnames(Biomasses[seq(2,dim(Biomasses)[2],3)])
ns0 <- length(Species.names)
Bio.mean <- Biomasses[,seq(2,dim(Biomasses)[2],3)]
Bio.min <- Biomasses[,seq(3,dim(Biomasses)[2],3)]
Bio.max <- Biomasses[,seq(4,dim(Biomasses)[2],3)]

# handling of primary production as Import
Bio.min[,1] <- Bio.min[,1]*0
Import <- matrix(nrow = ny, ncol = ns0, data = 0)
colnames(Import) <- Species.names
Import[,1] <- Bio.mean$PP

Landings=read.delim(file = 'Landings.mini.csv',header = TRUE, sep=';')
if(dim(Biomasses)[1]!=dim(Landings)[1]){
  cat('not the same number of years in the biomass and landing datasets')
}
if(dim(Biomasses)[2]!=dim(Landings)[2]){
  cat('not the same number of species in the biomass and landing datasets')
}
Land.mean <- Landings[,seq(2,dim(Biomasses)[2],3)]
Land.min <- Landings[,seq(3,dim(Biomasses)[2],3)]*0.9
Land.max <- Landings[,seq(4,dim(Biomasses)[2],3)]*1.1

# Construct a new species 'fishery' and update relevant input data
# The idea is that the landings will accumulate in the fishery compartment
# through the years
Species.names=c(Species.names,"Fishery")
ns <- length(Species.names)
Bio.min$Fishery.min <- cumsum(rowSums(Land.min))+sum(Land.min[1,])
Bio.mean$Fishery <- cumsum(rowSums(Land.mean))+sum(Land.mean[1,])
Bio.max$Fishery.max <- cumsum(rowSums(Land.max))+sum(Land.max[1,])
Land.min$Fishery.min <- 0
Land.max$Fishery.max <- 0
Import=cbind(Import,rep(0,ny)); colnames(Import)[ns]='Fishery'

# loading input parameters ------------------------------------------------
# This must include parameters on assimilation efficiency (gamma), 
# digestibility (kappa), other losses (mu), inertia (rho), satiation (sigma)
# and constraints on flows (who eat whom). From these, we can derive alpha, G, H and K.

Input.params=read.delim(file = 'Input.parameters.mini.csv',header = TRUE, sep=';')
if((dim(Input.params)[2]-1)!=ns0){
  cat('not the same number of species in the biomass and parameters datasets')
}

Gamma <- as.numeric(Input.params[1,2:(ns0+1)])
Kappa <- as.numeric(Input.params[2,2:(ns0+1)])
Mu <- as.numeric(Input.params[3,2:(ns0+1)])
Rho <- as.numeric(Input.params[4,2:(ns0+1)])
Sigma <- as.numeric(Input.params[5,2:(ns0+1)])

# update input parameters to include the fishery
Gamma[ns] <- 1      # assimilation efficiency 100%
Kappa[ns] <- 0      # not digestible
Mu[ns] <- 0         # no losses
Rho[ns] <- 100      # no inertia
Sigma[ns] <- 1000   # (almost) never satiated

Diets <- read.delim(file = 'Trophic.flows.mini.csv',header = TRUE, sep = ';') # species trophic interactions
# add the Fishery
Fished.species <- names(Land.mean)[colSums(Land.mean)!=0]
Diets <- rbind(Diets,data.frame(Prey=Fished.species,Predator=rep('Fishery',length(Fished.species))))

# Construction of the diet matrix
# Note: species interaction matrices are always constructed with prey in rows and predators in columns
# When unfolded (i.e. vectorized) the Fluxes are taken first along the prey axis (going down the rows) 
# and then along the predator axis (going to the right in the columns)
x1 <- data.frame(sp.num=1:ns,sp.name=Species.names)
x2 <- merge(Diets,x1,by.x = "Prey", by.y = "sp.name",all.x = TRUE)
x3 <- merge(x2,x1,by.x="Predator",by.y="sp.name",all.x=TRUE)
DietMatrix <- matrix(nrow = ns,ncol = ns,data = 0)
DietMatrix[x3[,3]+(x3[,4]-1)*ns] <- 1
DietVector <- as.vector(DietMatrix)
w=which(DietMatrix==1,arr.ind = TRUE)
DietNames <- paste(Species.names[w[,1]],'-->',Species.names[w[,2]],sep='')

# Construction of the min and max fluxes, using landings as input
F.min <- rep(0,length=ns*ns*ny) # No fluxes can be lower than zero
F.min[rep(c(rep(0,(ns*(ns-1))),rep(1,ns)),ny)==1] <- as.vector(t(as.matrix(Land.min)))
F.max <- rep(max(max(Biomasses)),length=ns*ns*ny) # No fluxes can be greater than the largest biomass ever observed
F.max[rep(c(rep(0,(ns*(ns-1))),rep(1,ns)),ny)==1] <- as.vector(t(as.matrix(Land.max)))

# deriving G, H and K -----------------------------------------------------
# so that the Mater equation can be written in the form:
# B(i,t+1) - B(i,t) = Sumj(G(j,i)F(j,i,t) - Sumj(K(i)F(i,j,t)) + K(i)Import(i,t) - H(i)B(i,t)
H <- 1-exp(-Mu)
K <- H/Mu
K[Mu==0]=1 # special case when Mu=0 (fishery)
G <- (Kappa%*%t(K*Gamma))
G[,ns] <- 1 # special case for the fishery (no digestibility issue)
G <- G*DietMatrix

# Construction of matrix N (note (e) from equation 6) --------------------- C'est ici que �a merde?
IE <- diag(ns)                                     # identity matrix of size E (number of species)
Nijk <- matrix(data=0,nrow = ns,ncol = ns^2)       # interaction betwee species i and flux j->k
for (i in 1:ns){                                   # loop on species (i)
  deltaij <- matrix(data=0,nrow = ns,ncol = ns)    # initialise matrix delta(i,j)
  deltaij[i,] <- 1                                 # delta(i,j)=1 if the prey (row) is equal to i
  deltaik <- t(deltaij)                            # delta(i,k)=1 if the predator (col) is equal to i
  Ni <- G*deltaik-K[i]*deltaij                     # compute the sub-matrix N(j,k) for species i
  Nijk[i,] <- as.vector(Ni)                        # linearlise it and store it into the matrix N(i,jk)
}

# Screen summary of input data and parameters ------------------------------------
names(Gamma) <- Species.names
names(Kappa) <- Species.names
names(Mu) <- Species.names
names(Rho) <- Species.names
names(Sigma) <- Species.names
names(H) <- Species.names
names(K) <- Species.names
colnames(G) <- Species.names
rownames(G) <- Species.names
Nijk2 <- t(Nijk)
Nijk2 <- Nijk2[DietVector==1,]
rownames(Nijk2) <- DietNames
colnames(Nijk2) <- Species.names
print('Bio.min')
print(Bio.min)
print('Bio.max')
print(Bio.max)
print('B0')
print(B0)
print('Land.min')
print(Land.min)
print('Land.max')
print(Land.max)
print('Gamma')
print(Gamma)
print('Kappa')
print(Kappa)
print('Sigma')
print(Sigma)
print('Mu')
print(Mu)
print('G')
print(G)
print('H')
print(H)
print('K')
print(K)
print('N')
print(Nijk2)

# Matrix formulation of the dynamics (Appendix B.3) -----------------------
# Construct L1, L2, and M so that B=L.F+M (equation 14)
# Construction of matrix L1 (equation 15)
L0 <- matrix(data = 0,nrow = ny,ncol = ny)          # initialise the L0 matrix (exponents required to compute L)
for (i in 1:ny){                                    # loop on years
  L0[(row(L0)-i)==col(L0)] <- i                   # assign exponent
}
L <- matrix(data=0,nrow = ns*ny,ncol = ns*ns*ny) # 
for (j in 1:ny){                                # loop on columns (years)
  for (i in j:ny){                              # loop on rows (years)
    Lij <- ((IE-IE*H)%^%L0[i,j])%*%Nijk                # construct the submatrix L(i,j)
    L[((i-1)*ns+1):(i*ns),((j-1)*(ns^2)+1):(j*(ns^2))] <- Lij
  }
}
L=L[,rep(DietVector,ny)==1] # only retain the columns of L for the posible fluxes
rownames(L)=rep(Species.names,ny)
colnames(L)=rep(DietNames,ny)

# Construction of vector M (equation 16 revised to include import)
# M=M1+M2, with M1 the part related to B0 and M2 the part related to imports
B0 <- Biomasses[1,seq(2,dim(Biomasses)[2],3)]
B0[1]=Biomasses[1,1]*K[1]
B0$Fishery <- 0
M1 <- matrix(data=0,nrow = dim(Nijk)[1]*ny,ncol = 1)
for (i in 1:ny){                                    # loop on rows (years)
  M1[((i-1)*ns+1):(i*ns)] <- as.vector(((1-H)^i)*t(B0))
}

# Construction of Matrix M2 (added)
L2 <- matrix(data=0,nrow = ns*ny,ncol = ns*ny) 
for (j in 1:ny){                                # loop on columns (years)
  for (i in j:ny){                              # loop on rows (years)
    Lij <- ((IE-IE*H)%^%L0[i,j])%*%(IE*K)                # construct the submatrix L(i,j)
    L2[((i-1)*ns+1):(i*ns),((j-1)*ns+1):(j*ns)] <- Lij
  }
}
M2 <- L2%*%as.vector(t(Import))

M <- M1+M2

# Matrix formulation of the constraints (Appendix B.4) --------------------
# so that all constraints can be expressed in the form A.F<=B
# constraint #1: flows are positive F>=0 (i.e. -F<=0)
# !!! Note that this constraint is redundant with constraint #8
#A1 <- diag(ns*ns*ny)*-1 # identity matrix of size fluxes * time-steps (with value -1)
#b1 <- rep(0,ns*ns*ny) # vector of zeroes of length fluxes * time-steps
#A1 <- A1[rep(DietVector,ny)==1,]; A1 <- A1[,rep(DietVector,ny)==1] # reduce A1 to possible fluxes only
#row.names(A1) <- colnames(L); colnames(A1) <- colnames(L);
#b1 <- b1[rep(DietVector,ny)==1] # reduce B1 to possible fluxes only
#names(b1) <- row.names(A1)
A1=NULL
b1=NULL

# constraint #2: biomasses are positive L.F+M>=0 (i.e. -L.F<=M) 
# !!! Note that this constraint is irrelevant as it is superceeded by the constraint #4
A2 <- -L # the L matrix is not specified properly
b2 <- M1; row.names(b2) <- rownames(A2)
#A2=NULL
#b2=NULL

# constraint #3: biomasses are bounded above: L.F+M<=Bio.max (i.e. L.F<=Bio.max-M)
A3 <- L
b3 <- as.vector(t((Bio.max)))-as.vector(M1); names(b3) <- rownames(A3)

# constraint #4: biomasses are bounded below: L.F+M>=Bio.min (i.e. -L.F<=M-Bio.min)
A4 <- -L
b4 <- as.vector(M1)-as.vector(t((Bio.min))); names(b4) <- rownames(A4)

# constraint #5: variations in biomass are bounded above: 
A5=NULL
b5=NULL

# constraint #6: variations in biomass are bounded below: 
A6=NULL
b6=NULL

# constraint #7: flows are bouned above: F<=F.max
A7 <- diag(ns*ns*ny)
b7 <- F.max
A7 <- A7[(rep(DietVector,ny)==1),];A7 <- A7[,rep(DietVector,ny)==1]
row.names(A7) <- colnames(L); colnames(A7) <- colnames(L);
b7 <- b7[(rep(DietVector,ny)==1)] # reduce A7 and b7 to possible flows only
names(b7) <- row.names(A7)

# constraint #8: flows are bouned below: F>=F.min (i.e. -F<=-F.min)
A8 <- diag(ns*ns*ny)*-1
b8 <- -F.min
A8 <- A8[(rep(DietVector,ny)==1),];A8 <- A8[,rep(DietVector,ny)==1]
row.names(A8) <- colnames(L); colnames(A8) <- colnames(L);
b8 <- b8[(rep(DietVector,ny)==1)] # reduce A8 and b8 to possible flows only
names(b8) <- row.names(A8)

# constraint #9: satiation. Sum(F(.,i))<=sigma*B(i,t)
VI<-as.matrix(rep(1:ns,byrow=F,each=ns))
c<-matrix(1:ns)
e<-matrix(1,1,ns*ns)
Fij<-c%*%e==d%*%t(VI)                         # Matrix product of c-e and d-VI
SumFij<-apply(Fij,c(1,2),as.numeric)
A9=matrix(0,nrow = ns*ny,ncol = ns*ns*ny)
for (i in 1:ny){
  for (j in i:ny){
    A9[((i-1)*ns+1):(i*ns),((i-1)*ns^2+1):(i*ns^2)] <- SumFij
  }
}
A9 <- A9[,rep(DietVector,ny)==1] # reduce A9 to possible fluxes only
row.names(A9) <- row.names(L)
colnames(A9) <- colnames(L)
# b9...this one is going to be tricky to write. 


# Combining all constraints together
A <- rbind(A1,A2,A3,A4,A5,A6,A7,A8)
b <- c(b1,b2,b3,b4,b5,b6,b7,b8)

# reducing the dimensions of A and b to well defined constraints
Alines <- apply(abs(A), 1,sum)>0                   # Sum up the matrix to keep the rows
Ap <- A[which(Alines==TRUE),]                      # Select the lines that do not contain zeroes only
bp <- as.matrix(b[which(Alines==T)])               # Select elements in the vector b that are left in matrix Ap

# sampling of the polytope ------------------------------------------------
F0<-chebycenter(rbind(A7,A8),c(b7,b8))                   # Need for a starting point x0 : chebycenter method to define it
F0<-chebycenter(rbind(A1,A2,A7,A8),c(b1,b2,b7,b8))                   # Need for a starting point x0 : chebycenter method to define it

F0<-chebycenter(rbind(A1,A3,A7,A8),c(b1,b3,b7,b8))                   # Need for a starting point x0 : chebycenter method to define it

F0<-chebycenter(Ap,bp)                   # Need for a starting point x0 : chebycenter method to define it
# check that F0 meets the constraints
print(prod(Ap%*%F0<=bp))

#Fsample<-cpgs(10,Ap,bp,F0)              # Sample with Gibbs algorithm 100 vectors of flows


# saving results ----------------------------------------------------------


