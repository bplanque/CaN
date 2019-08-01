# modelling Chance and Necessity: CaN
# An R attempt to build a foodweb assessment model based on CaN principles
#
# Benjamin Planque - August 2019


# initialisation ----------------------------------------------------------
graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

require(expm)                                 # Loading the expm library -- require to compute matrix powers
# library(abind)                              # Loading the abind package -- Bindind arrays
# library(arrayhelpers)                       # Loading the arrayhelpers package -- Functions for array modification and shaping
# library(devtools)                           # Loading the devtools package -- Tools for package development
# library(ggdendro)                           # Loading the ggdendro package -- Package to plot dendrograms with ggplot
# library(ggplot2)                            # Loading the ggplot package -- Graphics package
# library(inline)                             # Loading the inline package -- Read C++ functions called in R
# library(LIM)                                # Loading the LIM package -- The LIM library is used for the sampling of flows -- Soetaert and van Oevelen (2014)
# library(linprog)                            # Loading the linprog package -- Linear programming and optimization
# library(memisc)                             # Loading the memisc package -- 
# library(NbClust)                            # Loading the NbClust package -- Package to define the optimal number of cluster in a dataset
# library(plyr)                               # Loading the plyr package -- Dataframe manipulation
# library(pracma)                             # Loading the pracma package -- Mathematical tools and functions
# library(Rcpp)                               # Loading the Rcpp package -- Integration of R and C++
# library(RcppEigen)                          # Loading the RcppEigen package -- Linear algebra on C++ integrated in R
# library(rgl)                                # Loading the rgl package -- 3d graphs
# library(rstudioapi)                         # Loading the rstudioapi package -- Directory and file interactive choice
# library(scales)                             # Loading the scales package -- Graphical attributes package
# library(tidyverse)                          # Loading the tidyverse package -- Set of packages for data manipulations

# loading input data ------------------------------------------------------

# Biomass data - this must include estimates + min and max values (constraints)
# Catches/landings - this must include estimates + min and max values (constraints)

Biomasses <- read.delim(file = 'Biomasses.csv',header = TRUE, sep=';')
Years <- Biomasses[,1]
ny <- length(Years)
Species.names <- colnames(Biomasses[seq(2,dim(Biomasses)[2],3)])
ns <- length(Species.names)
Bio.min <- Biomasses[,seq(3,dim(Biomasses)[2],3)]
Bio.max <- Biomasses[,seq(4,dim(Biomasses)[2],3)]
B0 <- Biomasses[1,seq(2,dim(Biomasses)[2],3)]

Landings=read.delim(file = 'Landings.csv',header = TRUE, sep=';')
if(dim(Biomasses)[1]!=dim(Landings)[1]){
  cat('not the same number of years in the biomass and landing datasets')
}
if(dim(Biomasses)[2]!=dim(Landings)[2]){
  cat('not the same number of species in the biomass and landing datasets')
}
Land.min <- Landings[,seq(3,dim(Biomasses)[2],3)]
Land.max <- Landings[,seq(4,dim(Biomasses)[2],3)]

# loading input parameters ------------------------------------------------
# This must include parameters on assimilation efficiency (gamma), 
# digestibility (kappa), other losses (mu), inertia (rho), satiation (sigma)
# and constraints on flows (who eat whom). From these, we can derive alpha, G, H and K.

Input.params=read.delim(file = 'Input.parameters.csv',header = TRUE, sep=';')
if((dim(Input.params)[2]-1)!=ns){
  cat('not the same number of species in the biomass and parameters datasets')
}

Gamma <- as.numeric(Input.params[1,2:(ns+1)])
Kappa <- as.numeric(Input.params[2,2:(ns+1)])
Mu <- as.numeric(Input.params[3,2:(ns+1)])
Rho <- as.numeric(Input.params[4,2:(ns+1)])
Sigma <- as.numeric(Input.params[5,2:(ns+1)])
# alpha added by hand for now
Alpha <- c(100,100,100,0.25,100,100,0.5,0.5)

Diets <- read.delim(file = 'Trophic.flows.csv',header = TRUE, sep = ';')
# Construction of the diet matrix
# Note: species interaction matrices are always constructed with prey in rows and predators in columns
x1 <- data.frame(sp.num=1:ns,sp.name=Species.names)
x2 <- merge(Diets,x1,by.x = "Prey", by.y = "sp.name",all.x = TRUE)
x3 <- merge(x2,x1,by.x="Predator",by.y="sp.name",all.x=TRUE)
DietMatrix <- matrix(nrow = ns,ncol = ns,data = 0)
DietMatrix[x3[,3]+(x3[,4]-1)*ns] <- 1
DietVector <- as.vector(DietMatrix)


# deriving G, H and K -----------------------------------------------------
# so that the Mater equation can be written in the form:
# B(i,t+1) − B(i,t) = Sumj(G(j,i)F(j,i,t)−Sumj(K(i)F(i,j,t)) − H(i)B(i,t)
H <- 1-exp(-Mu)
K <- H/Mu
G <- (Kappa%*%t(K*Gamma))*DietMatrix

# construction of the polytope --------------------------------------------
# from the above, build matrix A and vector b so that A.f=<b, with f, the  
# vector of all trophic flows at all time steps

# matrix version of H (note (c) from equation 6)
IE <- diag(ns)                                     # identity matrix of size E (number of species)
H <- IE*H                                          # matrix with diagonal terms Hi 

# Construction of matrix J (note (e) from equation 6)
Jijk <- matrix(data=0,nrow = ns,ncol = ns^2)       # interaction betwee species i and flux j->k
for (i in 1:ns){                                   # loop on species (i)
  deltaij <- matrix(data=0,nrow = ns,ncol = ns)    # initialise matrix delta(i,j)
  deltaij[i,] <- 1                                 # delta(i,j)=1 if the prey (row) is equal to i
  deltaik <- t(deltaij)                            # delta(i,k)=1 if the predator (col) is equal to i
  Ji <- G*deltaik-K[i]*deltaij                     # compute the sub-matrix J(j,k) for species i
  Jijk[i,] <- as.vector(Ji)                        # linearlise it and store it into the matrix J(i,jk)
}

# Matrix formulation of the dynamics (Appendix B.3)
# Construct L, and M so that B=L.F+M (equation 14)

# Construction of matrix L (equation 15)
L0 <- matrix(data = 0,nrow = ny,ncol = ny)        # initialise the L0 matrix (exponents required to compute L)
for (i in 1:(ny-1)){                            # loop on years
  L0[(row(L0)-i)==col(L0)] <- i-1                 # assign exponent
}
L <- matrix(data=0,nrow = dim(Jijk)[1]*ny,ncol = dim(Jijk)[2]*ny) 
for (j in 1:(ny-1)){                              # loop on columns (years)
  for (i in (j+1):ny){                          # loop on rows (years)
    Lij <- ((IE-H)%^%L0[i,j])%*%Jijk                 # construct the submatrix L(i,j)
    L[((i-1)*ns+1):(i*ns),((i-1)*(ns^2)+1):(i*(ns^2))] <- Lij
  }
}

# Construction of vector M (equation 16)
M <- matrix(data=0,nrow = dim(Jijk)[1]*ny,ncol = 1)
for (i in 1:ny){                          # loop on rows (years)
  M[((i-1)*ns+1):(i*ns)] <- ((IE-H)%^%(i-1))%*%t(as.vector(B0))
}

# Matrix formulation of the constraints (Appendix B.4)
# so that all constraints can be expressed in the form A.F<=B
# constraint #1: flows are positive F>=0
A1 <- diag(ns*ns*ny)*-1
B1 <- rep(0,ns*ns*ny)
A1 <- A1[rep(DietVector,ny),]; A1 <- A1[,rep(DietVector,ny)]
# That looks ok but the problem is that the landings are not included in the flow matrices above

# constraint #2: biomasses are positive L.F+M>=0
# constraint #3: biomasses are bounded above: L.F+M<=Bio.max
# constraint #4: biomasses are bounded below: L.F+M>=Bio.min
# constraint #5: variations in biomass are bounded above: 
# constraint #6: variations in biomass are bounded below: 
# constraint #7: some flows (landings) are bouned above: F<=Land.max
# constraint #8: some flows (landings) are bouned below: F>=Land.min

# sampling of the polytope ------------------------------------------------


# saving results ----------------------------------------------------------


