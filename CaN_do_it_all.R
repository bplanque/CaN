# modelling Chance and Necessity: CaN
# An R attempt to build a foodweb assessment model based on CaN principles
#
# Benjamin Planque - August 2019


# initialisation ----------------------------------------------------------
graphics.off()                              # clear all graphical windows
cat("\014")                                 # clear console
rm(list=ls())                               # clear the work environment

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

Biomasses=read.delim(file = 'Biomasses.csv',header = TRUE, sep=';')
Years=Biomasses[,1]
ny=length(Years)
Species.names=colnames(Biomasses[seq(2,dim(Biomasses)[2],3)])
ns=length(Species.names)
Bio.min=Biomasses[,seq(3,dim(Biomasses)[2],3)]
Bio.max=Biomasses[,seq(4,dim(Biomasses)[2],3)]

Landings=read.delim(file = 'Landings.csv',header = TRUE, sep=';')
if(dim(Biomasses)[1]!=dim(Landings)[1]){
  cat('not the same number of years in the biomass and landing datasets')
}
if(dim(Biomasses)[2]!=dim(Landings)[2]){
  cat('not the same number of species in the biomass and landing datasets')
}
Land.min=Landings[,seq(3,dim(Biomasses)[2],3)]
Land.max=Landings[,seq(4,dim(Biomasses)[2],3)]

# loading input parameters ------------------------------------------------
# This must include parameters on assimilation efficiency (gamma), 
# digestibility (kappa), other losses (mu), inertia (rho), satiation (sigma)
# and constraints on flows (who eat whom). From these, we can derive alpha, G, H and K.

Input.params=read.delim(file = 'Input.parameters.csv',header = TRUE, sep=';')
if((dim(Input.params)[2]-1)!=ns){
  cat('not the same number of species in the biomass and parameters datasets')
}

Gamma=as.numeric(Input.params[1,2:(ns+1)])
Kappa=as.numeric(Input.params[2,2:(ns+1)])
Mu=as.numeric(Input.params[3,2:(ns+1)])
Rho=as.numeric(Input.params[4,2:(ns+1)])
Sigma=as.numeric(Input.params[5,2:(ns+1)])
# alpha added by hand for now
Alpha=c(100,100,100,0.25,100,100,0.5,0.5)

Diets <- read.delim(file = 'Trophic.flows.csv',header = TRUE, sep = ';')
# Construction of the diet matrix
x1=data.frame(sp.num=1:ns,sp.name=Species.names)
x2=merge(Diets,x1,by.x = "Prey", by.y = "sp.name",all.x = TRUE)
x3=merge(x2,x1,by.x="Predator",by.y="sp.name",all.x=TRUE)
DietMatrix=matrix(nrow = ns,ncol = ns,data = 0)
DietMatrix[x3[,3]+(x3[,4]-1)*ns]=1
DietVector=as.vector(DietMatrix)

# deriving G, H and K
H=1-exp(-Mu)
K=H/Mu
G=(Kappa%*%t(H*Gamma/Mu))*DietMatrix

# construction of the polytope --------------------------------------------
# from the above, build matrix A and vector b so that A.f=<b, with f, the  
# vector of all trophic flows at all time steps

# Construction of matrix L
Lijk=matrix(data=0,nrow = ns,ncol = ns^2) # je suis perdu!
for (i in 1:ns){ # loop on species (i)
  deltaij=deltaik=matrix(data = 0,nrow = ns,ncol = ns)
  deltaij[,i]=1
  deltaik[]
  Ljk={}
}

IE=diag(ns) # identity matrix of size E (number of species)
L0=matrix(data = NA,nrow = ny,ncol = ny)
for (i in 1:(ny-1)){
  L0[(row(L0)-i)==col(L0)]<-i-1
}


# Construction of matrix M


# sampling of the polytope ------------------------------------------------


# saving results ----------------------------------------------------------


