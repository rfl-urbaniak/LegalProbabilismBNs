#source("measures/Rafals2Var.R")

#Define confirmation measure Z
# Z <- function(posterior,prior){
#   d <- posterior - prior
# #  if (prior == 1 || prior == 0) 0 else   #This is if we want to define extreme cases by hand
#   ifelse(posterior >= prior, d/(1-prior), d/prior)
# }


Z <- function(posterior,prior){
  d <- posterior - prior
  ifelse(prior == posterior, 0, ifelse(posterior > prior, d/(1-prior), d/prior))
}




#Define Fitelson's measure

Fitelson <- function(AIfB,AIfNotB,BIfA,BIfNotA){
  (((AIfB - AIfNotB)/(AIfB + AIfNotB))
  +
  ((BIfA - BIfNotA)/(BIfA + BIfNotA)))/2
}


#Roche's average mutual firmness

Roche <- function(AIfB,BIfA){
  (AIfB + BIfA)/2 
}


#Olsson's relative overlap

Olsson <- function(AAndB,AOrB){
  AAndB/AOrB
}

#Shogenji's coherence

Shogenji <- function(AAndB,A,B){
  AAndB/(A * B)
}


#Douven and Meij's mesure

Douven <- function(AIfB,BIfA,A,B){
  (AIfB - A + BIfA - B)/2
}


