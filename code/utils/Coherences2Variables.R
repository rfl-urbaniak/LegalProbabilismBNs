source("measures/CoherenceMeasures2Var.R")

Coherences2Var <- function(A,B,AAndB,AOrB,AIfNotB,BIfNotA,BIfA){
  #defined probabilities
  AIfB <- BIfA * A / B
  Prob <- AAndB #question: why isn't this redundant? why not use AAndB later?
  #Calculate coherences
  Ra <- RA(AIfB, BIfA,AAndB)
  Raf <- RAF(AIfB,BIfA,A,B)
  Fi <- Fitelson(AIfB,AIfNotB,BIfA,BIfNotA)
  Ro <- Roche(AIfB,BIfA)
  Ol <- Olsson(AAndB,AOrB)
  Sh <- Shogenji(AAndB,A,B)
  Do <- Douven(AIfB,BIfA,A,B)
  Coherences <- c(Prob,Ra,Raf,Fi,Ro,Ol,Sh,Do)
  names(Coherences) <- c("Prior","RA","RAF","Fitelson","Roche","Olsson","Shogenji","DM")
  return(Coherences)
}

Coherences4Scenarios <- function(nameA, nameB,
                                 A, B, AAndB, AOrB,
                                 NotBIfNotA, BIfA){
  #defined probabilities
  NotA <- 1-A
  NotB <- 1-B
  NotAIfNotB <- NotBIfNotA * NotA / NotB
  AIfNotB <- 1 - NotAIfNotB
  AIfB <- BIfA * A / B
  NotBIfA <- AIfNotB * NotB / A
  BIfNotA <- 1 - NotBIfNotA
  NotAIfB <- BIfNotA * NotA / B
  NotBIfNotA <- 1 - BIfNotA
  NotAAndB <- B - AAndB
  AAndNotB <- A - AAndB
  NotAAndNotB <- 1 - AOrB
  NotAOrB <- NotA + B - NotAAndB
  AOrNotB <- A + NotB - AAndNotB
  NotAOrNotB <- NotA + NotB - NotAAndNotB
  
  Scenarios <- list()
  
  #A and B
  Scenarios[[paste(nameA,"and",nameB)]] <- Coherences2Var(A,B,AAndB,AOrB,
    AIfNotB,BIfNotA,BIfA)
  
  #Not A and B
  Scenarios[[paste("Not",nameA,"and",nameB)]] <- Coherences2Var(NotA,B,NotAAndB,NotAOrB, NotAIfNotB,BIfA,BIfNotA)
  
  #A and Not B  
  Scenarios[[paste(nameA,"and Not",nameB)]] <- Coherences2Var(A,NotB,AAndNotB,AOrNotB,AIfB,NotBIfNotA,NotBIfA)
  
  #Not A and Not B
  Scenarios[[paste("Not",nameA,"and Not",nameB)]] <- Coherences2Var(NotA,NotB,NotAAndNotB,NotAOrNotB,NotAIfB,NotBIfA,NotBIfNotA)
  
  ScenariosTable <- as.data.frame(t(as.data.frame(Scenarios)))
  ScenariosTable
}
