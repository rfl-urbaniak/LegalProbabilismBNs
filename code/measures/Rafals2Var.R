## These are my weird attempts at better definions of binary coherence
#Define RA
RA <- function(AIfB,BIfA,BAndA){
  mean(c(Z(AIfB,BAndA),Z(BIfA,BAndA)))
}

#Defina RA extended
RAE <- function(AIfB,BIfA,BAndA,A,B){
  mean(c(Z(AIfB,BAndA),Z(BIfA,BAndA),
         Z(AIfB,A),Z(BIfA,B)))
}

#Define RAFixed
RAF <- function(AIfB,BIfA,A,B){
  mean(c(Z(AIfB,A),Z(BIfA,B)))
}
