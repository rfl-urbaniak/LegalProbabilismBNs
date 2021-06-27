library(tidyverse)
library(plyr)
library(rje)
library(bnlearn)
library(utils)
library(useful)



#use this:
#P(A v B) =  P(A & B) + P(A & ~ B) + P(~A, B)


OlssonCoherenceForBNs <- function(BN, narrationNodes, states){
  
startTime <- proc.time()
JN <- compile(as.grain(BN))

# this calculates the probability of the disjunction of narration nodes in a BN
options <- c("1","0")
fullTable <- as.data.frame(expand.grid(rep(list(options),length(narrationNodes))))
names(fullTable) <- narrationNodes

#adding disjunction
for (valuation in 1:nrow(fullTable)){
result <- logical(length(narrationNodes))
    for(node in 1:length(narrationNodes)){
  result[node] <- states[node] == fullTable[valuation,node]
  }
fullTable$disjunction[valuation] <- max(result)
}

# pick rows where the disjunction is true
summands <- fullTable[fullTable$disjunction==1,]

#summands

#generates a table of all prior joint probabilities
PriorJoints <- querygrain(JN,nodes=narrationNodes,type="joint")
#PriorJoints

# get the narration prior
steps <- numeric(length(narrationNodes))
for(i in 1:length(narrationNodes)){
  steps[i] <- paste(narrationNodes[i], "=", "\"",states[i],"\"")
}
steps<- gsub(" ", "", steps, fixed = TRUE)
final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
narrationPrior <- eval(parse(text=final))


#picks the priors for the description states in which the disjunction is true
for (valuation in 1:nrow(summands)){
steps <- numeric(length(narrationNodes))
for(i in 1:length(narrationNodes)){
  steps[i] <- paste(narrationNodes[i], "=", "\"", summands[valuation,i],"\"")
}
steps<- gsub(" ", "", steps, fixed = TRUE)
final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
summands$probs[valuation] <-  eval(parse(text=final))
}


stopTime <- proc.time()
elapsedTime <- stopTime - startTime
Result <- list("Narration prior" = narrationPrior, "Disjunction probability" = sum(summands$probs), "Olsson coherence"  = narrationPrior/sum(summands$probs),
               "Computation time" = elapsedTime)
return(Result)
}






