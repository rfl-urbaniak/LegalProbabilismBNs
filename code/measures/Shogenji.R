



ShogenjiCoherenceForBNs <- function (BN, narrationNodes, states){
startTime <- proc.time() #start measuring computation time
JN <- compile(as.grain(BN)) #convert to junction tree

#generates a table of all prior joint probabilities
PriorJoints <- querygrain(JN,nodes=narrationNodes,type="joint")
#PriorJoints

# calculate the prior of the conjunction of all narration nodes
steps <- numeric(length(narrationNodes))
for(i in 1:length(narrationNodes)){
  steps[i] <- paste(narrationNodes[i], "=", "\"",states[i],"\"")
}
steps<- gsub(" ", "", steps, fixed = TRUE)
final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
narrationPrior <- eval(parse(text=final))
#narrationPrior

#calculate individual priors of the narration nodes
individualProbs <- numeric(length(narrationNodes))
for(i in 1:length(narrationNodes)){
individualProbs[i] <- querygrain(JN,nodes=narrationNodes[i])[[1]][states[i]]
}



stopTime <- proc.time()
elapsedTime <- stopTime - startTime
Result <- list("Narration prior" = narrationPrior, "Node probabilities" = individualProbs, "Shogenji coherence"  = narrationPrior/(prod(individualProbs)),
               "Computation time" = elapsedTime)
Result
}

