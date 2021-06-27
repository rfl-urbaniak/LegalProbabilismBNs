OlssonGeneralizedCoherence <- function(BN, narrationNodes, states){

startTime <- proc.time() #start measuring computation time

sets <- NonSingletonSubsets(narrationNodes)
OlssonCoherences <- numeric(length(sets))
for (i in 1:length(sets)){
subset <- sets[[i]]
values <- states[which(narrationNodes %in% subset)]
OlssonCoherences[i] <- OlssonCoherenceForBNs(BN,subset,values)[[3]]
}
OlssonGeneralizedScore <- mean(OlssonCoherences)

stopTime <- proc.time()
elapsedTime <- stopTime - startTime

return(list("Olsson generalized score" = OlssonGeneralizedScore, "Computation time" = elapsedTime))
}

#test
#OlssonGeneralizedCoherence(BirdBN, c("B","G","P"), c("1","1","1"))

#OlssonGeneralizedCoherence(BirdBN, c("B","G"), c("1","1"))

#OlssonGeneralizedCoherence(BirdBN, c("B","P"), c("1","1"))

