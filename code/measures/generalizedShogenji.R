ShogenjiGeneralizedCoherence <- function(BN, narrationNodes, states){
  
  startTime <- proc.time() #start measuring computation time
 
   sets <- NonSingletonSubsets(narrationNodes)
  ShogenjiCoherences <- numeric(length(sets))
  for (i in 1:length(sets)){
    subset <- sets[[i]]
    values <- states[which(narrationNodes %in% subset)]
    ShogenjiCoherences[i] <- ShogenjiCoherenceForBNs(BN,subset,values)[[3]]
  }
  ShogenjiGeneralizedScore <- mean(ShogenjiCoherences)
    
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  
  return(list("Shogenji generalized score" = ShogenjiGeneralizedScore,
              "Computation time" = elapsedTime))
}


#ShogenjiGeneralizedCoherence(BirdBN, c("B","G","P"), c("1","1","1"))

#ShogenjiGeneralizedCoherence(BirdBN, c("B","G"), c("1","1"))

#ShogenjiGeneralizedCoherence(BirdBN, c("B","P"), c("1","1"))





