### Ok, now let's make the whole thing a function
RAcoherenceForBNs <- function(BN,narrationNodes,states){
  startTime <- proc.time()
  #compiles the BN as a junction tree
  JN <- compile(as.grain(BN))
  
  #generates a table of all prior joint probabilities
  PriorJoints <- querygrain(JN,nodes=narrationNodes,type="joint")
  #PriorJoints
  
  #picks the prior of a narration specified by nodes and states
  steps <- numeric(length(narrationNodes))
  for(i in 1:length(narrationNodes)){
    steps[i] <- paste(narrationNodes[i], "=", "\"",states[i],"\"")
  }
  steps<- gsub(" ", "", steps, fixed = TRUE)
  final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
  narrationPrior <- eval(parse(text=final))
  
  
  #now for each node i, pick the i-th node with its corresponding state as specified, 
  #add as evidence and create a posterior junction tree, save the posteriors to a vector
  narrationPosteriors <- numeric(length(states))
  
  for(node in 1:length(narrationNodes)){
    
    JNwithEvidence <- setEvidence(JN,nodes = narrationNodes[node],
                                  states = states[node])  
    
    PosteriorJoints <-  querygrain(JNwithEvidence,nodes = narrationNodes,type="joint")
    nodesPosterior <- narrationNodes[-node]
    statesPosterior <- states[-node]
    
    stepsPosterior <- numeric(length(nodesPosterior))
    for(i in 1:length(nodesPosterior)){
      stepsPosterior[i] <- paste(nodesPosterior[i], "=", "\"",statesPosterior[i],"\"")
    }
    stepsPosterior<- gsub(" ", "", stepsPosterior, fixed = TRUE)
    finalPosterior <- paste("PosteriorJoints[",paste(stepsPosterior,collapse=","),"]",sep="")
    narrationPosteriors[node] <- eval(parse(text=finalPosterior))
  }
  #narrationPosteriors
  
  
  Zs <- numeric(length(narrationNodes))
  for(z in 1:length(narrationNodes)){
    Zs[z] <- Z(posterior = narrationPosteriors[z], prior = narrationPrior)
  }
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  Result <- list("Narration Prior" = narrationPrior, "Narration Posteriors" = narrationPosteriors, "RA coherence" = mean(Zs), "Z values" = Zs,  "Computation time" = elapsedTime)
Result
}
