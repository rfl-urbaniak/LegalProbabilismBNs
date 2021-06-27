differenceConfirmation <- function (hIfe, h){
     hIfe - h
    }

#_______________________________
#Roche measure
#arguments will include
# - BN
# - narrationNodes
# - states of narrationNodes

RocheCoherenceForBNs <- function (BN, narrationNodes, states){
  
  startTime <- proc.time() #start measuring computation time
  JN <- compile(as.grain(BN)) #convert to junction tree
  
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
  #narrationPrior
  
  # initiates the results table
  combinations <- DisjointPairs(narrationNodes) #sets of nodes to test
  combinations$PriorL <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PriorR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrLR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrRL <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$DiffLR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$DiffRL <- numeric(FitelsonSize(length(narrationNodes)))
  
  
  #these should be identical
  #FitelsonSize(length(narrationNodes)) == nrow(combinations)
  
  
  #calculate probabilities
  for (j in 1:nrow(combinations)){   #go through rows
    left <- combinations$Var1[j][[1]]  #pick the two vectors from combinations
    right <- combinations$Var2[j][[1]]
    
    # this sorts the vectors and their states in the order in which they occur in the narration nodes; this is needed to read off probs from joint probability tables later on.
    leftOrdered <- narrationNodes[sort(sapply(left,placeOfNode, 
                                              narrationNodes = narrationNodes),  decreasing=FALSE)]
    
    leftStatesOrdered <- states[sort(sapply(left,placeOfNode, 
                                            narrationNodes = narrationNodes),  decreasing=FALSE)]
    
    rightOrdered <- narrationNodes[sort(sapply(right,placeOfNode, 
                                               narrationNodes = narrationNodes), decreasing=FALSE)]
    
    rightStatesOrdered <- states[sort(sapply(right,placeOfNode, 
                                             narrationNodes = narrationNodes),  decreasing=FALSE)]
    
    #we will need priors of both vectors
    PriorL <- querygrain(JN,nodes=leftOrdered,type="joint")
    PriorR <- querygrain(JN,nodes=rightOrdered,type="joint")
    
    #we also need to instantiate  left and ask about  right
    JNleft <- setEvidence(JN, nodes = leftOrdered, states = leftStatesOrdered)
    JNleftTable <- querygrain(JNleft, nodes = rightOrdered,type = "joint")
    
    # and instantiate right and ask about left
    JNright <- setEvidence(JN, nodes = rightOrdered, states = rightStatesOrdered)
    JNrightTable <- querygrain(JNright, nodes = leftOrdered,type = "joint")
    
    #calculate the prior of left
    #everywhere: this mess with steps etc. turns a vector of states like c("1","0") for nodes like c("A","B")
    #into a way of picking the right value from the  right joint probability table.
    stepsL <- numeric(length(leftOrdered))
    for(i in 1:length(leftOrdered)){
      stepsL[i] <- paste(leftOrdered[i], "=", "\"",leftStatesOrdered[i],"\"")
    }
    stepsL <- gsub(" ", "", stepsL, fixed = TRUE)
    final <- paste("PriorL[",paste(stepsL,collapse=","),"]",sep="")
    eval(parse(text=final))
    combinations$PriorL[j] <- eval(parse(text=final))
    
    
    #calculate prior of right
    stepsR <- numeric(length(rightOrdered))
    for(i in 1:length(rightOrdered)){
      stepsR[i] <- paste(rightOrdered[i], "=", "\"",rightStatesOrdered[i],"\"")
    }
    stepsR <- gsub(" ", "", stepsR, fixed = TRUE)
    final <- paste("PriorR[",paste(stepsR,collapse=","),"]",sep="")
    combinations$PriorR[j] <- eval(parse(text=final))
    
    
    #calculate right if left
    stepsLR <- numeric(length(rightOrdered))
    stepsLR
    for(i in 1:length(rightOrdered)){
      stepsLR[i] <- paste(rightOrdered[i], "=", "\"",rightStatesOrdered[i],"\"")
    }
    stepsLR <- gsub(" ", "", stepsLR, fixed = TRUE)
    final <- paste("JNleftTable[",paste(stepsLR,collapse=","),"]",sep="")
    eval(parse(text=final))
    combinations$PrLR[j] <- eval(parse(text=final))
    
    
    #calculate left if right
    stepsRL <- numeric(length(leftOrdered))
    for(i in 1:length(leftOrdered)){
      stepsRL[i] <- paste(leftOrdered[i], "=", "\"",leftStatesOrdered[i],"\"")
    }
    stepsRL <- gsub(" ", "", stepsRL, fixed = TRUE)
    final <- paste("JNrightTable[",paste(stepsRL,collapse=","),"]",sep="")
    combinations$PrRL[j] <- eval(parse(text=final))
    
    
            
    
   
    combinations$DiffLR[j] <- combinations$PrLR[j] - combinations$PriorR[j]
    
    #symmetrically in the opposite direction
    combinations$DiffRL[j] <- combinations$PrRL[j] - combinations$PriorL[j]
  }
  
  #now stop the calculations and the timer, and spit out everything
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  Result <- list("Narration Prior" = narrationPrior, "Pairs and measures" = combinations, "Roche coherence"  = mean(c(combinations$DiffLR,combinations$DiffRL)),
                 "Computation time" = elapsedTime)
  Result
}


