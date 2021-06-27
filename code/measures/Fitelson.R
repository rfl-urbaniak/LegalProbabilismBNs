FitelsonConfirmation <- function (eIfh, eIfnh, status = "Ind"){
  if (status == "Ent") {
    1
  } else if (status == "Ref") {
    -1 
  } else { 
    (eIfh - eIfnh) /  (eIfh + eIfnh)
    }
  }


#_______________________________
#Fitelson measure
#arguments will include
# - BN
# - narrationNodes
# - states of narrationNodes
# - logical status left-to-right (a vector from "Ind","Ent","Ref"of length Fitelson(length(narrationNodes)))
# - logical status right-to-left (a similar vector)
# - to build the last two, first run 
# DisjointPairs(narrationNodes) 
# then build the two vectors by hand, unless all pairs are logically independent

FitelsonCoherenceForBNs <- function (BN, narrationNodes, states,
                        statusLR = rep("Ind",FitelsonSize(length(narrationNodes))),
                        statusRL = rep("Ind",FitelsonSize(length(narrationNodes)))
                          ){
    
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
  combinations$statusLR <- statusLR
  combinations$statusRL <- statusRL
  combinations$PriorL <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PriorR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrLR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrRL <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrNLR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$PrNRL <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$FitLR <- numeric(FitelsonSize(length(narrationNodes)))
  combinations$FitRL <- numeric(FitelsonSize(length(narrationNodes)))
  
  
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
    
    
    #for right-if-not left we use Bayes theorem:
    # P(R|~L) =   P(~L|R) P(R)/ P(~L)
    # =    (1 - P(L|R))  P(R) /  (1 - P(L) )
    combinations$PrNLR[j] <- ((1 - combinations$PrRL[j]) * combinations$PriorR[j]) /
      (1 - combinations$PriorL[j])         
    
    #symmetrically for left-if-not right:
    combinations$PrNRL[j] <- ((1 - combinations$PrLR[j]) * combinations$PriorL[j]) /
      (1 - combinations$PriorR[j])         
    
    #finally, Fitelson measures! notice left support for right, so we 
    #look at the probablity of left (e) given right (h), of left (e) given not right (~h), 
    #keeping track of the logical relationship between left (e) and right (h).
    combinations$FitLR[j] <- FitelsonConfirmation(combinations$PrRL[j],combinations$PrNRL[j],combinations$statusLR[j])
    
    #symmetrically in the opposite direction
    combinations$FitRL[j] <- FitelsonConfirmation(combinations$PrLR[j],combinations$PrNLR[j],combinations$statusRL[j])
    
  }
  
  #now stop the calculations and the timer, and spit out everything
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  Result <- list("Narration Prior" = narrationPrior, "Pairs and measures" = combinations, "Fitelson coherence"  = mean(c(combinations$FitLR,combinations$FitRL)),
                 "Computation time" = elapsedTime)
  Result
}


