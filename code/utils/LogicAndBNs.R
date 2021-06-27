
#__________________________________________________________________
#produces dataframe with joint probabilities of nodes of interest
StateDescriptionsDF <- function(BN, narrationNodes){
  
  JN <- compile(as.grain(BN))  
  options <- c("1","0")
  fullTable <- as.data.frame(expand.grid(rep(list(options),length(narrationNodes))))
  names(fullTable) <- narrationNodes
  #generates a table of all prior joint probabilities
  PriorJoints <- querygrain(JN,nodes=narrationNodes,type="joint")
  #PriorJoints
  
   for (valuation in 1:nrow(fullTable)){
    #picks the prior of a narration specified by nodes and states
    steps <- numeric(length(narrationNodes))
    for(i in 1:length(narrationNodes)){
      steps[i] <- paste(narrationNodes[i], "=", "\"", fullTable[valuation,i],"\"")
    }
    steps<- gsub(" ", "", steps, fixed = TRUE)
    final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    fullTable$probs[valuation] <-  eval(parse(text=final))
   }
  return(fullTable)
}


#_________________________________________
#Calculates the probability that at least one of listed nodes have the corresponding state

DisjunctionBN <- function (BN, narrationNodes, states){

JN <- compile(as.grain(BN))

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

for (valuation in 1:nrow(summands)){
  #picks the prior of a narration specified by nodes and states
  steps <- numeric(length(narrationNodes))
  for(i in 1:length(narrationNodes)){
    steps[i] <- paste(narrationNodes[i], "=", "\"", summands[valuation,i],"\"")
  }
  steps<- gsub(" ", "", steps, fixed = TRUE)
  final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
  summands$probs[valuation] <-  eval(parse(text=final))
}

return(sum(summands$probs))
}



#this turned out unnecessary; I keep it just in case
#negateAllStates <- function (states){
#  nStates <- ifelse(states == "1", "0", "1")
#  return(nStates)
#}


#update BN

Nodes <- "W1"
States <- "1"

updateBN <- function(BN,Nodes,States){
  JN <- compile(as.grain(BN))
  JNu <- setEvidence(JN, nodes = Nodes, states = States)
  BNu <- as.bn.fit(JNu,including.evidence = TRUE)
  return(BNu)
}

