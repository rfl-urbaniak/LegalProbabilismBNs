#These will be used in Fitelson & DouvenMeijs


# Powerset without the empty set
NonEmptySubsets <- function(vector){
  power <- powerSet(vector)
  subsets <- power[-1]
  subsets
}


#test
# NonEmptySubsets(1:5)


#All pairs of subsets of a vector
PairsOfSubsets <- function (vector) {
  as.data.frame(expand.grid(NonEmptySubsets(vector),NonEmptySubsets(vector)))
}


#test
# pairs <- PairsOfSubsets(1:3)
# pairs

NonSingletonSubsets <- function(narrationNodes){
  subsets <- powerSet(narrationNodes)
  idsToRemove <- sapply(subsets, function(i) length(i) < 2)
  nonSingletons <- subsets[!idsToRemove]
  return(nonSingletons)
}



#this does the same, but eliminates pairs repeated in a different order
UniquePairsOfSubsets <- function (vector) {
  p <- PairsOfSubsets(vector)  
  p$code <- 1:nrow(p)
  right <- logical(nrow(p))
  left <- logical(nrow(p))
  pick <- logical(nrow(p))
  for(i in 1:nrow(p)){
    for (j in 1:nrow(p)){
      right[j] <- compare.list(p$Var1[i],p$Var2[j])
      left[j] <- compare.list(p$Var2[i],p$Var1[j])
      double <- which(right & left == TRUE)
      pick[i] <- min(i,double)
    }
  }
  return(p[unique(pick),])
}



#test
# UniquePairsOfSubsets(1:3)




#list of unique disjoint pairs of subsets of a vector
DisjointPairs <- function (vector){
  pairs <- UniquePairsOfSubsets(vector)
  common <- numeric(nrow(pairs))
  for(i in 1:nrow(pairs)){
    common[i] <- length(intersect(pairs[i,1][[1]],pairs[i,2][[1]]))
  }
  disjoint <- pairs[common==0,]
  disjoint
}


#test
# DisjointPairs(c("A","B","C"))
# DisjointPairs(1:3)


### How the number of pairs involved grows
# Let set S have k elements. A pair of disjoint subsets is formed by examining each
# element of S and deciding whether to put it into one subset, into the other subset, 
# or into neither.   Three choices for each of k elements gives us the term 3^k as the
# count of all possible ordered pairs of disjoint subsets that can be formed.
# To count only unordered pairs, group pairs into equivalence classes. Each pair is 
# one of 2! permutations, except the pair ⟨∅,∅⟩. ( The null set is disjoint with itself).
# But we don't want this case. Also, we don't want to count those which inolve an empty 
# set, and there are 2^k-1 of those. 

#this calculates by formula the number of pairs needed
FitelsonSize <- function(k){
  (3^k -1)/2 - (2^k-1)
} 
# test
# FitelsonSize(1:5)
#This is used in vector limits in some definitions in Fitelson and Douven & Meijs



#inspect
# x <- 1:20
# y <- FitelsonSize(x)
# ggplot()+geom_line(aes(x=x,y=y))
# FitelsonSize(10)


# utilities for BNs that associate nodes with their states in a given narration etc.
placeOfNode <- function(node, narrationNodes){
  which(narrationNodes == node)
}


#find the value assigned to the node by the narration
stateOfNode <- function(node, narrationNodes, states){
  states[placeOfNode(node, narrationNodes)]
}


selectStatesOfaSet <- function (set){
  states[which(narrationNodes %in% set)]
}



#find states of a node 
findStates <- function(node){
  eval(parse(text = paste('dimnames(BN$',node,'$prob)$',node, sep= "")))
}


findStatesBN <- function (node, BNname = "BN"){
  eval(parse(text = paste('dimnames(',BNname,'$',node,'$prob)$',node, sep= "")))
}


#find parents of a node

parents <- function(net,node){
  bnlearn::parents(net,node)
}


#finds joint priors


FindPriorJN <- function(JN,nodes,states){
PriorJoints <- querygrain(JN,nodes=nodes,type="joint")
    PriorJoints <- aperm(PriorJoints, nodes)
    steps <- numeric(length(nodes))
    for(rn in 1:length(nodes)){
      steps[rn] <- paste(nodes[rn], "=", "\"",states[rn],"\"")
    }
    steps<- gsub(" ", "", steps, fixed = TRUE)
    final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    noquote(final)
    prior <- eval(parse(text=final))
    return(prior)
  }
  
  
  
#   ##OUTDATED
#   PriorJoints <- querygrain(JN,nodes=nodes,type="joint")  
#   steps <- numeric(length(nodes))
#   for(i in 1:length(nodes)){
#     steps[i] <- paste(nodes[i], "=", "\"",states[i],"\"")
#   }
#   steps<- gsub(" ", "", steps, fixed = TRUE)
#   final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
#   prior <- eval(parse(text=final))
#   return(prior)
# }



FindPriorBN  <- function(BN,nodes,states){
  JN <- compile(as.grain(BN)) 
  PriorJoints <- querygrain(JN,nodes=nodes,type="joint")
  PriorJoints <- aperm(PriorJoints, nodes)
  steps <- numeric(length(nodes))
  for(rn in 1:length(nodes)){
    steps[rn] <- paste(nodes[rn], "=", "\"",states[rn],"\"")
  }
  steps<- gsub(" ", "", steps, fixed = TRUE)
  final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
  noquote(final)
  prior <- eval(parse(text=final))
  return(prior)
}







