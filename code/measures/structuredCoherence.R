
Z <- function(posterior,prior){
  d <- posterior - prior
  ifelse(prior == posterior, 0, ifelse(posterior > prior, d/(1-prior), d/prior))
}




structuredNoSD <- function(BN, narrationNodes, states){
  startTime <- proc.time() #start measuring computation time
  JN <- compile(as.grain(BN))
  
  #find non-root nodes
  parented <- unique(arcs(BN)[,2])
  
  parented
  
  #assign parents to parented
  parentList <- list()
  for(node in 1:length(parented)){
    parentList[[node]] <- bnlearn::parents(BN,parented[node])  
  }
  
  
  #initiate lists of results
  expConfFull <- list()
  ECSanteS <- numeric(length(parented))
  
  
  for (i in 1:length(parented)){
    consequent <- parented[i]
    
    
    consequentStates <- if (consequent %in% narrationNodes){
      stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)    
    } else {
      findStates(consequent)
    }
    
    
    antecedents  <- parentList[i]
    
    
    antecedentStates <- list()
    for(a in 1:length(antecedents[[1]])){
      antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
        stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes, 
                    states = states)    
      } else if (!(antecedents[[1]][a] %in% parented)){
        eval(parse(text = paste('dimnames(BN$',antecedents[[1]][a],'$prob)', 
                                sep= "")))[[1]]
      } else {
        findStates(antecedents[[1]][a])
      }
    }
    
    
    
    #start the outcome table
    variants <-  expand.grid(c(list(consequentStates),antecedentStates))
    colnames(variants) <- c(consequent,antecedents[[1]])
    
    
    #add the prior of consequent, to be used for Z calculation
    pr <- numeric(nrow(variants))
    for (s in 1:nrow(variants)){
      pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
    }
    variants <- cbind(variants, priorCons = pr)
    
    variants
    
    #posteriors 
    posteriors <- numeric(nrow(variants))
    for (row in 1:nrow(variants)){
      JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
      posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
    }
    variants <- cbind(variants,posteriors)
    
    
    #variants
    #priors for conjunctions
    #priorJoint <- numeric(nrow(variants))
    # 
    # for (row in 1:nrow(variants)){
    #   rowNodes <- as.vector(unlist(c(consequent,antecedents)))
    #   rowStates <- as.vector(unlist(variants[row,1:(length(antecedents[[1]])+1)]))
    #   PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
    #   PriorJoints <- aperm(PriorJoints, rowNodes)
    #   steps <- numeric(length(rowNodes))
    #   for(rn in 1:length(rowNodes)){
    #     steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
    #   }
    #   steps<- gsub(" ", "", steps, fixed = TRUE)
    #   final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    #   prior <- eval(parse(text=final))
    #   priorJoint[row] <- prior
    # }
    # 
    # variants$PriorJoint <- priorJoint
    # 
    #priors for the antecedents
    
    
    priorAnte <- numeric(nrow(variants))
    
    for (row in 1:nrow(variants)){
      rowNodes <- as.vector(unlist(c(antecedents)))
      rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
      PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
      PriorJoints <- aperm(PriorJoints, rowNodes)
      steps <- numeric(length(rowNodes))
      for(rn in 1:length(rowNodes)){
        steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
      }
      steps<- gsub(" ", "", steps, fixed = TRUE)
      final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
      noquote(final)
      prior <- eval(parse(text=final))
      priorAnte[row] <- prior
    }
    variants$PriorAnte <- priorAnte
    
    # #scaled weights for conjunction
    # if(sum(variants$PriorJoint) > 0){
    #   variants$WeightsCon <-  variants$PriorJoint / sum(variants$PriorJoint)
    # } else {
    #   variants$WeightsCon <-  1/nrow(variants)
    # }
    # 
    
    
    #scaled weights for antecedent
    if(sum(variants$PriorAnte) > 0){
      variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
    } else {
      variants$WeightsAnte <-  1/nrow(variants)
    }
    
    
    #Z measures
    variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
    #variants$Zweighted <- variants$Z * variants$Weights
    #variants$ZweightedCon <- variants$Z * variants$PriorJoint
    #variants$ZweightedAnte <- variants$Z * variants$PriorAnte
    #variants$ZscaledCon <- variants$Z * variants$WeightsCon
    variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
    
    
    
    expConfFull[[i]] <- list( "Consequent node" = consequent,
                              "Options & calculations" = variants,
                              "ECSante scaled" = sum(variants$ZscaledAnte))
    
    ECSanteS[i] <- sum(variants$ZscaledAnte)
  }
  
  
  populationSD <- function( vector ){
    sqrt(sum((vector - mean(vector))^2)/(length(vector)))
  }
  
  structuredScore <- function(expConf)  if (min(expConf) <= 0) {
    (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
  } else {
    (mean(expConf) - populationSD(expConf))
  }
  
  structuredScoreSquared <- function(ECS)  if (min(ECS) <= 0) {
    (mean(ECS) - populationSD(ECS)^2) * (min(ECS) +1) - min(ECS)^2
  } else {
    (mean(ECS) - populationSD(ECS)^2)
  }
  
  
  structuredScoreNoSD <- function(ECS)  if (min(ECS) <= 0) {
    (mean(ECS)) * (min(ECS) +1) - min(ECS)^2
  } else {
    (mean(ECS))
  }
  
  
  
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  
  return(list("Full calculations" = expConfFull, 
              "ECSanteS" = ECSanteS,
              "structured Coherence Ante scaled" = structuredScore(ECSanteS),
              "structuredSquared" = structuredScoreSquared(ECSanteS),
              "structuredNoSD" = structuredScoreNoSD(ECSanteS),
              "Computation time" = elapsedTime))
}

















#----------------------------


structuredCoherenceUpdated <- function(BN, narrationNodes, states){
  startTime <- proc.time() #start measuring computation time
  JN <- compile(as.grain(BN))
  
  #find non-root nodes
  parented <- unique(arcs(BN)[,2])
  
  parented
  
  #assign parents to parented
  parentList <- list()
  for(node in 1:length(parented)){
    parentList[[node]] <- bnlearn::parents(BN,parented[node])  
  }
  
  
  #initiate lists of results
  expConfFull <- list()
  ECSanteS <- numeric(length(parented))
  
  
  for (i in 1:length(parented)){
    consequent <- parented[i]
    
    
    consequentStates <- if (consequent %in% narrationNodes){
      stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)    
    } else {
      findStates(consequent)
    }
    
    
    antecedents  <- parentList[i]
    

    antecedentStates <- list()
    for(a in 1:length(antecedents[[1]])){
      antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
        stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes, 
                    states = states)    
      } else if (!(antecedents[[1]][a] %in% parented)){
        eval(parse(text = paste('dimnames(BN$',antecedents[[1]][a],'$prob)', 
                                sep= "")))[[1]]
      } else {
        findStates(antecedents[[1]][a])
      }
    }
  
    
    
    #start the outcome table
    variants <-  expand.grid(c(list(consequentStates),antecedentStates))
    colnames(variants) <- c(consequent,antecedents[[1]])
    
    
    #add the prior of consequent, to be used for Z calculation
    pr <- numeric(nrow(variants))
    for (s in 1:nrow(variants)){
      pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
    }
    variants <- cbind(variants, priorCons = pr)
    
    variants
    
    #posteriors 
    posteriors <- numeric(nrow(variants))
    for (row in 1:nrow(variants)){
      JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
      posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
    }
    variants <- cbind(variants,posteriors)
    
    
    variants
    #priors for conjunctions
    #priorJoint <- numeric(nrow(variants))
    # 
    # for (row in 1:nrow(variants)){
    #   rowNodes <- as.vector(unlist(c(consequent,antecedents)))
    #   rowStates <- as.vector(unlist(variants[row,1:(length(antecedents[[1]])+1)]))
    #   PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
    #   PriorJoints <- aperm(PriorJoints, rowNodes)
    #   steps <- numeric(length(rowNodes))
    #   for(rn in 1:length(rowNodes)){
    #     steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
    #   }
    #   steps<- gsub(" ", "", steps, fixed = TRUE)
    #   final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    #   prior <- eval(parse(text=final))
    #   priorJoint[row] <- prior
    # }
    # 
    # variants$PriorJoint <- priorJoint
    # 
    #priors for the antecedents
    
    
    priorAnte <- numeric(nrow(variants))
    
    for (row in 1:nrow(variants)){
      rowNodes <- as.vector(unlist(c(antecedents)))
      rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
      PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
      PriorJoints <- aperm(PriorJoints, rowNodes)
      steps <- numeric(length(rowNodes))
      for(rn in 1:length(rowNodes)){
        steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
      }
      steps<- gsub(" ", "", steps, fixed = TRUE)
      final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
      noquote(final)
      prior <- eval(parse(text=final))
      priorAnte[row] <- prior
    }
    variants$PriorAnte <- priorAnte
    
    # #scaled weights for conjunction
    # if(sum(variants$PriorJoint) > 0){
    #   variants$WeightsCon <-  variants$PriorJoint / sum(variants$PriorJoint)
    # } else {
    #   variants$WeightsCon <-  1/nrow(variants)
    # }
    # 
    
    
    #scaled weights for antecedent
    if(sum(variants$PriorAnte) > 0){
      variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
    } else {
      variants$WeightsAnte <-  1/nrow(variants)
    }
    
    
    #Z measures
    variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
    #variants$Zweighted <- variants$Z * variants$Weights
    #variants$ZweightedCon <- variants$Z * variants$PriorJoint
    #variants$ZweightedAnte <- variants$Z * variants$PriorAnte
    #variants$ZscaledCon <- variants$Z * variants$WeightsCon
    variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
    
    
    
    expConfFull[[i]] <- list( "Consequent node" = consequent,
                              "Options & calculations" = variants,
                              "ECSante scaled" = sum(variants$ZscaledAnte))
    
    ECSanteS[i] <- sum(variants$ZscaledAnte)
  }
  
  
  populationSD <- function( vector ){
    sqrt(sum((vector - mean(vector))^2)/(length(vector)))
  }
  
  structuredScore <- function(expConf)  if (min(expConf) <= 0) {
    (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
  } else {
    (mean(expConf) - populationSD(expConf))
  }
  
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  
  return(list("Full calculations" = expConfFull, 
              "ECSanteS" = ECSanteS,
              "structured Coherence Ante scaled" = structuredScore(ECSanteS),
              "Computation time" = elapsedTime))
}



























#__________________THE REAL DEAL

# narrationNodes <-  c("O2","TT","C")
# states <- rep("1",3)
# BN <- D4BN
# 
# structuredCoherence(BN,narrationNodes,states)


structuredCoherence <- function(BN, narrationNodes, states){
  startTime <- proc.time() #start measuring computation time
  JN <- compile(as.grain(BN))
  
  #find non-root nodes
  parented <- unique(arcs(BN)[,2])
  
  parented
  
  #assign parents to parented
  parentList <- list()
  for(node in 1:length(parented)){
    parentList[[node]] <- bnlearn::parents(BN,parented[node])  
  }
  
  
  #initiate lists of results
  expConfFull <- list()
  ECSanteS <- numeric(length(parented))
  
  
  
  for (i in 1:length(parented)){
    consequent <- parented[i]
    
    consequentStates <- if (consequent %in% narrationNodes){
      stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)    
    } else {
      findStates(consequent)
    }
    
    antecedents  <- parentList[i]
    antecedents
    
    antecedentStates <- list()
    for(a in 1:length(antecedents[[1]])){
      antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
        stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes, states = states)    
      } else {
        findStates(antecedents[[1]][a])
      }
    }
    
    #start the outcome table
    variants <-  expand.grid(c(list(consequentStates),antecedentStates))
    colnames(variants) <- c(consequent,antecedents[[1]])
    
    #add the prior of consequent, to be used for Z calculation
    pr <- numeric(nrow(variants))
    for (s in 1:nrow(variants)){
      pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
    }
    variants <- cbind(variants, priorCons = pr)
    
    #posteriors 
    posteriors <- numeric(nrow(variants))
    for (row in 1:nrow(variants)){
      JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
      posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
    }
    variants <- cbind(variants,posteriors)
    
    #priors for conjunctions
    #priorJoint <- numeric(nrow(variants))
    # 
    # for (row in 1:nrow(variants)){
    #   rowNodes <- as.vector(unlist(c(consequent,antecedents)))
    #   rowStates <- as.vector(unlist(variants[row,1:(length(antecedents[[1]])+1)]))
    #   PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
    #   PriorJoints <- aperm(PriorJoints, rowNodes)
    #   steps <- numeric(length(rowNodes))
    #   for(rn in 1:length(rowNodes)){
    #     steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
    #   }
    #   steps<- gsub(" ", "", steps, fixed = TRUE)
    #   final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    #   prior <- eval(parse(text=final))
    #   priorJoint[row] <- prior
    # }
    # 
    # variants$PriorJoint <- priorJoint
    # 
    #priors for the antecedents
    
    
    priorAnte <- numeric(nrow(variants))
    
    for (row in 1:nrow(variants)){
      rowNodes <- as.vector(unlist(c(antecedents)))
      rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
      PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
      PriorJoints <- aperm(PriorJoints, rowNodes)
      steps <- numeric(length(rowNodes))
      for(rn in 1:length(rowNodes)){
        steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
      }
      steps<- gsub(" ", "", steps, fixed = TRUE)
      final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
      noquote(final)
      prior <- eval(parse(text=final))
      priorAnte[row] <- prior
    }
    variants$PriorAnte <- priorAnte
    
    
    # #scaled weights for conjunction
    # if(sum(variants$PriorJoint) > 0){
    #   variants$WeightsCon <-  variants$PriorJoint / sum(variants$PriorJoint)
    # } else {
    #   variants$WeightsCon <-  1/nrow(variants)
    # }
    # 
    
    
    #scaled weights for antecedent
    if(sum(variants$PriorAnte) > 0){
      variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
    } else {
      variants$WeightsAnte <-  1/nrow(variants)
    }
    
    
    #Z measures
    variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
    #variants$Zweighted <- variants$Z * variants$Weights
    #variants$ZweightedCon <- variants$Z * variants$PriorJoint
    #variants$ZweightedAnte <- variants$Z * variants$PriorAnte
    #variants$ZscaledCon <- variants$Z * variants$WeightsCon
    variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
    
    
    
    expConfFull[[i]] <- list( "Consequent node" = consequent,
                              "Options & calculations" = variants,
                              "ECSante scaled" = sum(variants$ZscaledAnte))
    
    ECSanteS[i] <- sum(variants$ZscaledAnte)
  }
  
  #expConf
  
  populationSD <- function( vector ){
    sqrt(sum((vector - mean(vector))^2)/(length(vector)))
  }
  
  structuredScore <- function(expConf)  if (min(expConf) <= 0) {
    (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
  } else {
    (mean(expConf) - populationSD(expConf))
  }
  
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  
  return(list("Full calculations" = expConfFull, 
              "ECSanteS" = ECSanteS,
              "structured Coherence Ante scaled" = structuredScore(ECSanteS),
              "Computation time" = elapsedTime))
}

































#______________________
### This is with conjunction etc. 

structuredCoherence3 <- function(BN, narrationNodes, states){
startTime <- proc.time() #start measuring computation time
JN <- compile(as.grain(BN))

#find non-root nodes
parented <- unique(arcs(BN)[,2])

#assign parents to parented
parentList <- list()
for(node in 1:length(parented)){
  parentList[[node]] <- bnlearn::parents(BN,parented[node])  
}


#initiate lists of results
expConfFull <- list()
ECScon <- numeric(length(parented))
ECSante <- numeric(length(parented))
ECSconS <- numeric(length(parented))
ECSanteS <- numeric(length(parented))



for (i in 1:length(parented)){
  consequent <- parented[i]
  
  consequentStates <- if (consequent %in% narrationNodes){
    stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)    
  } else {
    findStates(consequent)
  }
  
  antecedents  <- parentList[i]
  antecedents
  
  antecedentStates <- list()
  for(a in 1:length(antecedents[[1]])){
    antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
      stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes, states = states)    
    } else {
      findStates(antecedents[[1]][a])
    }
  }
  
  #start the outcome table
  variants <-  expand.grid(c(list(consequentStates),antecedentStates))
  colnames(variants) <- c(consequent,antecedents[[1]])
  
  #add the prior of consequent, to be used for Z calculation
  pr <- numeric(nrow(variants))
  for (s in 1:nrow(variants)){
    pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
  }
  variants <- cbind(variants, priorCons = pr)
  
  #posteriors 
  posteriors <- numeric(nrow(variants))
  for (row in 1:nrow(variants)){
    JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
    posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
  }
  variants <- cbind(variants,posteriors)
  
  #priors for conjunctions
  priorJoint <- numeric(nrow(variants))
  
  for (row in 1:nrow(variants)){
    rowNodes <- as.vector(unlist(c(consequent,antecedents)))
    rowStates <- as.vector(unlist(variants[row,1:(length(antecedents[[1]])+1)]))
    PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
    PriorJoints <- aperm(PriorJoints, rowNodes)
  steps <- numeric(length(rowNodes))
    for(rn in 1:length(rowNodes)){
      steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
    }
    steps<- gsub(" ", "", steps, fixed = TRUE)
    final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    prior <- eval(parse(text=final))
    priorJoint[row] <- prior
  }
  
  variants$PriorJoint <- priorJoint
  
  #priors for the antecedents
  priorAnte <- numeric(nrow(variants))
  
  for (row in 1:nrow(variants)){
    rowNodes <- as.vector(unlist(c(antecedents)))
    rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
    PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
    PriorJoints <- aperm(PriorJoints, rowNodes)
    steps <- numeric(length(rowNodes))
    for(rn in 1:length(rowNodes)){
      steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
    }
    steps<- gsub(" ", "", steps, fixed = TRUE)
    final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
    noquote(final)
    prior <- eval(parse(text=final))
    priorAnte[row] <- prior
  }
  variants$PriorAnte <- priorAnte
  
  
  #scaled weights for conjunction
  if(sum(variants$PriorJoint) > 0){
  variants$WeightsCon <-  variants$PriorJoint / sum(variants$PriorJoint)
  } else {
  variants$WeightsCon <-  1/nrow(variants)
  }
  
  #scaled weights for antecedent
  if(sum(variants$PriorAnte) > 0){
    variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
  } else {
    variants$WeightsAnte <-  1/nrow(variants)
  }
  
  
  #Z measures
  variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
  #variants$Zweighted <- variants$Z * variants$Weights
  variants$ZweightedCon <- variants$Z * variants$PriorJoint
  variants$ZweightedAnte <- variants$Z * variants$PriorAnte
  variants$ZscaledCon <- variants$Z * variants$WeightsCon
  variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
  
  
  
  expConfFull[[i]] <- list( "Consequent node" = consequent,
                        "Options & calculations" = variants, "ECScon" = sum(variants$ZweightedCon),
                        "ECSante" = sum(variants$ZweightedAnte), "ECScon scaled" = sum(variants$ZscaledCon),
                        "ECSante scaled" = sum(variants$ZscaledAnte))
  
  ECScon[i]   <-  sum(variants$ZweightedCon)
  ECSante[i]  <- sum(variants$ZweightedAnte) 
  ECSconS[i]  <- sum(variants$ZscaledCon)
  ECSanteS[i] <- sum(variants$ZscaledAnte)
}

#expConf

populationSD <- function( vector ){
  sqrt(sum((vector - mean(vector))^2)/(length(vector)))
}

structuredScore <- function(expConf)  if (min(expConf) <= 0) {
  (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
} else {
  (mean(expConf) - populationSD(expConf))
}

stopTime <- proc.time()
elapsedTime <- stopTime - startTime

return(list("Full calculations" = expConfFull, "ECScon" = ECScon, "ECSante" = ECSante,
            "ECSconS" = ECSconS,
            "ECSanteS" = ECSanteS,
            "structured Coherence Con" = structuredScore(ECScon), 
            "structured Coherence Ante" = structuredScore(ECSante),
            "structured Coherence Con scaled" = structuredScore(ECSconS),
            "structured Coherence Ante scaled" = structuredScore(ECSanteS),
            "Computation time" = elapsedTime))
}



#_____________________ CLEAN UP THE DROPPED ONES, ADD WEIGHING ACCORDING TO NARRATION

#BN <- robbersLargeBN
#narrationNodes <- c("P","R")
#states <- c("1","1")


# 
# robbersTable2 <- CoherencesTable2(robbersLargeBN,
#                                   scenariosList = list(c("P","R"),c("P","R"),c("P","R")),
#                                   statesList   = list(c("1","1"),c("1","0"),c("0","1")),
#                                   exampleName = "Robbers"
# )


#structuredCoherence2(robbersLargeBN,c("MIsP","MIsR"),c("1","1"))


# 
# 
# structuredCoherence2 <- function(BN, narrationNodes, states){
#   startTime <- proc.time() #start measuring computation time
# 
#   JN <- compile(as.grain(BN))
#   
#   #find non-root nodes
#   parented <- unique(arcs(BN)[,2])
#   
#   parented
#   
#   #assign parents to parented
#   parentList <- list()
#   for(node in 1:length(parented)){
#     parentList[[node]] <- parents(BN,parented[node])  
#   }
#   
#   parentList
#   
#   #initiate lists of results
#   expConfFull <- list()
#   ECSante <- numeric(length(parented))
#   ECSanteS <- numeric(length(parented))
#   ECSanteNarr <- numeric(length(parented))
#   ECSanteNarrS <- numeric(length(parented))
#   
#   
#   parented
#   
#   for (i in 1:length(parented)){
# 
#  
#     
#   
#     consequent <- parented[i]
#     
#     consequent
#     
#     consequentStates <- if (consequent %in% narrationNodes){
#       stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)    
#     } else {
#       findStates(consequent)
#     }
#     
#     
#     antecedents  <- parentList[i]
#     antecedents
#     
#     antecedentStates <- list()
#     for(a in 1:length(antecedents[[1]])){
#       antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
#         stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes, states = states)    
#       } else {
#         findStates(antecedents[[1]][a])
#       }
#     }
#     
#     #start the outcome table
#     variants <-  expand.grid(c(list(consequentStates),antecedentStates))
#     colnames(variants) <- c(consequent,antecedents[[1]])
#     
#     variants
#     
#     #add the prior of consequent, to be used for Z calculation
#     pr <- numeric(nrow(variants))
#     for (s in 1:nrow(variants)){
#         pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
#       }
#     variants <- cbind(variants, priorCons = pr)
#     variants
#     
#     #posteriors 
#     posteriors <- numeric(nrow(variants))
#     for (row in 1:nrow(variants)){
#       JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
#       posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
#     }
#     variants <- cbind(variants,posteriors)
#     
#     variants
#     
#     
#     #priors for the antecedents
#     priorAnte <- numeric(nrow(variants))
#     
#     for (row in 1:nrow(variants)){
#       rowNodes <- as.vector(unlist(c(antecedents)))
#       rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
#       PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
#       PriorJoints <- aperm(PriorJoints, rowNodes)
#       steps <- numeric(length(rowNodes))
#       for(rn in 1:length(rowNodes)){
#         steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
#       }
#       steps<- gsub(" ", "", steps, fixed = TRUE)
#       final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
#       noquote(final)
#       prior <- eval(parse(text=final))
#       priorAnte[row] <- prior
#     }
#     variants$PriorAnte <- priorAnte
#     
#     
#     #variants
#     
#     
#     #priors for the antecedents given narrations
#     priorAnteNarr <- numeric(nrow(variants))
#     
#     for (row in 1:nrow(variants)){
#       
#       
#       
#       rowNodes <- as.vector(unlist(c(antecedents)))
#       
#       rowNodes
#       
#       rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
#       
#       rowStates
#       
#       JNnarr <-  setEvidence(JN, nodes = narrationNodes, states = states)
#       
#      # JNnarr
#       
#       PriorJoints <- querygrain(JNnarr,nodes=rowNodes,type="joint")
#       
#     
#     #  PriorJoints
#       
#     #  is.numeric(PriorJoints)
#       
#       loose <- which(! rowNodes  %in% narrationNodes)
#       
#      # length(loose)
#       
#       PriorJoints[rowStates[loose]]
#       
#       if(length(loose) == 0){
#         priorAnteNarr[row] <- 1
#       } else if   (is.numeric(PriorJoints)){
#         priorAnteNarr[row]  <- PriorJoints[rowStates[loose]]       #PriorJoints[1]
#       } else {  
#       relevantNodes <-   rowNodes[which(rowNodes %in% names(PriorJoints))]
# 
#       relevantStates <- rowStates[which(rowNodes %in% names(PriorJoints))]
#       
#       relevantNodes
#       relevantStates
#       
#       
#       
#       arPriors <- as.array(PriorJoints)
#       
#      
#       
#       PriorJoints <-   aperm(arPriors[[1]], relevantNodes)
#       steps <- numeric(length(relevantNodes))
#       for(rn in 1:length(relevantNodes)){
#         steps[rn] <- paste(relevantNodes[rn], "=", "\"",relevantStates[rn],"\"")
#       }
#       steps<- gsub(" ", "", steps, fixed = TRUE)
#       final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
#       prior <- eval(parse(text=final))
#       priorAnteNarr[row] <- prior
#       }
#       
#       priorAnteNarr
#       
#       
#       
#       length(loose) == 0
#     
#       rowNodes
#       rowStates
#       narrationNodes
#       
#       
#       tight <- which( narrationNodes %in% rowNodes)
#       
#       states[tight]
#       
#         
#       PriorJoints[rowStates[loose]]
#       
#       if(length(loose) == 0){
#           1
#       } else if   (is.numeric(PriorJoints) & length(PriorJoints) >= 1){
#           priorAnteNarr[row] <-       PriorJoints[rowStates[loose]]
#     #PriorJoints[1]
#       } else {
#           PriorJoints <- aperm(PriorJoints, rowNodes)
#           steps <- numeric(length(rowNodes))
#           for(rn in 1:length(rowNodes)){
#               steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
#               }
#           steps<- gsub(" ", "", steps, fixed = TRUE)
#           final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
#           prior <- eval(parse(text=final))
#           priorAnteNarr[row] <- prior
#       }
#     }
#     variants$PriorAnteNarr <- priorAnteNarr
#     
#     variants
#     
#     
#   
#     #scaled weights for antecedent
#     if(sum(variants$PriorAnte) > 0){
#       variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
#     } else {
#       variants$WeightsAnte <-  1/nrow(variants)
#     }
#     
#     
#     sum(variants$PriorAnteNarr, na.rm = TRUE) > 0
#     
#     #scaled weights for antecedent given conjunction
#     if(sum(variants$PriorAnteNarr, na.rm = TRUE) > 0){
#       variants$WeightsNarr <-  variants$PriorAnteNarr / sum(variants$PriorAnteNarr)
#     } else {
#       variants$WeightsNarr <-  1/nrow(variants)
#     }
#     
#     
#     
#     #Z measures
#     variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
#     #variants$Zweighted <- variants$Z * variants$Weights
#     #variants$ZweightedCon <- variants$Z * variants$PriorJoint
#     variants$ZweightedAnte <- variants$Z * variants$PriorAnte
#     #variants$ZscaledCon <- variants$Z * variants$WeightsCon
#     variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
#     variants$ZweightedAnteNarr <- variants$Z * variants$PriorAnteNarr
#     variants$ZweightedAnteNarrS <- variants$Z * variants$WeightsNarr
#     
#     
#     variants
#     
#     
#     expConfFull[[i]] <- list( "Consequent node" = consequent,
#                               "Options & calculations" = variants,
#                               "ECSante" = sum(variants$ZweightedAnte),
#                               "ECSante scaled" = sum(variants$ZscaledAnte),
#                               "ECSanteNarr" = sum(variants$ZweightedAnteNarr),
#                               "ECSanteNarrS" =  sum(variants$ZweightedAnteNarrS)
#                               )
#     ECSante[i] <-  sum(variants$ZweightedAnte)
#     ECSanteS[i]  <-  sum(variants$ZscaledAnte)
#     ECSanteNarr[i] <- sum(variants$ZweightedAnteNarr)
#     ECSanteNarrS[i] <-  sum(variants$ZweightedAnteNarrS)
#   }
#   
#   
#   ECSanteNarr[is.na(ECSanteNarr)] <- -1
#   ECSanteNarrS[is.na(ECSanteNarr)] <- -1
#   
#     
# 
#   populationSD <- function( vector ){
#     sqrt(sum((vector - mean(vector))^2)/(length(vector)))
#   }
#   
#   
#   structuredScore <- function(expConf)  if (min(expConf) <= 0) {
#     (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
#   } else {
#     (mean(expConf) - populationSD(expConf))
#   }
#   
#   stopTime <- proc.time()
#   elapsedTime <- stopTime - startTime
#   
#   return(list("Full calculations" = expConfFull,  
#               "ECSante" = ECSante,
#               "ECSanteS" = ECSanteS,
#               "ECSanteNarr" = ECSanteNarr,
#               "ECSanteNarrS" = ECSanteNarrS, 
#               "structured Coherence Ante" = structuredScore(ECSante),
#               "structured Coherence Ante scaled" = structuredScore(ECSanteS),
#               "structured Coherence Narr" = structuredScore(ECSanteNarr), 
#               "structured Coherence Narr scaled" = structuredScore(ECSanteNarrS),
#               "Computation time" = elapsedTime))
# }
# 
# 
# 
# # PRgreaterPnR <- robbersTable[1,] > robbersTable[2,] 
# # #PRgreaterPnR
# # 
# # PRgreaterNeutral <- robbersTable[1,] > neutralPoints
# # 
# # robbersResults <- as.data.frame(rbind(PRgreaterPnR,PRgreaterNeutral))
# # 
# # rownames(robbersResults) <- c("Robbers: PR$>$P$\\neg$R","Robbers: PR$>$neutral")
# # 
# # 
# # save(robbersResults,file="calculations/RdataObjects/robbersResults.Rda")
# 
# 
# BN <- RegularBN
# narrationNodes <- c("T","TF")
# states <- c("1","1")
# 
# 
# 
# 
# 
structuredCoherenceNarration <- function(BN, narrationNodes, states){
  startTime <- proc.time() #start measuring computation time
  JN <- compile(as.grain(BN))
  JNnarration <- setEvidence(JN,nodes = narrationNodes, states = states)


  #find non-root nodes
  parented <- unique(arcs(BN)[,2])

  parented

  #assign parents to parented
  parentList <- list()
  for(node in 1:length(parented)){
    parentList[[node]] <- bnlearn::parents(BN,parented[node])
  }

  #initiate lists of results
  expConfFull <- list()
  ECSanteS <- numeric(length(parented))
  ECSnarrS <- numeric(length(parented))



  for (i in 1:length(parented)){
    consequent <- parented[i]


    consequentStates <- if (consequent %in% narrationNodes){
      stateOfNode(node = consequent, narrationNodes = narrationNodes, states = states)
    } else {
      findStates(consequent)
    }



    antecedents  <- parentList[i]


    antecedentStates <- list()
    for(a in 1:length(antecedents[[1]])){
      antecedentStates[[a]] <- if (antecedents[[1]][a] %in% narrationNodes){
        stateOfNode(node = antecedents[[1]][a], narrationNodes = narrationNodes,
                    states = states)
      } else if (!(antecedents[[1]][a] %in% parented)){
        eval(parse(text = paste('dimnames(BN$',antecedents[[1]][a],'$prob)',
                                sep= "")))[[1]]
      } else {
        findStates(antecedents[[1]][a])
      }
    }




    #start the outcome table
    variants <-  expand.grid(c(list(consequentStates),antecedentStates))
    colnames(variants) <- c(consequent,antecedents[[1]])




    #add the prior of consequent, to be used for Z calculation
    pr <- numeric(nrow(variants))
    for (s in 1:nrow(variants)){
      pr[s] <- as.numeric(querygrain(JN, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[s]",sep=""))))])
    }
    variants <- cbind(variants, priorCons = pr)

    variants

    #posteriors
    posteriors <- numeric(nrow(variants))
    for (row in 1:nrow(variants)){
      JNtemp <- setEvidence(JN, nodes = antecedents[[1]], states = as.vector(unlist(variants[row,antecedents[[1]],drop= FALSE]))  )
      posteriors[row]   <-   querygrain(JNtemp, nodes = consequent)[[1]][as.character(eval(parse(text = paste("variants$",consequent,"[row]",sep=""))))]
    }
    variants <- cbind(variants,posteriors)


    variants



    priorAnte <- numeric(nrow(variants))

    for (row in 1:nrow(variants)){
      rowNodes <- as.vector(unlist(c(antecedents)))
      rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
      PriorJoints <- querygrain(JN,nodes=rowNodes,type="joint")
      PriorJoints <- aperm(PriorJoints, rowNodes)
      steps <- numeric(length(rowNodes))
      for(rn in 1:length(rowNodes)){
        steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
      }
      steps<- gsub(" ", "", steps, fixed = TRUE)
      final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
      noquote(final)
      prior <- eval(parse(text=final))
      priorAnte[row] <- prior
    }
    variants$PriorAnte <- priorAnte

    variants

    priorAnteNarr <- numeric(nrow(variants))

    for (row in 1:nrow(variants)){
      rowNodes <- as.vector(unlist(c(antecedents)))
      rowStates <- as.vector(unlist(variants[row,2:(length(antecedents[[1]])+1)]))
      PriorJoints <- querygrain(JNnarration,nodes=rowNodes,type="joint")
      PriorJoints
      PriorJoints <- if(length(PriorJoints) == 1){
         PriorJoints
        } else {
           aperm(PriorJoints, rowNodes)
        }
      steps <- numeric(length(rowNodes))
      for(rn in 1:length(rowNodes)){
        steps[rn] <- paste(rowNodes[rn], "=", "\"",rowStates[rn],"\"")
      }
      steps<- gsub(" ", "", steps, fixed = TRUE)
      final <- paste("PriorJoints[",paste(steps,collapse=","),"]",sep="")
      noquote(final)
      prior <- eval(parse(text=final))
      priorAnteNarr[row] <- ifelse(length(PriorJoints) == 1,PriorJoints,prior)
    }
    variants$PriorAnteNarr <- priorAnteNarr

    variants


    #scaled weights for antecedent
    if(sum(variants$PriorAnte) > 0){
      variants$WeightsAnte <-  variants$PriorAnte / sum(variants$PriorAnte)
    } else {
      variants$WeightsAnte <-  1/nrow(variants)
    }



    #scaled weights for narrationstructuredCoherenceUpdated(W12BN,A,c("1","1"))
    if(sum(variants$PriorAnteNarr) > 0){
      variants$WeightsAnteNarr <-  variants$PriorAnteNarr / sum(variants$PriorAnteNarr)
    } else {
      variants$WeightsAnteNarr <-  1/nrow(variants)
    }





    #Z measures
    variants$Z <- Z(posterior = variants$posteriors, prior = variants$priorCons)
    #variants$Zweighted <- variants$Z * variants$Weights
    #variants$ZweightedCon <- variants$Z * variants$PriorJoint
    #variants$ZweightedAnte <- variants$Z * variants$PriorAnte
    #variants$ZscaledCon <- variants$Z * variants$WeightsCon
    variants$ZscaledAnte <- variants$Z * variants$WeightsAnte
    variants$ZscaledNarr <- variants$Z * variants$WeightsAnteNarr




    expConfFull[[i]] <- list( "Consequent node" = consequent,
                              "Options & calculations" = variants,
                              "ECSante scaled" = sum(variants$ZscaledAnte),
                              "ECSnarr scaled" = sum(variants$ZscaledNarr))

    ECSanteS[i] <- sum(variants$ZscaledAnte)
    ECSnarrS[i] <- sum(variants$ZscaledNarr)
  }


  populationSD <- function( vector ){
    sqrt(sum((vector - mean(vector))^2)/(length(vector)))
  }

  structuredScore <- function(expConf)  if (min(expConf) <= 0) {
    (mean(expConf) - populationSD(expConf)) * (min(expConf) +1) - min(expConf)^2
  } else {
    (mean(expConf) - populationSD(expConf))
  }

  
  structuredScoreSquared <- function(ECS)  if (min(ECS) <= 0) {
    (mean(ECS) - populationSD(ECS)^2) * (min(ECS) +1) - min(ECS)^2
  } else {
    (mean(ECS) - populationSD(ECS)^2)
  }
  
  
  structuredScoreNoSD <- function(ECS)  if (min(ECS) <= 0) {
    (mean(ECS)) * (min(ECS) +1) - min(ECS)^2
  } else {
    (mean(ECS))
  }
  
  
  
  stopTime <- proc.time()
  elapsedTime <- stopTime - startTime
  
  return(list("Full calculations" = expConfFull, 
              "ECSanteS" = ECSanteS,
              "ECSnarrS" = ECSnarrS,
              "structuredNarrated" = structuredScore(ECSnarrS),
              "structuredSquared" = structuredScoreSquared(ECSnarrS),
              "structuredNoSD" = structuredScoreNoSD(ECSnarrS),
              "Computation time" = elapsedTime))
# 
# 
#   
# 
#   return(list("Full calculations" = expConfFull,
#               "ECSanteS" = ECSanteS,
#               "ECSnarrS" = ECSnarrS,
#               "structured Coherence Ante scaled" = structuredScore(ECSanteS),
#               "structured Coherence Narr scaled" = structuredScore(ECSnarrS),
#               "Computation time" = elapsedTime))
}
# 
# 
# 
# 
# 
# 



