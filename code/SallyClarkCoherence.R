install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", 
                 repos = NULL, type = "source")
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "Rgraphviz"))
BiocManager::install(c("graph", "RBGL", "Rgraphviz"))
install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", 
                 repos = NULL, type = "source")


install.packages("gRbase", dependencies=TRUE); 
install.packages("gRain", dependencies=TRUE); 
install.packages("gRim", dependencies=TRUE)
library(bnlearn)
library(Rgraphviz)
library(gRain)

getwd()
source("utils//CombinationsBN.R")
source("utils//CptCreate.R")
source("utils//LogicAndBNs.R")
source("utils//kableCPTs.R")
source("measures//structuredCoherence.R")
source("utils//CoherenceTables.R")
source("measures//Fitelson.R")
source("measures//DouvenMeijs.R")
source("measures//generalizedOlsson.R")
source("measures//Olsson.R")
source("measures//generalizedShogenji.R")
source("measures//Shogenji.R")
source("measures//Roche.R")
#source("measures//RA.R")
Z <- function(posterior,prior){
  d <- posterior - prior
  ifelse(prior == posterior, 0, ifelse(posterior > prior, d/(1-prior), d/prior))
}

source("cptCreate.R")


#define the structure of the Sally Clark BN
SallyClarkDAG <- model2network("[Abruising|Acause][Adisease|Acause][Bbruising|Bcause][Bdisease|Bcause][Acause][Bcause|Acause][NoMurdered|Acause:Bcause][Guilty|NoMurdered]")

#plot 
graphviz.plot(SallyClarkDAG)


#CPTs as used in Fenton & al.
AcauseProb <-prior.CPT("Acause","SIDS","Murder",0.921659)
AbruisingProb <- single.CPT("Abruising","Acause","Yes","No","SIDS","Murder",0.01,0.05)
AdiseaseProb <- single.CPT("Adisease","Acause","Yes","No","SIDS","Murder",0.05,0.001)
BbruisingProb <- single.CPT("Bbruising","Bcause","Yes","No","SIDS","Murder",0.01,0.05)
BdiseaseProb <- single.CPT("Bdisease","Bcause","Yes","No","SIDS","Murder",0.05,0.001)
BcauseProb <- single.CPT("Bcause","Acause","SIDS","Murder","SIDS","Murder",0.9993604,1-0.9998538)

#E goes first; order: last variable through levels, second last, then first
NoMurderedProb <- array(c(0, 0, 1, 0, 1, 0, 0,1,0,1,0,0), dim = c(3, 2, 2),dimnames = list(NoMurdered = c("both","one","none"),Bcause = c("SIDS","Murder"), Acause = c("SIDS","Murder")))

#this one is definitional
GuiltyProb <-  array(c( 1,0, 1,0, 0,1), dim = c(2,3),dimnames = list(Guilty = c("Yes","No"), NoMurdered = c("both","one","none")))

# Put CPTs together
SallyClarkCPT <- list(Acause=AcauseProb,Adisease = AdiseaseProb,
                      Bcause = BcauseProb,Bdisease=BdiseaseProb,
                      Abruising = AbruisingProb,Bbruising = BbruisingProb,
                      NoMurdered = NoMurderedProb,Guilty=GuiltyProb)

# join with the DAG to get a BN
SallyClarkBN <- custom.fit(SallyClarkDAG,SallyClarkCPT)

graphviz.chart(SallyClarkBN,type="barprob", scale = c(0.7,1.3), main = "Priors in the Sally Clark case")



SCnodes <- c("Acause","Bcause")

SCstates <- list(c("SIDS","SIDS"), c("Murder", "Murder"), c("SIDS", "Murder"), c("Murder", "SIDS"))

SCtable <- CoherencesTable(SallyClarkBN, scenariosList = SCnodes,
          statesList = c("SIDS","SIDS"), exampleName = "Sally Clark")





penguinsTable <- rbind(penguinsTableBGP,penguinsTableBG,penguinsTableBP)


