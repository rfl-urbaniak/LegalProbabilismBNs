# Load all libraries 
library(ggplot2)
library(gridExtra)
library(bnlearn)
library(knitr)
library(kableExtra)
library(gRain)
library(reshape2)
library(tidyverse)
library(plyr)
library(rje)
library(bnlearn)
library(utils)
library(latex2exp)
library(useful)
library(tidyverse)
library(stringr)
library(plot3D)





#________ Load all utility function
source("utils//CombinationsBN.R")
source("utils//CptCreate.R")
source("utils//LogicAndBNs.R")
source("utils//kableCPTs.R")


# #_______ Load confirmation measures
source("measures//CoherenceMeasures2Var.R")
source("measures//DouvenMeijs.R")
source("measures//Fitelson.R")
source("measures//Olsson.R")
source("measures//RA.R")
source("measures//Roche.R")
source("measures//Shogenji.R")
source("measures//generalizedShogenji.R")
source("measures//generalizedOlsson.R")
source("measures//structuredCoherence.R")
# 
# 
# #______ Load the BNs used as examples
source("bns//Books.R")
source("bns//Depth.R")
source("bns//Dodecahedron.R")
source("bns//Dunnit.R")
source("bns//Penguins.R")
source("bns//Robbers.R")
source("bns//Witness.R")
source("bns//JapaneseSwords.R")
source("bns//Beatles.R")

# #___Load the calculations
# 
source("utils//CoherenceTables.R")
# 
#  source("calculations/BeatlesCoherence.R")
# source("calculations/BooksCoherence.R")
# source("calculations/DepthCoherence.R")
# source("calculations/DodecahedronCoherence.R")
# source("calculations/DunnitCoherence.R")
# source("calculations/JapaneseSwordsCoherence.R")
# source("calculations/PenguinsCoherence.R")
# source("calculations/RobbersCoherence.R")
# source("calculations/WitnessCoherence.R")
# # 

#source("calculations/AllCoherenceTables.R")

