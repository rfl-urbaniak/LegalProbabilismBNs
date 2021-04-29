
# Hypothesis - evidence
HE <- model2network("[Evidence|Hypothesis][Hypothesis]")

graphviz.plot(HE)

Hypothesis.prob <-  array(c(0.01, 0.99), dim = 2, dimnames = list(HP =  c("guilty","not guilty")))

Evidence.prob <- array(c( 0.99999, 0.00001, 0, 1), dim = c(2,2),dimnames = list(Evidence = c("DNA match","no match"), Hypothesis = c("guilty","not guilty")))

HEcpt <- list(Hypothesis=Hypothesis.prob,Evidence=Evidence.prob)
HEbn = custom.fit(HE,HEcpt)


# Hypothesis - two pieces of evidence

HEEdag <- model2network("[H][W|H][DNA|H]")
graphviz.plot(HEE.dag)

Hprob <- array(c(0.01, 0.99), dim = 2, 
                dimnames = list(h = c("murder","nomurder")))

Wprob <- array(c( 0.7, 0.3, 0.4, 0.6), dim = c(2,2),dimnames = list(W= c("seen","notseen"), H = c("murder","nomurder")))

DNAprob <- array(c( 1, 0, 0.001, 0.999), dim = c(2,2),
                  dimnames = list(DNA =c("dnamatch","nomatch"),
                                  H = c("murder","nomurder")))



HEEcpt <- list(H=Hprob,W=Wprob,DNA = DNAprob)

HEEbn <- custom.fit(HEEdag,HEEcpt)

# Calculating probabilities in HEE with exact inference 
library(gRain)



#convert to a junction tree for calculations
junction <- compile(as.grain(HEEbn))

#update with match and seen
junctionMS <- setEvidence(junction, nodes = c("DNA","W"), states = c("dnamatch","seen") )
querygrain(junctionMS)$H


#update with match but not seen
junctionMN <- setEvidence(junction, nodes = c("DNA","W"), states = c("dnamatch","notseen"))
querygrain(junction.mn)$H



#update with no match
junctionNOMATCH <- setEvidence(junction, nodes = c("DNA"), states = c("nomatch"))
querygrain(junctionNOMATCH)$H


#convert to a BNs with propagation
HEEms <- as.bn.fit(junctionMS, including.evidence = TRUE) 
HEEmn <- as.bn.fit(junctionMN, including.evidence = TRUE) 
HEEnomatch <- as.bn.fit(junctionNOMATCH, including.evidence = TRUE) 



#plot updated BNs
graphviz.chart(HEEms, grid = FALSE, type = "barprob",  scale = c(1.5,1.5), 
               main="DNA match and  witness evidence")



graphviz.chart(HEEmn, grid = FALSE, type = "barprob",  scale = c(1.5,1.5), 
               main="DNA match and negative witness evidence")



graphviz.chart(HEEnomatch, grid = FALSE, type = "barprob",  scale = c(1.5,1.5), 
               main="DNA match and  no witness evidence")


#cause - consequence

CauseCon <-  model2network("[Cause][Consequence|Cause]")
graphviz.plot(CauseCon)

#measurement
Measurement <- model2network("[Accuracy][Actual value][Observed value|Accuracy:Actual value]")
graphviz.plot(Measurement)

#synthesis
Definitional <- model2network("[Distance][Time][Velocity|Distance:Time]")
graphviz.plot(Definitional)

#evidence accuracy
evidenceAccuracy <-  model2network("[Accuracy of evidence][Excess alcohol level][Evidence for excess|Accuracy of evidence:Excess alcohol level]")
graphviz.plot(evidenceAccuracy)


#opportunity

opportunity <- model2network("[H1][H2|H1][A1][A2][E1|H1:A1][E2|H1:A2]")
graphviz.plot(opportunity)


#dependency
cameras <- model2network("[H][D][C1|H][C2|H:C1:D]")
graphviz.plot(cameras, layout = "neato")


#alibi
alibi <- model2network("[S present][S guilty|S present][Alibi accuracy|S guilty][Alibi|Alibi accuracy:S present]")
graphviz.plot(alibi)


#mixed DNA
DNA789 <- model2network("[S is C1][S is C2][Genotype of C1|S is C1][Genotype of C2|S is C2][S is the source|S is C1:S is C2][(7,8,9) found|Genotype of C1:Genotype of C2]")
graphviz.plot(DNA789)
