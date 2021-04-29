
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
graphviz.chart(HEEms, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  witness evidence")



graphviz.chart(HEEmn, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  negative witness evidence")



graphviz.chart(HEEnomatch, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  no witness evidence")


