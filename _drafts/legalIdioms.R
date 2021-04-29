
# Hypothesis - evidence
HE <- model2network("[Evidence|Hypothesis][Hypothesis]")

graphviz.plot(HE)

Hypothesis.prob <-  array(c(0.01, 0.99), dim = 2, dimnames = list(HP =  c("guilty","not guilty")))

Evidence.prob <- array(c( 0.99999, 0.00001, 0, 1), dim = c(2,2),dimnames = list(Evidence = c("DNA match","no match"), Hypothesis = c("guilty","not guilty")))

HEcpt <- list(Hypothesis=Hypothesis.prob,Evidence=Evidence.prob)
HEbn = custom.fit(HE,HEcpt)


# Hypothesis - two pieces of evidence

HEE.dag <- model2network("[H][W|H][DNA|H]")
graphviz.plot(HEE.dag)

H.prob <- array(c(0.01, 0.99), dim = 2, 
                dimnames = list(h = c("murder","no.murder")))

W.prob <- array(c( 0.7, 0.3, 0.4, 0.6), dim = c(2,2),dimnames = list(W= c("seen","not.seen"), H = c("murder","no.murder")))

DNA.prob <- array(c( 1, 0, 0.001, 0.999), dim = c(2,2),
                  dimnames = list(DNA =c("dna.match","no.match"),
                                  H = c("murder","no.murder")))

HEE.cpt <- list(H=H.prob,W=W.prob,DNA = DNA.prob)

HEEbn <- custom.fit(HEE.dag,HEE.cpt)
