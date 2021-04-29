blah
================

### Set-up

The examples are given in R code, which we intertwine with additional explanations, which the reader is welcome to skip if they look familiar.

First, you need to install the relevant R libraries. This is a bit tricky, because some of them have to be installed through BiocManager. One way to go is this:

``` r
install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")
install.packages("BiocManager")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("graph", "Rgraphviz"))

install.packages("https://www.bnlearn.com/releases/bnlearn_latest.tar.gz", repos = NULL, type = "source")
```

Then load the libraries we use (you need to have them installed first):

``` r
library(bnlearn)
library(Rgraphviz)
library(knitr)
library(kableExtra)
```

While Bayes's Theorem is of immense use when it comes to calculating various conditional probabilities, if we're interested in the interaction of multiple hypotheses at various levels and multiple pieces of evidence, calculations quickly become inconvenient, to say the least. Moreover, if such considerations are to be presented to a fact-finder, it is rather unlikely that they would be transparent and easily understood. Luckily, a tool exist to make such tasks easier. Bayesian networks (BNs) can be used for a fairly transparent and computationally more manageable evalation and presentation of the interaction of multiple pieces of evidence and hypotheses. We'll start with a general presentation of BNs, and then go over a few main methods of employing BNs in presentation, aggregation and evaluation of evidence in legal fact-finding.

``` r
cancer1 <- empty.graph(nodes = c("PS","SH","S","C"))
cancer1.arcs <- matrix(c("PS", "SH",
                   "PS", "S",
                   "SH", "C",
                    "S", "C"),
                 byrow = TRUE, ncol = 2,
                 dimnames = list(NULL, c("from", "to")))
arcs(cancer1) = cancer1.arcs
graphviz.plot(cancer1)
```

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/cancer-1.png" style="display: block; margin: auto;" />

This can be achieved in a somewhat simpler manner like this:

``` r
cancer2DAG <- model2network("[SH|PS][PS][S|PS][C|SH:S]")
graphviz.plot(cancer2DAG)
```

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/cancerSimpler-1.png" style="display: block; margin: auto;" />

Adding probabilties to the DAG we already have can be achieved easily using a bunch of wrappers we wrote for root nodes, and single- and double-parented nodes. First, load the script, then, remember that `E` stands for the child, and `H` for a parent.

``` r
source("cptCreate.R")

#tables for separate nodes
PSprob <- priorCPT("PS",prob1 = .3) 

Sprob <-  singleCPT(eNode = "S", hNode = "PS", probEifHS1 = .4 , probEifHS2 =  .2)
SHprob <- singleCPT(eNode = "SH", hNode = "PS", probEifHS1 = .8, probEifHS2 = .3)

Cprob <- doubleCPT(eNode= "C", h1Node = "S", h2Node = "SH",
                   probEifH1S1H2S1 = .6,
                   probEifH1S1H2S2 = .4,
                   probEifH1S2H2S1 = .1,
                   probEifH1S2H2S2 = .01)

#put them together and add to the DAG to create a BN
cancerCPT <- list(PS = PSprob, S = Sprob, SH = SHprob, C = Cprob)
cancerBN <- custom.fit(cancer2DAG,cancerCPT)

#display the CPTs for all the nodes
cancerBN
```

    ## 
    ##   Bayesian network parameters
    ## 
    ##   Parameters of node C (multinomial distribution)
    ## 
    ## Conditional probability table:
    ##  
    ## , , SH = 1
    ## 
    ##    S
    ## C      1    0
    ##   1 0.60 0.10
    ##   0 0.40 0.90
    ## 
    ## , , SH = 0
    ## 
    ##    S
    ## C      1    0
    ##   1 0.40 0.01
    ##   0 0.60 0.99
    ## 
    ## 
    ##   Parameters of node PS (multinomial distribution)
    ## 
    ## Conditional probability table:
    ##  PS
    ##   1   0 
    ## 0.3 0.7 
    ## 
    ##   Parameters of node S (multinomial distribution)
    ## 
    ## Conditional probability table:
    ##  
    ##    PS
    ## S     1   0
    ##   1 0.4 0.2
    ##   0 0.6 0.8
    ## 
    ##   Parameters of node SH (multinomial distribution)
    ## 
    ## Conditional probability table:
    ##  
    ##    PS
    ## SH    1   0
    ##   1 0.8 0.3
    ##   0 0.2 0.7

``` r
library(DiagrammeR)
#devtools::install_github("rich-iannone/DiagrammeR")
#devtools::install_github("rich-iannone/DiagrammeRsvg")
library(DiagrammeRsvg)
library(magrittr)
library(rsvg)

graphviz.chart(cancerBN, grid = FALSE, type = "barprob", layout = "neato", scale = c(1,1.5),  main="marginal probabilities for the cancer BN") 
```

    ## Loading required namespace: gRain

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/cancerBarchart-1.png" width="100%" style="display: block; margin: auto;" />

A quantitative BN (further on, simply BN) is a DAG with CPTs, and we say that it *represents* a probabilistic measure $ $ quantitatively just in case they are compatible, and $ $ agrees with its assignment of CPTs. It can be shown that any quantitative BN represents a unique probabilistic measure. However, any probabilistic measure can be represented by multiple BNs.

Here's a sketch of the argument for the last claim above, if you're interested. Skip this passage if you aren't. Take any permutation of the RVs under consideration, obtaining ![X\_1, \\dots, X\_k](https://latex.codecogs.com/png.latex?X_1%2C%20%5Cdots%2C%20X_k "X_1, \dots, X_k"). For each ![i](https://latex.codecogs.com/png.latex?i "i") find a minimal subset ![P\_i](https://latex.codecogs.com/png.latex?P_i "P_i") such that ![I\_{\\mathsf{P} {}}(\\{X\_1,\\dots,x\_{i-1}\\},X\_i\\vert P\_i)](https://latex.codecogs.com/png.latex?I_%7B%5Cmathsf%7BP%7D%20%7B%7D%7D%28%5C%7BX_1%2C%5Cdots%2Cx_%7Bi-1%7D%5C%7D%2CX_i%5Cvert%20P_i%29 "I_{\mathsf{P} {}}(\{X_1,\dots,x_{i-1}\},X_i\vert P_i)") --- that is, a minimal set of RVs which are earlier in the squence, which makes ![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i") independent of all the (other) RVs which are earlier in the squence. Such a set always exists, in the worst-case scenario, it is the set of all ![X\_1,\\dots,X\_{i-1}](https://latex.codecogs.com/png.latex?X_1%2C%5Cdots%2CX_%7Bi-1%7D "X_1,\dots,X_{i-1}"). Next, make ![P\_i](https://latex.codecogs.com/png.latex?P_i "P_i") the parents of ![X\_i](https://latex.codecogs.com/png.latex?X_i "X_i") in the DAG and copy the values of $ $ into the CPTs.

One reccuring pattern captures the relation between a hypothesis and a piece of evidence, the idea being that it is whether the hypothesis is true that caused the (non-)occurence of the evidence.

``` r
HE <- model2network("[Evidence|Hypothesis][Hypothesis]")
graphviz.plot(HE)
```

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/EHbn-1.png" width="100%" style="display: block; margin: auto;" />

For instance, ![H](https://latex.codecogs.com/png.latex?H "H") might take values from the range of ![1-40](https://latex.codecogs.com/png.latex?1-40 "1-40"), the distance in meters from which the gun has been shot, and ![E](https://latex.codecogs.com/png.latex?E "E") might be a continuous variable representing the density of gun shot residues. (This example also indicates that RVs don't have to be binary for legal applications.)

Another example, this time with binary variables, takes ![H](https://latex.codecogs.com/png.latex?H "H") to be the claim that the supect is guilty and ![E](https://latex.codecogs.com/png.latex?E "E") the presence of a DNA match with a crime scene stain. One way the CPTs could look like in this case is this:

``` r
Hypothesis.prob <-  array(c(0.01, 0.99), dim = 2, dimnames = list(HP =  c("guilty","not guilty")))

Evidence.prob <- array(c( 0.99999, 0.00001, 0, 1), dim = c(2,2),dimnames = list(Evidence = c("DNA match","no match"),
      Hypothesis = c("guilty","not guilty")))

HEcpt <- list(Hypothesis=Hypothesis.prob,Evidence=Evidence.prob)
HE_bn = custom.fit(HE,HEcpt)
Evidence.frame <- as.data.frame(Evidence.prob)
row.names(Evidence.frame) <- c("DNA match","No match")
kable(Evidence.frame,col.names = c("guilty","not guilty"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
guilty
</th>
<th style="text-align:right;">
not guilty
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
DNA match
</td>
<td style="text-align:right;">
0.99999
</td>
<td style="text-align:right;">
0
</td>
</tr>
<tr>
<td style="text-align:left;">
No match
</td>
<td style="text-align:right;">
0.00001
</td>
<td style="text-align:right;">
1
</td>
</tr>
</tbody>
</table>
The true power of BNs, however, appears when we go beyond a simple two-node situations for which calculations can be done by hand. For instance, imagine we have two pieces of evidence: a DNA match, and a witness testimony. The DAG might look like this:

``` r
HEE.dag <- model2network("[H][W|H][DNA|H]")
graphviz.plot(HEE.dag)
```

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/unnamed-chunk-5-1.png" style="display: block; margin: auto;" />

The CPTs can, for instance, as follows:

``` r
HEEdag <- model2network("[H][W|H][DNA|H]")
Hprob <- array(c(0.01, 0.99), dim = 2, 
                dimnames = list(h = c("murder","nomurder")))

Wprob <- array(c( 0.7, 0.3, 0.4, 0.6), dim = c(2,2),dimnames = list(W= c("seen","notseen"), H = c("murder","nomurder")))

DNAprob <- array(c( 1, 0, 0.001, 0.999), dim = c(2,2),
                  dimnames = list(DNA =c("dnamatch","nomatch"),
                                  H = c("murder","nomurder")))



HEEcpt <- list(H=Hprob,W=Wprob,DNA = DNAprob)

HEEbn <- custom.fit(HEEdag,HEEcpt)

kable(as.data.frame(Hprob),col.names="Pr(H)")
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
Pr(H)
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
murder
</td>
<td style="text-align:right;">
0.01
</td>
</tr>
<tr>
<td style="text-align:left;">
nomurder
</td>
<td style="text-align:right;">
0.99
</td>
</tr>
</tbody>
</table>
``` r
W.frame <- as.data.frame(Wprob)
row.names(W.frame) <- c("W=seen","W=notseen")
kable(W.frame,col.names = c("H=murder","H=nomurder"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
H=murder
</th>
<th style="text-align:right;">
H=nomurder
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
W=seen
</td>
<td style="text-align:right;">
0.7
</td>
<td style="text-align:right;">
0.4
</td>
</tr>
<tr>
<td style="text-align:left;">
W=notseen
</td>
<td style="text-align:right;">
0.3
</td>
<td style="text-align:right;">
0.6
</td>
</tr>
</tbody>
</table>
``` r
DNA.frame <- as.data.frame(DNAprob)
row.names(DNA.frame) <- c("DNA=match","DNA=nomatch")
kable(DNA.frame,col.names = c("H=murder","H=nomurder"))
```

<table>
<thead>
<tr>
<th style="text-align:left;">
</th>
<th style="text-align:right;">
H=murder
</th>
<th style="text-align:right;">
H=nomurder
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
DNA=match
</td>
<td style="text-align:right;">
1
</td>
<td style="text-align:right;">
0.001
</td>
</tr>
<tr>
<td style="text-align:left;">
DNA=nomatch
</td>
<td style="text-align:right;">
0
</td>
<td style="text-align:right;">
0.999
</td>
</tr>
</tbody>
</table>
The CPT for the hypothesis contains the prior probability that a murder has been commited by the suspect. The CPTs for the other variables include (made up) probabilities of a DNA match and of a witness seeing the suspect near the crime scene at an appropriate time conditional on various states of the murder hypothesis: ![\\mathsf{P} (\\textrm{W=seen}\\vert \\textrm{H=murder})=0.7, \\mathsf{P} (\\textrm{W=seen}\\vert \\textrm{H=no.murder})=0.4](https://latex.codecogs.com/png.latex?%5Cmathsf%7BP%7D%20%28%5Ctextrm%7BW%3Dseen%7D%5Cvert%20%5Ctextrm%7BH%3Dmurder%7D%29%3D0.7%2C%20%5Cmathsf%7BP%7D%20%28%5Ctextrm%7BW%3Dseen%7D%5Cvert%20%5Ctextrm%7BH%3Dno.murder%7D%29%3D0.4 "\mathsf{P} (\textrm{W=seen}\vert \textrm{H=murder})=0.7, \mathsf{P} (\textrm{W=seen}\vert \textrm{H=no.murder})=0.4") etc.

``` r
library(bnlearn)
library(gRain)
junction <- compile(as.grain(HEEbn))

junctionMS <- setEvidence(junction, nodes = c("DNA","W"), states = c("dnamatch","seen") )
querygrain(junctionMS)$H
```

    ## H
    ##     murder   nomurder 
    ## 0.94645754 0.05354246

``` r
junctionMN <- setEvidence(junction, nodes = c("DNA","W"), states = c("dnamatch","notseen"))
querygrain(junctionMN)$H
```

    ## H
    ##    murder  nomurder 
    ## 0.8347245 0.1652755

``` r
junctionNOMATCH <- setEvidence(junction, nodes = c("DNA"), states = c("nomatch"))
querygrain(junctionNOMATCH)$H
```

    ## H
    ##   murder nomurder 
    ##        0        1

``` r
HEEms <- as.bn.fit(junctionMS, including.evidence = TRUE) 
HEEmn <- as.bn.fit(junctionMN, including.evidence = TRUE) 
HEEnomatch <- as.bn.fit(junctionNOMATCH, including.evidence = TRUE) 

graphviz.chart(HEEms, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  witness evidence")

graphviz.chart(HEEmn, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  negative witness evidence")

graphviz.chart(HEEnomatch, grid = FALSE, type = "barprob",  scale = c(2,2), 
               main="marginal probabilities after DNA match and  no witness evidence")
```

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/unnamed-chunk-8-1.png" width="100%" style="display: block; margin: auto;" /><img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/unnamed-chunk-8-2.png" width="100%" style="display: block; margin: auto;" /><img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/unnamed-chunk-8-3.png" width="100%" style="display: block; margin: auto;" />
