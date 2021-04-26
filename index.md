---
layout: page
title: Bayesian Networks for the Legal Probabilism SEP entry
output:
  md_document:
    variant: markdown_github
    preserve_yaml: true
---

### Set-up

The examples are given in R code, which we intertwine with additional
explanations, which the reader is welcome to skip if they look familiar.

First, you need to install the relevant R libraries. This is a bit
tricky, because some of them have to be installed through BiocManager.
One way to go is this:

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





### Bayesian Networks: a crashcourse


While Bayes's Theorem is of immense use when it comes to calculating various conditional probabilities, if we're interested in the interaction of multiple hypotheses at various levels and multiple pieces of evidence, calculations quickly become inconvenient, to say the least. Moreover, if such considerations are to be presented to a fact-finder, it is rather unlikely that they would be transparent and easily understood. Luckily, a tool exist to make such tasks easier. Bayesian networks (BNs) can be used for a fairly transparent and computationally more manageable evalation and presentation of the interaction of multiple pieces of evidence and hypotheses. We'll start with a general presentation of BNs, and then go over a few main methods of employing BNs in presentation, aggregation and evaluation of evidence in legal fact-finding.  


A \emph{random variable} (RV) $X$ is a function from the elements of a sample space into $\mathbb{R}$, the set of real numbers. For instance, if our sample space is the set of all potential outcomes of tossing a fair coin four times (each such outcome can be represented as a sequence, for instance  $HHHT$, or $HTHT$), $X$ can be the number of heads among the tosses.  

Given a probability measure $P$, two events $A$ and $B$ are conditionally independent given another event $C$, $I_{P}(A,B\vert C)$, just in case $P(A\wedge B\vert C) = P(A\vert C)P(B \vert C)$. Conditional and unconditional independence don't have to coincide. If you toss twice a coin which is fair with probability $\frac{1}{2}$, and $\frac{3}{4}$ biased towards heads with probability $\frac{1}{2}$, the result of the second toss is not independent of the first one. After all, if the first result is heads, this increases the probability that the coin is biased, and so increases the probability of heads in the second toss. On the other hand, conditionally on knowledge whether the coin is fair, the results are independent. If the coin is fair, the probability of heads in the second toss is $\frac{1}{2}$ and if the coin is biased, it is $\frac{3}{4}$, no matter what the first result was. And in the opposite direction, indepedence can disappear when we condition. Say I have two friends, Alice and Peter, who call me regularly, but they decide to do so independently. Then, whether they call in five minutes is independent. Now, suppose the phone rings. Conditional on the phone ringing, I know that if it isn't Alice, it's Peter, and so the identities of the callers are no longer independent.





Two RVs $X$ and $Y$ are conditionally independent given another RV $Z$, $I_{P}(X,Y\vert Z)$ just in case for any combination of values of these RVs $x,y,z$ it is the case that $I_{P}(X=x \wedge Y=y \vert Z=z)$ (notice: $X,Y$ and $Z$ are RVs, while $x,y$ and $z$ are some particular values they can take). The notion naturally generalizes to sets of RVs.
Often, instead of saying things like $P(X_1 = x_1\wedge Y_5=y_5 \vert Z_3=z_3)$ we'll rather say $P(x_1,y_5\vert z_3)$.  


Now, if we have $n$ RVs, even if we assume for simplicity that they're binary (that is, they can take only one of two values), there are $2^n$ possible combinations of values they could take, and so a direct description of a probability measure for them would require $2^n-1$ numbers. This would be a highly unfeasible method of specifying a probability distribution for a set of random variables.



Moreover, even if we had specified the joint probability distribution for all our combinations of values of Rvs $X, Y, Z$, using it wouldn't be the most efficient way of calculating conditional probabilities or the probability that a certain selected RV takes a certain particular value. For instance, we would have to rely on:

$$ xP(x_1\vert y_1)  =  \frac{P{x_1,y_1}}{P{y_1}} \\ &  = \frac{\sum_{i}P(x_1,y_1,Z=z_i)}{
\sum_{i,j}P(X=x_j,y_1,Z=Z_i)}$$

\noindent in which calculations we'd have to travel through all possible values of $Z$ and $X$ -- this would become even less feasible as the number of RVs and their possible values increase. With 100 binary RVs we'd need 2^{99} terms in the sum in the denominator, so it seems that to be able to calculate a single conditional probability we'd have to elicit quite a few uncoditional ones.



Instead, we start with representing dependencies between RVs in such a set by means of a *directed acyclic graph* (DAG). A DAG is a collection of *nodes* (called also *vertices*) -- think of them as corresponding to the RVs, *directed edges* (also called *arcs*; they  can be thought of as ordered pairs of nodes), such that there is no sequence of nodes $v_0,\dots, v_k$ with edges from $v_i$ to $v_{i+1}$ for $0\leq i\leq k-1$ with $v_0=v_k$.\footnote{Sometimes it is also required that the graph should be connected: that for any two nodes there is an undirected path between them.} A *qualitative BN* (QBN) is a DAG with nodes labeled by RVs. Here's one example of a QBN:



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

<img src="https://rfl-urbaniak.github.io/LegalProbabilismBNs/images/unnamed-chunk-3-1.png" style="display: block; margin: auto;" />

With the intended reading:

-----  ------------------------------
RV     Proposition
-----  ----------------------------
PS     At least one parent smokes.

SH     The subject is a second-hand smoker.

S      The subject smokes.

C      The subject develops cancer.
----   ------------------------------------
