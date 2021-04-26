---
layout: page
#title: Bayesian Networks for the Legal Probabilism SEP entry
output:
  md_document:
    variant: markdown_github
    preserve_yaml: true
---

The examples are given in R code, which we intertwine with additional
explanations, which the reader is welcome to skip if they look familiar.


### Set-up


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


A *random variable* (RV) $X$ is a function from the elements of a sample space into $\mathbb{R}$, the set of real numbers. For instance, if our sample space is the set of all potential outcomes of tossing a fair coin four times (each such outcome can be represented as a sequence, for instance  $HHHT$, or $HTHT$), $X$ can be the number of heads among the tosses.  

Given a probability measure $P$, two events $A$ and $B$ are conditionally independent given another event $C$, $I_{P}(A,B\vert C)$, just in case $P(A\wedge B\vert C) = P(A\vert C)P(B \vert C)$. Conditional and unconditional independence don't have to coincide. If you toss twice a coin which is fair with probability $\frac{1}{2}$, and $\frac{3}{4}$ biased towards heads with probability $\frac{1}{2}$, the result of the second toss is not independent of the first one. After all, if the first result is heads, this increases the probability that the coin is biased, and so increases the probability of heads in the second toss. On the other hand, conditionally on knowledge whether the coin is fair, the results are independent. If the coin is fair, the probability of heads in the second toss is $\frac{1}{2}$ and if the coin is biased, it is $\frac{3}{4}$, no matter what the first result was. And in the opposite direction, indepedence can disappear when we condition. Say I have two friends, Alice and Peter, who call me regularly, but they decide to do so independently. Then, whether they call in five minutes is independent. Now, suppose the phone rings. Conditional on the phone ringing, I know that if it isn't Alice, it's Peter, and so the identities of the callers are no longer independent.





Two RVs $X$ and $Y$ are conditionally independent given another RV $Z$, $I_{P}(X,Y\vert Z)$ just in case for any combination of values of these RVs $x,y,z$ it is the case that $I_{P}(X=x \wedge Y=y \vert Z=z)$ (notice: $X,Y$ and $Z$ are RVs, while $x,y$ and $z$ are some particular values they can take). The notion naturally generalizes to sets of RVs.
Often, instead of saying things like $P(X_1 = x_1\wedge Y_5=y_5 \vert Z_3=z_3)$ we'll rather say $P(x_1,y_5\vert z_3)$.  


Now, if we have $n$ RVs, even if we assume for simplicity that they're binary (that is, they can take only one of two values), there are $2^n$ possible combinations of values they could take, and so a direct description of a probability measure for them would require $2^n-1$ numbers. This would be a highly unfeasible method of specifying a probability distribution for a set of random variables.



Moreover, even if we had specified the joint probability distribution for all our combinations of values of Rvs $X, Y, Z$, using it wouldn't be the most efficient way of calculating conditional probabilities or the probability that a certain selected RV takes a certain particular value. For instance, we would have to rely on:

$$xP(x_1\vert y_1)  =  \frac{P{x_1,y_1}}{P{y_1}}   = \frac{\sum_{i}P(x_1,y_1,Z=z_i)}{
\sum_{i,j}P(X=x_j,y_1,Z=Z_i)}$$

 in which calculations we'd have to travel through all possible values of $Z$ and $X$ -- this would become even less feasible as the number of RVs and their possible values increase. With 100 binary RVs we'd need 2^{99} terms in the sum in the denominator, so it seems that to be able to calculate a single conditional probability we'd have to elicit quite a few uncoditional ones.



Instead, we start with representing dependencies between RVs in such a set by means of a *directed acyclic graph* (DAG). A DAG is a collection of *nodes* (called also *vertices*) -- think of them as corresponding to the RVs, *directed edges* (also called *arcs*; they  can be thought of as ordered pairs of nodes), such that there is no sequence of nodes $v_0,\dots, v_k$ with edges from $v_i$ to $v_{i+1}$ for $0\leq i\leq k-1$ with $v_0=v_k$. (Sometimes it is also required that the graph should be connected: that for any two nodes there is an undirected path between them.) A *qualitative BN* (QBN) is a DAG with nodes labeled by RVs.

### An example with R

Here's one example of a QBN:



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

<table style="width:50%;">
<colgroup>
<col style="width: 9%" />
<col style="width: 40%" />
</colgroup>
<thead>
<tr class="header">
<th style="text-align: left;">RV</th>
<th style="text-align: left;">Proposition</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td style="text-align: left;">PS</td>
<td style="text-align: left;">At least one parent smokes.</td>
</tr>
<tr class="even">
<td style="text-align: left;">SH</td>
<td style="text-align: left;">The subject is a second-hand smoker.</td>
</tr>
<tr class="odd">
<td style="text-align: left;">S</td>
<td style="text-align: left;">The subject smokes.</td>
</tr>
<tr class="even">
<td style="text-align: left;">C</td>
<td style="text-align: left;">The subject develops cancer.</td>
</tr>
</tbody>
</table>


 The *ancestors* of $C$ are all the other RVs, the  parents of $C$ are $SH$ and $S$, the descentants of $PS$ are all the other RVs, and the children of $PS$ are $SH$ and $S$.

The edges, intuitively, are meant to capture direct influence between RVs. The role that such direct influence plays is that in a BN any RV is conditionally independent of its nondescentants (including ancestors), given its parents. If this is the case for a given probabilistic measure $P{}$ and a given DAG $\mathsf{G}$, we say that $P{}$ is compatible with $\mathsf{G}$, or that it satisfies the *Markov condition*.







## Whence simplicity?

This is a somewhat more technical explanation of how BNs help in reducing complexity. An uninterested reader can skip it.

First, the *chain rule* tells us that $P(A\wedge B) = P(A\vert B)P(B)$. Its application to RVs (say the RVs in G are $X_1,\dots X_n$) yields:


$$
P(x_1,\cdots,x_n)  = P(x_n\vert x_1,\cdots,x_{n-1})\times
\\ = P(x_{n-1}\vert x_1,\cdots,x_{n-2})\times \\  \times \cdots \times P(x_2 \vert x_1) P(x_1)
$$

So, if $P$ is compatible with $\mathsf{G}$, we don't have to represent it directly by listing all the $2^n-1$ values. Instead, the joint probability $P(X_1\dots,X_n)$ (note: this is really an assignment of probability values to *all* possible combinations of the values of these RVs), can be represented using the conditional probabilities on the right-hand side, and moreover, for each conditional probability of an RVs $X$ given some other RVs, non-parents of $X$ can be removed from the condition, since RVs are independent of them.
Let's slow down and take a look at the argument. Pick an ancestral ordering $X_1, \dots, X_n$, of the RVs, that is, an ordering in which if $Z$ is a descendant of $X$, $Z$ follows $Y$ in the ordering. Take any selection of values of these variables, $x_1, \dots, x_n$. Let $\mathsf{pa_i}$ be the set containing all the values of the parents of $X_i$ that belong to this sequence. Since this is an ancestral ordering, the parents have to occur before $x_i$. We need to prove

$$
P(x_n, x_{n-1}, \dots, x_1) = P(x_n\vert \mathsf{pa_n})P(x_{n-1}\vert \mathsf{pa_{n-1}})\cdots P(x_1)
$$

We prove it by induction on the length of the sequence. The basic case comes for free. Now assume:

$$
P(x_i, x_{i-1}, \dots, x_1 = P(x_i\vert \mathsf{pa_i})P(x_{i-1}\vert \mathsf{pa_{i-1}})\cdots P(x_1)$$

 we need to show:

$$P(x_{i+1}, x_{i}, \dots, x_1) = P(x_{i+1}\vert \mathsf{pa_{i+1}})P(x_{i}\vert \mathsf{pa_{i}})\cdots P(x_1)\label{eq:markoproof3}$$

One option is that $P(x_i,x_{i-1},\dots,x)=0$. Then, also  $P(x_{i+1}, x_{i}, \dots, x_1)=0$, and by the induction hypothesis, there is a $1\leq k\leq i$ such that $P(x_k\vert \mathsf{pa_k})=0$, and so also the right-hand side   equals 0 and so the claim holds.

Another option is that  $P(x_i,x_{i-1},\dots,x)\neq 0$.
Then we reason:

$$
P(x_{i+1}, x_{i}, \dots, x_1) = P(x_{i+1}\vert x_i, \dots, x_1) P(x_i, \dots, x_1)\\
 = P(x_{i+1}\vert \mathsf{pa_{i+1}}) P(x_i, \dots, x_1)\\
 = P(x_{i+1}\vert \mathsf{pa_{i+1}})  P(x_{i}\vert \mathsf{pa_{i}})\cdots P(x_1)
 $$


 The first step is by the chain rule. The second is by the Markov  condition and the fact that we employed an ancestral ordering. The third one uses an identity we already obtained. This ends the proof.


Now, why  does the product of CPTs yield a joint probability distribution satisfying the Markov condition, if we're dealing with discrete RVs?


The argument for the discrete case is as follows. Fix an ancestral ordering of the RVs and definite  joint probability in terms of the conditional probabilities determined by the CPTs as we already did. We need to show it indeed is a probabilistic measure. Since the conditional probabilities are between 0 and 1, clearly so is their product. Then we need to show that the sum of $P(x_1,x_2,\dots,x_n)$ for all combinations of values of $X_1,X_2,\dots, X_n$ is 1. This can be shown by moving summation symbols around, the calculations start as follows:

$$
\sum_{x_1,\dots,x_n}P(x_1,\dots, x_n)  =  \\  = \sum_{x_1,\dots,x_n} P(x_n\vert \mathsf{pa_n})P(x_{n-1}\vert \mathsf{pa_{n-1}})\cdots P(x_1)\\
= \sum_{x_1,\dots,x_{n-1}}\left[\sum_{x_n}P(x_n\vert \mathsf{pa_n})\right]
P(x_{n-1}\vert \mathsf{pa_{n-1}})\cdots P(x_1)\\
= \sum_{x_1,\dots,x_{n-1}}P(x_1,\dots, x_{n-1})\left[1\right] \dots
P(x_{n-1}\vert \mathsf{pa_{n-1}})\cdots P(x_1)\\
$$

 when we continue this way, all the factors reduce to 1.


To show that the Markov condition holds in the resulting joint distribution we have to show that for any $x_k$ we have $P(x_k\vert \mathsf{nd_k})= P(x_k \vert \mathsf{pa_k})$, where $nd_k$ is any combination of values for all the nondescendants of $X_k$, $\mathsf{ND_k}$. So take any $k$, order the nodes so that all and only the nondescendants of $X_k$ precede it. Let $\widehat{x}_k,\widehat{nd}_k$ be some particular values of $X_k$ and $\mathsf{ND_k}$. $\mathsf{d_k}$ ranges over combinations of values of the descendants of $X_K$. The reasoning goes:

$$
P(\widehat{x}_k\vert \widehat{nd}_k)
= \frac{P(\widehat{x}_k,\widehat{nd}_k)}{P(\widehat{nd}_k)}\\
= \frac{\sum_{\mathsf{d_k}}P(\widehat{x}_1,
  \dots,\widehat{x}_k,x_{k+1},\dots, x_n)}{
\sum_{\mathsf{d_k}\cup\{x_k\}}P(\widehat{x}_1,
  \dots,\widehat{x}_{k-1},x_{k},\dots, x_n)}
\\
= \frac{\sum_{\mathsf{d_k}}P(x_n\vert \mathsf{pa_n})\cdots
P(x_{k+1}\vert \mathsf{pa_{k+1}})P(\widehat{x}_{k}\vert \mathsf{\widehat{pa}_{k}})\cdots P(\widehat{x}_1)
}{
\sum_{\mathsf{d_k}\cup\{x_k\}}P(x_n\vert \mathsf{pa_n})\cdots
P(x_{k}\vert \mathsf{pa_{k}})P(\widehat{x}_{k-1}
  \vert \mathsf{\widehat{pa}_{k-1}})\cdots P(\widehat{x}_1)
} \\
= \frac{P(\widehat{x}_{k}
  \vert \mathsf{\widehat{pa}_{k}})
  \cdots P(\widehat{x}_1)\sum_{\mathsf{d_k}}P(x_n\vert \mathsf{pa_n})\cdots
P(x_{k+1}\vert \mathsf{pa_{k+1}})
}{
P(\widehat{x}_{k-1}
  \vert \mathsf{\widehat{pa}_{k-1}})\cdots P(\widehat{x}_1)\sum_{\mathsf{d_k}\cup\{x_k\}}P(x_n\vert \mathsf{pa_n})\cdots
P(x_{k}\vert \mathsf{pa_{k}})
}
$$

Now, the sums in the fractions sum both to 1 (for reasons clear from the previous step of the proof), and

 $$
 P(\widehat{x}_{k-1}
  \vert \mathsf{\widehat{pa}_{k-1}}
  )\cdots P(\widehat{x}_1)
  $$

are both in the numerator and the denominator, so we are left with

$$
P(\widehat{x}_{k}
    \vert \mathsf{\widehat{pa}_{k}})
$$

as desired.
