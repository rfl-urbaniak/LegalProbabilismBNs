---
layout: page
title: blah
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

While Bayes’s Theorem is of immense use when it comes to calculating
various conditional probabilities, if we’re interested in the
interaction of multiple hypotheses at various levels and multiple pieces
of evidence, calculations quickly become inconvenient, to say the least.
Moreover, if such considerations are to be presented to a fact-finder,
it is rather unlikely that they would be transparent and easily
understood. Luckily, a tool exist to make such tasks easier. Bayesian
networks (BNs) can be used for a fairly transparent and computationally
more manageable evalation and presentation of the interaction of
multiple pieces of evidence and hypotheses. We’ll start with a general
presentation of BNs, and then go over a few main methods of employing
BNs in presentation, aggregation and evaluation of evidence in legal
fact-finding.

A (RV) *X* is a function from the elements of a sample space into ℝ, the
set of real numbers. For instance, if our sample space is the set of all
potential outcomes of tossing a fair coin four times (each such outcome
can be represented as a sequence, for instance *H**H**H**T*, or
*H**T**H**T*), *X* can be the number of heads among the tosses.

Given a probability measure *P*, two events *A* and *B* are
conditionally independent given another event *C*,
*I*<sub>*P*</sub>(*A*, *B*\|*C*), just in case
*P*(*A* ∧ *B*\|*C*) = *P*(*A*\|*C*)*P*(*B*\|*C*). Conditional and
unconditional independence don’t have to coincide. If you toss twice a
coin which is fair with probability $\\frac{1}{2}$, and $\\frac{3}{4}$
biased towards heads with probability $\\frac{1}{2}$, the result of the
second toss is not independent of the first one. After all, if the first
result is heads, this increases the probability that the coin is biased,
and so increases the probability of heads in the second toss. On the
other hand, conditionally on knowledge whether the coin is fair, the
results are independent. If the coin is fair, the probability of heads
in the second toss is $\\frac{1}{2}$ and if the coin is biased, it is
$\\frac{3}{4}$, no matter what the first result was. And in the opposite
direction, indepedence can disappear when we condition. Say I have two
friends, Alice and Peter, who call me regularly, but they decide to do
so independently. Then, whether they call in five minutes is
independent. Now, suppose the phone rings. Conditional on the phone
ringing, I know that if it isn’t Alice, it’s Peter, and so the
identities of the callers are no longer independent.

Two RVs *X* and *Y* are conditonally independent given another RV *Z*,
*I*<sub>*P*</sub>(*X*, *Y*\|*Z*) just in case for any combination of
values of these RVs *x*, *y*, *z* it is the case that
*I*<sub>*P*</sub>(*X* = *x* ∧ *Y* = *y*\|*Z* = *z*) (notice: *X*, *Y*
and *Z* are RVs, while *x*, *y* and *z* are some particular values they
can take). The notion naturally generalizes to sets of RVs. Often,
instead of saying things like
*P*(*X*<sub>1</sub> = *x*<sub>1</sub> ∧ *Y*<sub>5</sub> = *y*<sub>5</sub>\|*Z*<sub>3</sub> = *z*<sub>3</sub>)
we’ll rather say *P*(*x*<sub>1</sub>, *y*<sub>5</sub>\|*z*<sub>3</sub>).

Now, if we have *n* RVs, even if we assume for simplicity that they’re
binary (that is, they can take only one of two values), there are
2<sup>*n*</sup> possible combinations of values they could take, and so
a direct description of a probability measure for them would require
2<sup>*n*</sup> − 1 numbers. This would be a highly unfeasible method of
specifying a probability distribution for a set of random variables.

Moreover, even if we had specified the joint probability distribution
for all our combinations of values of Rvs *X*, *Y*, *Z*, using it
wouldn’t be the most efficient way of calculating conditional
probabilities or the probability that a certain selected RV takes a
certain particular value. For instance, we would have to rely on: in
which calculations we’d have to travel through all possible values of
*Z* and *X* – this would become even less feasible as the number of RVs
and their possible values increase. With 100 binary RVs we’d need 2^{99}
terms in the sum in the denominator, so it seems that to be able to
calculate a single conditional probability we’d have to elicit quite a
few uncoditional ones.

Instead, we start with representing dependencies between RVs in such a
set by means of a (DAG). A DAG is a collection of (called also ) – think
of them as corresponding to the RVs, (also called ; they can be thought
of as ordered pairs of nodes), such that there is no sequence of nodes
*v*<sub>0</sub>, …, *v*<sub>*k*</sub> with edges from *v*<sub>*i*</sub>
to *v*<sub>*i* + 1</sub> for 0 ≤ *i* ≤ *k* − 1 with
*v*<sub>0</sub> = *v*<sub>*k*</sub>. A (QBN) is a DAG with nodes labeled
by RVs. Here’s one example of a QBN:

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

Looking at our examples, one possible combination of CPTs is:

In contrast, BN representation significantly simplifies the information
that needs to be stored and calculated. If RVs are binary and each node
has at most *k* parents, at most 2<sup>*k*</sup>*n* values are needed to
determine the distribution.

Continuing our example, we now only need to assign prior probabilities
to parentless nodes and *conditional probability tables* (CPTs) to
edges. If *X* has *m* possible states and is the parent of *Y* with *n*
possible states, the CPTs for their connection is an *m* × *n* table
containing conditional probabilities for $\\pr(Y= y\\vert X=x)$. In
general, once we have defined CPTs for each node given its parents, the
product of these conditional probabilities yields a joint probability
distribution satisfying the Markov condition, if we’re dealing with
discrete RVs (the claim also holds for normally distributed RVs, but not
universally for any type of continuous RVs).

Whence simplicity?
------------------

This is a somewhat more technical explanation of how BNs help in
reducing complexity. An uninterested reader can skip ahead.

The tells us that $\\pr(A\\wedge B) = \\pr(A\\vert B)\\pr(B)$. Its
application to RVs (say the RVs in G are
*X*<sub>1</sub>, …*X*<sub>*n*</sub>) yields:

So, if $\\pr$ is compatible with G, we don’t have to represent it
directly by listing all the 2<sup>*n*</sup> − 1 values. Instead, the
joint probability $\\pr(X\_1\\dots,X\_n)$ (note: this is really an
assignment of probability values to possible combinations of the values
of these RVs), can be represented using the conditional probabilities on
the right-hand side of , and moreover, for each conditional probability
of an RVs *X* given some other RVs, non-parents of *X* can be removed
from the condition, since RVs are independent of them. Let’s slow down
and take a look at the argument. Pick an ancestral ordering
*X*<sub>1</sub>, …, *X*<sub>*n*</sub>, of the RVs, that is, an ordering
in which if *Z* is a descendant of *X*, *Z* follows *Y* in the ordering.
Take any selection of values of these variables,
*x*<sub>1</sub>, …, *x*<sub>*n*</sub>. Let pa<sub>i</sub> be the set
cotaining all the values of the parents of *X*<sub>*i*</sub> that belong
to this sequence. Since this is an ancestral ordering, the parents have
to occur before *x*<sub>*i*</sub>. We need to prove We prove it by
induction on the length of the sequence. The basic case comes for free.
Now assume: we need to show: One option is that
$\\pr(x\_i,x\_{i-1},\\dots,x)=0$. Then, also
$\\pr(x\_{i+1}, x\_{i}, \\dots, x\_1)=0$, and by the induction
hypothesis, there is a 1 ≤ *k* ≤ *i* such that
$\\pr(x\_k\\vert \\mathsf{pa\_k})=0$, and so also the right-hand side of
equals 0 and so holds.

Another option is that $\\pr(x\_i,x\_{i-1},\\dots,x)\\neq 0$. Then we
reason:

The first step is by the chain rule. The second is by the Markov
condition and the fact that we employed an ancestral ordering. The third
one uses . This ends the proof.

Now, why does the product of CPTs yield a joint probability distribution
satisfying the Markov condition, if we’re dealing with discrete RVs?
