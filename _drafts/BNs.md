---
layout: page
output:
  md_document:
    #variant: #markdown_github
    preserve_yaml: true
---

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

Given a probability measure $\\pr$, two events *A* and *B* are
conditionally independent given another event *C*,
$I\_{\\pr}(A,B\\vert C)$, just in case
$\\pr(A\\et B\\vert C) = \\pr(A\\vert C)\\pr(B \\vert C)$. Conditional
and unconditional independence don’t have to coincide. If you toss twice
a coin which is fair with probability $\\nicefrac{1}{2}$, and
$\\nicefrac{3}{4}$ biased towards heads with probability
$\\nicefrac{1}{2}$, the result of the second toss is not independent of
the first one. After all, if the first result is heads, this increases
the probability that the coin is biased, and so increases the
probability of heads in the second toss. On the other hand,
conditionally on knowledge whether the coin is fair, the results are
independent. If the coin is fair, the probability of heads in the second
toss is $\\nicefrac{1}{2}$ and if the coin is biased, it is
$\\nicefrac{3}{4}$, no matter what the first result was. And in the
opposite direction, indepedence can disappear when we condition. Say I
have two friends, Alice and Peter, who call me regularly, but they
decide to do so independently. Then, whether they call in five minutes
is independent. Now, suppose the phone rings. Conditional on the phone
ringing, I know that if it isn’t Alice, it’s Peter, and so the
identities of the callers are no longer independent.

Two RVs *X* and *Y* are conditonally independent given another RV *Z*,
$I\_{\\pr}(X,Y\\vert Z)$ just in case for any combination of values of
these RVs *x*, *y*, *z* it is the case that
$I\_{\\pr}(X=x \\et Y=y \\vert Z=z)$ (notice: *X*, *Y* and *Z* are RVs,
while *x*, *y* and *z* are some particular values they can take). The
notion naturally generalizes to sets of RVs. Often, instead of saying
things like $\\pr(X\_1 = x\_1\\et Y\_5=y\_5 \\vert Z\_3=z\_3)$ we’ll
rather say $\\pr(x\_1,y\_5\\vert z\_3)$.

Now, if we have *n* RVs, even if we assume for simplicity that they’re
binary (that is, they can take only one of two values), there are
2<sup>*n*</sup> possible combinations of values they could take, and so
a direct description of a probability measure for them would require
2<sup>*n*</sup> − 1 numbers. This would be a highly unfeasible method of
specifying a probability distribution for a set of random variables.

Optional material starts

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

Optional material ends

Instead, we start with representing dependencies between RVs in such a
set by means of a (DAG). A DAG is a collection of (called also ) – think
of them as corresponding to the RVs, (also called ; they can be thought
of as ordered pairs of nodes), such that there is no sequence of nodes
*v*<sub>0</sub>, …, *v*<sub>*k*</sub> with edges from *v*<sub>*i*</sub>
to *v*<sub>*i* + 1</sub> for 0 ≤ *i* ≤ *k* − 1 with
*v*<sub>0</sub> = *v*<sub>*k*</sub>. A (QBN) is a DAG with nodes labeled
by RVs. Here’s one example of a QBN:
