# Implementing zkSTARKs

This document is intended to explain the implementation strategy of zkSTARKs over prime fields. Many of the terms of the original paper (e.g. algebraic routing) will be replaced with meanings what these procedures actually do.

FRI procedure will not be covered by this document as there is a highest quality article publically [available](https://eccc.weizmann.ac.il/report/2017/134/download/).

AIR representation will be explained, but there will be no simple tooling given for generation of the AIR constraints. This is a separate and challenging task to be solved.

Soundness analysis and discussions of DEEP-* variants is also not covered by this document.

## 1. Start from AIR

AIR is a representation of the program (or statement) that is being proven as a set of states of some registry machine during an execution. From here we use $T$ (capital) as a total length of the execution trace over $W$ (capital) registers. Each step of the trace will be labeled as either $t$ or $t_i$ and individual registers during the execution will be labeled as either $w$ or $w_j$. When individual registed at the execution step needs to be addressed it's labeled as a tuple $(t_i, w_j)$.

For a full execution trace two sequential states of all the registers $W$ at the steps $t_i, t_{i+1}$ must satisfy a set of polynomial constraints $P$ (capital). Emphasize: each consequative set of registers $W$ at execution steps $t_i, t_{i+1}$ must satifly ALL individual polynomial constraints $p_k$ from the set $P$.

More informally:
```
for p in P {
    for (t, t_next) in pairs(T) {
        p(W_at_t, W_at_t_next) == 0;
    }
}
```

Once again, designing the constraints set $P$ and what "meaning" each of the registers $W$ has during the execution of the program at each of the execution steps $T$ is not covered in this article. Interested readers can try to implement toolings for efficient AIR representation generation for something more difficult than Fibonacci sequence themselves. One good point to start would be to use FSM (finite state machine) paradigm.

## 2. ARP (algebraic routing and placement)

If reader is familiar with a series of Vitalik's posts about STARKs, at this step we should do what Vitalik has called "placing witness values at the roots of unity". Let's go through this in details.

We start from far - from the problem of reduction of R1CS to QAP that is used in SNARKS, for example in Groth16 proof system. Given a series of constraints in a form $<a, b> = c$, where $a$, $b$ and $c$ are linear combinations (one can threat them as sparse vectors) and $<,>$ operator is elementwise multiplication, one can reduce is to the problem of satisfiability of polynomial equation $A(x)*B(x) = C(x)$ for some set of $x$. Usually this set is chosed to be a radix-2 domain of roots of unity in a prime field of a size $D$, with a generator $\omega$ with $\omega^{D} = 1$, $x_i = \omega^{i}$. It allows to find polynomials $A(x)$, etc. using FFT. Also, if $A(x)*B(x) = C(x)$ holds for all $x$, then it's a multiple of some polynomial $Z(x) = (x-x_1)(x-x_2)...(x-x_{|D|})$, with divisibility being proven by presenting a polynomial $H(x) = (A(x)*B(x) - C(x))/Z(x)$.

Now one can go into ARP that is in some sense close to the logic above:

Predicates and notations
- $T$ and $t_i$ as defined above, same for constraints and register states.
- There is a set of register states $(t_i, w_j)$ that represents a correct execution trace. Obviously these values must satisfy some constraints.
- Discussion of incorporation of boundary conditions is skipped in this document, at least in the iniital version.
- We assume all values used in the problem to be elements of either prime field $F_q$ or a domain $D$ if this field that forms a multipliactive group. For elements of such domain a term "roots if unity" is used interchangably.

Ok, an overview of the problem:
- We have registers, their states and polynomial constraints that must be satisfied at EVERY pair of sequential states in a trace.
- We want to reduce it some way to the form that is good for a FRI.
- We have a powerfull trick of polynomial divisibility as reminded by the SNARKs.
- States of the registers $(t_i, w_j)$ are just field elements.

Now a tricky change of terms:
- We call "witness" not a full set of register states $(t_i, w_j)$, but a POLYNOMIAL $f_{witness}(x)$ over the field $F_q$. Formally $f_{witness}(x) \in F_{q}[x]$.
- We use a set of "masking coefficients" (that we later find to be just elements of a multiplicative subgroup) that are used to *pick a specific $x$ for evaluation of $f_{witness}(x)$*. A set of masking coefficients is labeled as $M_k$ and has a length of $2|W|$ (so we can mask every register from the two adjustent state if we want). In principle each constraints may have individual set of masking coefficients, so there is a $k$ index.

For this explainer we first make a statement, that for a correct witness polynomial $f_{witness}(x)$ and a constraint $p_k$ with masking coefficients $M_k$ (use capital notating ignoring that there may be potentially different masks for different constraints, and label them coefficients as $m_0, ..., m_{|M|-1}$ just for ease of notations) there must hold
$$
p_{k}(f_{witness}(x*m_{k,0}), f_{witness}(x*m_{k,1}), ..., f_{witness}(x*m_{k,|M|-1})) = 0
$$

for every $x$ for which $Q_{k}(x) = 0$, where $Q_{k}(x)$ is constraint specific polynomial.

Tough? Yes...

What we have described above is an ARP instance. Now let's try to make a step from AIR representation to ARP instance.

- Let's chose a multiplicative group of a size $|T|*|W|$ with a generator $\omega$, with $\omega^{|T|*|W|} = 1$.
- Let's find $f_{witness}(x)$ (or define, or guess it) from the requirement that $f_{witness}(\omega^{i*|W| + j}) = (t_i, w_j)$. Notation here is quite abused, but important parts to remember are:
  - $|W|$ is a number of registers
  - $t_i$ is a number of the state in the execution trace
  - $(t_i, w_j)$ is a state of the register $w_j$ at the step $t_i$
  - $j$ is an index in a range $[0, |W|)$
  - an exponent for $\omega$ looks a lot like how one could represent indexing of a two-dimentional matrix with dimentions $|T|*|W|$ in the linear array...
- chose (or define, or guess) set of the masking coefficients $M_k$ to be $1, \omega^{1}, ..., \omega^{2|W|-1}$
- chose (or define, or guess, you know) a vanishing polynomial $Q_k(x)$ to vanish at values $1, \omega^{|W|}, ..., \omega^{(|T|-1)*|W|}$

## 3. ALI

Let's say our intuition is corrent, and let's just stuff all those quesses and requirements into the ARP instance:
- $Q_k(x)$ vanishes at $1, \omega^{|W|}, ..., \omega^{(|T|-1)*|W|}$
- for any of those values of $x$ where $Q_k(x)$ vanishes let's chose particular $x_i = \omega^{i|W|}$, so 
$$
f_{witness}(x_{i}*m_j) = \\
f_{witness}(x_{i}*\omega^{j}) = \\
f_{witness}(\omega^{i|W| + j}) = (t_i, w_j)
$$

A trick that we've just pulled: we guessed, or chosen a witness polynomial form that is basically lookups a register state $(t_i, w_j)$ when evaluated at some point $x \in F_q$. We have chosed masking coefficients and roots of $Q_k(x)$ to just be just right to make a system self-consistent. For readers more proficient in the question an explanation may be the following: we did a polynomial interpolation, and have chosen roots of $Q_k(x)$ to be roots of unity for domain of the size $|T|$, while $f_{witness}$ is defined over the domain of size $|T|*|W|$. Let's also define degree $d$, $deg(f_{witness}) < d$ to be an upper bound on degree of the witness.

Now let's observe another property: if 
$$
p_{k}(f_{witness}(x*m_{k,0}), f_{witness}(x*m_{k,1}), ..., f_{witness}(x*m_{k, |M|-1})) = 0
$$
at every $x$ where $Q_k(x) = 0$, then $Q_k(x)$ divides $p_{k}(f_{witness}(x*m_{k,0}),...)$, so one can calculate a quotient polynomial
$$
g_k(x) = \frac{p_{k}(f_{witness}(x*m_{k, 0}), f_{witness}(x*m_{k, 1}), ..., f_{witness}(x*m_{k, |M|-1}))}{Q_k(x)}
$$

if we define $d_c$ as a highest degree of all the polynomials $P_k$ (constraint polynomials), then $deg(g_k) < d*d_c$ (here $p_k(x)$ is *polynomial over polynomials over x*, do it's degree is larger).

Here should be a security proof, but let's simplify it and put as the following:

If there are polynomials $f(x)$ and $g(x)$ derived as defined above, then they are degree abounded, and their *evaluations* $f: D -> F$ (same for $g(x)$) at the corresponding domains $D$ (such evaluations are also called "interpolants" in original articles) are RS codes for a domain $D$, error rate $\rho$ and field $F$ $RS[F, D, \rho]$ where $\rho$ is determined from the degree $d < \rho |D|$.

## 4. Verification

Better description to follow, but the final goal of verifier is to ensure that $f(x)$ and $g(x)$ are indeed degree $d$ and $d*d_c$ bounded polynomials and also to check that constraints hold at some points of the evaluation domain.

## More tricks

To follow


