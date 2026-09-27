# What the programs check

This note is for readers who want to audit the finite part of the paper. Start with `example.sage` for an ordinary calculation. Run `verify.sage` for the entire parameter domain.

## Parameters and data

A potential counterexample is described by the field size $q$, redundancy $\rho$, generalized-weight order $t$, and covering radius $r$. The paper proves that only certain combinations require computation. The programs independently enumerate those combinations. The records in `proof_data.zip` contain proposed lengths, weight bounds, and puncturing sequences. They are **proposals to be checked**, not codes claimed to exist.

When first run, `verify.sage` extracts the archive into `certificates/`. This folder and the generated `outputs/` folder are ignored by Git. The archive contains 87 original input files, organized by the four computational parts of the proof.

## The two ways to exclude a tuple

For each necessary tuple the program checks a test length $N\ge\rho+t$ and proves that any putative counterexample has $n\ge N$. Then it checks either:

1. **Length:** $N$ exceeds the incidence/weighted-line upper bound $M(q,\rho,t,r)$ for a counterexample.
2. **Weight:** every $[N,N-\rho]_q$ code satisfies $d_t\le2r+2$. The paper's shortening lemma transfers this to the original length $n\ge N$.

A lower bound $n\ge N$ follows from the elementary bounds stated in the paper or from a strict exact ball-covering inequality. The programs re-evaluate that inequality and recompute $M$.

## Weight proof graphs

A graph node asserts $d_j\le b$ for every code with given length, field, and redundancy. Each node uses one of five rules: generalized Singleton, an ordinary Hamming/Griesmer distance bound, a syndrome-class ball inequality, fixed-radius iteration, or puncturing a subcode with minimum support. The puncturing rule checks **every possible integer support size** between the known bounds. A parent node may use only earlier checked children. Roots must cover every necessary tuple, and the independent audit checks graph reachability.

For example, the record $(q,\rho,t,r)=(2,15,3,6)$ has $N=68$ and $M=73$. Since $68\le73$, the length bounds are inconclusive. Its proof graph establishes $d_3\le14=2r+2$ for every binary $[68,53]$ code. Run `sage example.sage` to see the covering and length computations and the root node number.

## Puncturing sequences and tails

A numerical sequence lists pairs $(e_i,s_i)$. Starting at $(S,W)=(0,0)$, where $S$ is accumulated nullity and $W$ is the number of removed coordinates, a step checks a strict ball inequality at the residual parameters and updates $S\leftarrow S+s_i$ and $W\leftarrow W+e_i(s_i+1)$. It succeeds only if $S=t$ and $W\le2r+2$. Symbolic sequences check the paper's uniform inequalities for whole parameter ranges. Rational tail records check the value and nonnegative slope of affine inequalities at an initial radius, hence cover all larger radii.

## Detailed programs and coverage

| Program | Input inside `proof_data.zip` | Coverage |
| --- | --- | ---: |
| `verify_range.py` | `redundancy/rho_*.json` | 59,086 tuples; 537,416 graph nodes |
| `verify_fixed_gaps.py` | `fixed_gaps/gap_*.json` | 2,396,360 tuples; 11,803 nodes |
| `verify_extensions.py` | `extensions/order_*.json`, `symbolic.json` | 619,126 tuples; 1,022 symbolic paths |
| `verify_new.py` | `new_finite.json`, `tails.json` | 65,360 tuples; 454 paths; 28 rational tails |
| `independent_audit.py` | all relevant files | Independent count, inequality, chain, tail, and graph checks |

The ranges overlap. The checkers generate expected tuples independently and reject missing or duplicate input records. `independent_audit.py` imports none of the other programs and uses a different binomial-coefficient computation for balls. All proof comparisons use integers or rational numbers, never floating point.

`verify.sage` executes each program using SageMath's Python interpreter. If a check fails, it stops and gives the corresponding log filename in `outputs/`. The wrapper deliberately checks arithmetic only: source labels in the article must be reviewed separately after the final TeX reorganization.
