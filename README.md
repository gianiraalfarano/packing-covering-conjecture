# Computations for the generalized packing–covering conjecture

This repository accompanies the paper *A proof of the generalized packing–covering conjecture* by Gianira N. Alfarano, Giuseppe Marino, Alessandro Neri and Rocco Trombetti. It checks the finite computations used in the proof of $d_t(\mathcal C)\le 2R_t(\mathcal C)+2$.

## Run

With SageMath 10 or later installed, open a terminal in this folder and type:

```sh
sage verify.sage
```

The program unpacks `proof_data.zip` (using ZIP LZMA compression) into `certificates/`, checks every required case using exact arithmetic, and finishes with `ALL COMPUTATIONAL CHECKS PASSED` if all checks succeed. Logs are written to `outputs/`. Both folders are created by the run. No linear codes are enumerated.

To inspect the data files without running the check, unpack the archive with 7-Zip or with `sage -python -c "import zipfile; zipfile.ZipFile('proof_data.zip').extractall()"`. Some standard unzip tools cannot open LZMA-compressed archives.

The data file contains proposed bounds and puncturing sequences. The programs enumerate the parameter cases required by the paper themselves. Each case is either excluded by a direct inequality or matched to the necessary input record; every bound and step used from a record is checked with exact arithmetic.

Tested with SageMath 10.8 (passagemath 10.8.12); the full check took about 1.5 minutes.

For a small, readable example, type `sage example.sage`. It examines the case $(q,\rho,t,r)=(2,15,3,6)$; it does not replace the full check.

## Files

| File | Purpose |
| --- | --- |
| `verify.sage` | Runs all checks with one command. |
| `example.sage` | Works through one finite case, $(q,\rho,t,r)=(2,15,3,6)$, in SageMath. |
| `proof_data.zip` | Proposed test lengths, weight bounds and puncturing sequences. `verify.sage` unpacks it into `certificates/`: `redundancy/`, `fixed_gaps/`, `extensions/`, `new_finite.json`, and `tails.json`. |
| `code/verify_range.py` | Redundancy $\rho\le50$. |
| `code/verify_fixed_gaps.py` | Radius gaps $1\le r-t\le5$ (orders $3\le t\le131$). |
| `code/verify_extensions.py` | The 1,022 sequences for orders $t\ge32$ with $6\le r-t\le9$, and the ranges $r\le2t$ or $r-t\le15$ for orders $3\le t\le31$. |
| `code/verify_new.py` | The large-radius bounds for $3\le t\le31$ and the remaining finite set. |
| `code/independent_audit.py` | A second, separately written check of the parameter counts, sequences and large-radius bounds; it uses none of the other programs. |
| `docs/verification.md` | Mathematical explanation of the records and of each check. |
