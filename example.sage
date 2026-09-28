"""A small SageMath walkthrough of one finite root record.

Run from the repository directory: sage example.sage
It follows the example (q, rho, t, r) = (2, 16, 3, 6) of Appendix A.5 of
the paper step by step. It illustrates the inequalities only; verify.sage
checks every record.
"""

from sage.all import ZZ, binomial
from pathlib import Path
from zipfile import ZipFile
import json

if not Path("certificates/redundancy/rho_016.json").is_file():
    with ZipFile("proof_data.zip") as archive:
        archive.extractall(".")

data = json.loads(Path("certificates/redundancy/rho_016.json").read_text())
root = next(z for z in data["roots"] if z["parameters"] == [2, 16, 3, 6])
q, rho, t, r = map(ZZ, root["parameters"])
N, cap = ZZ(root["lower"]), ZZ(root["cap"])


def ball(alphabet, length, radius):
    return sum(binomial(length, i) * (alphabet - 1)**i
               for i in range(radius + 1))


def griesmer(q, k, d):
    return sum((d + q**i - 1) // q**i for i in range(k))


def rejects(q, n, h, d):
    """True if no [n, n-h]_q code has minimum distance d (Hamming,
    punctured Hamming for even binary d, Griesmer)."""
    if ball(q, n, (d - 1) // 2) > q**h:
        return True
    if q == 2 and d % 2 == 0 and ball(2, n - 1, (d - 2) // 2) > 2**(h - 1):
        return True
    return griesmer(q, n - h, d) > n


# Step 1. Length bounds. The exact covering bound gives n >= 112, while the
# length upper bound is M = 137, so the comparison is inconclusive.
n_cover = next(n for n in range(rho + t, 10**4)
               if ball(q**t, n, r) >= q**(t * rho))
m = 2*r - t + 2
incidence = [x + ((2*r + 1 - x) * (q**(rho - x) - 1)) // (q**(m - x) - 1)
             for x in range(m)]
dim = rho - m + 2
alpha, beta = divmod(t + 1, q + 1)
weighted = alpha * sum(q**j for j in range(dim))
if beta:
    weighted += 1 + (beta - 1) * sum(q**j for j in range(dim - 1))
M = min(incidence + [m - 2 + weighted])
assert n_cover == 112 and M == cap == 137

# Step 2. The record works at the (valid, weaker) test length N = 96.
assert N == 96 and ball(q**t, N - 1, r) < q**(t * rho)

# Step 3. Lemma 3.8 with s = 2 and radius 3: d_2 <= 9 for every [96,80]_2 code.
assert ball(q, N, 3) == 147537 > q**(rho + 2 - 1)
d2 = 3 * 3

# Step 4. For each possible value z of d_2, puncture a minimum-support
# 2-dimensional subcode: the residual code is [96-z, 78]_2, with redundancy
# 18-z. The ordinary bounds give an upper bound on its minimum distance.
bounds = []
for z in range(2, d2 + 1):
    n_res, h_res = N - z, rho - z + 2
    d1 = max(d for d in range(1, h_res + 2) if not rejects(q, n_res, h_res, d))
    bounds.append(d1)
    assert z + d1 <= 13
assert bounds == [6, 6, 6, 6, 5, 4, 4, 4]

# Step 5. Lemma 3.9 gives d_3 <= 13 <= 2r+2 at length 96; Lemma 3.11
# (shortening) transfers the bound to every length n >= 96.
assert root["kind"] == "weight" and root["bound"] == 13 <= 2*r + 2

print("Parameters (q, rho, t, r):", (q, rho, t, r))
print("Covering bound n >=", n_cover, "; length upper bound M =", M)
print("Test length N =", N, "; d_2 <=", d2)
print("d_1 bounds of the residual codes for z = 2..9:", bounds)
print("Hence d_3 <= 13 <= 2r+2 = 14 for every [96,80]_2 code.")
