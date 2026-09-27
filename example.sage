"""A small SageMath walkthrough of one finite root record.

Run from the repository directory: sage example.sage
This illustrates the inequalities; verify.sage checks every record.
"""

from sage.all import ZZ, binomial
from pathlib import Path
from zipfile import ZipFile
import json

if not Path("certificates/redundancy/rho_015.json").is_file():
    with ZipFile("proof_data.zip") as archive:
        archive.extractall(".")

data = json.loads(Path("certificates/redundancy/rho_015.json").read_text())
root = next(z for z in data["roots"] if z["parameters"] == [2, 15, 3, 6])
q, rho, t, r = map(ZZ, root["parameters"])
N, cap = ZZ(root["lower"]), ZZ(root["cap"])

def ball(alphabet, length, radius):
    return sum(binomial(length, i) * (alphabet - 1)**i
               for i in range(radius + 1))

# The covering inequality implies that a counterexample has n >= N.
assert ball(q**t, N - 1, r) < q**(t * rho)

# Recompute the incidence and weighted-line upper bounds on n.
m = 2*r - t + 2
incidence = [x + ((2*r + 1 - x) * (q**(rho - x) - 1))
             // (q**(m - x) - 1) for x in range(m)]
dim = rho - m + 2
alpha, beta = divmod(t + 1, q + 1)
weighted = alpha * sum(q**j for j in range(dim))
if beta:
    weighted += 1 + (beta - 1) * sum(q**j for j in range(dim - 1))
assert min(incidence + [m - 2 + weighted]) == cap

# N <= cap, so this root uses a generalized-weight proof graph.
assert N == 68 and cap == 73 and root["kind"] == "weight"
node = data["nodes"][root["node"]]
assert node["parameters"] == [2, 68, 15, 3]
assert node["value"] <= root["bound"] == 2*r + 2 == 14

print("Parameters (q, rho, t, r):", (q, rho, t, r))
print("Length lower bound and upper cap:", N, cap)
print("Weight proof node:", root["node"], node["kind"])
print("The full checker checks every child of that proof node.")
