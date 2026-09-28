"""Run all exact computational checks from the repository directory.

Usage: sage verify.sage
"""

from sage.all import ZZ, QQ
from pathlib import Path
from zipfile import ZipFile
import subprocess
import sys

assert ZZ(2)**10 == 1024 and QQ(1)/3 + QQ(2)/3 == 1
if not Path("proof_data.zip").is_file():
    raise SystemExit("Run this command from the repository directory.")

# Always unpack, so that the checked data are exactly those in proof_data.zip
# (a stale certificates/ folder from an older archive is overwritten).
with ZipFile("proof_data.zip") as archive:
    archive.extractall(".")

out = Path("outputs")
out.mkdir(exist_ok=True)
checks = [
    ("redundancy", "verify_range.py", "--start", "1", "--end", "50",
     "--directory", "certificates/redundancy"),
    ("first five gaps", "verify_fixed_gaps.py", "--start", "1", "--end", "5",
     "--directory", "certificates/fixed_gaps"),
    ("high-order boundary", "verify_extensions.py", "--directory",
     "certificates/extensions"),
    ("low-order completion", "verify_new.py"),
    ("independent audit", "independent_audit.py"),
]
for label, program, *arguments in checks:
    suffix = ".json" if program == "independent_audit.py" else ".txt"
    log = out / (program.removesuffix(".py") + suffix)
    with log.open("w") as file:
        result = subprocess.run([sys.executable, "code/" + program, *arguments],
                                stdout=file, stderr=subprocess.STDOUT)
    if result.returncode:
        raise SystemExit(f"FAILED: {label}. See {log}")
    print("PASS:", label, flush=True)

print("ALL COMPUTATIONAL CHECKS PASSED")
