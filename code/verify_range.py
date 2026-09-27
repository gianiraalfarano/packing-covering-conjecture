#!/usr/bin/env python3
"""Independent checker for proof DAGs produced by prove_range.py.

This checker does NOT import the search program, call its dynamic program, or
trust its reported upper bounds. It checks every local inference, every support
size in a residual inference, every length inequality, and completeness of the
finite parameter enumeration. All comparisons use integers.

Usage: python verify_range.py --end 50 --directory certificates
"""
from __future__ import annotations
import argparse
import json
from math import isqrt
from pathlib import Path
from typing import Any


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def ball(alphabet: int, length: int, radius: int) -> int:
    require(length >= 0 and radius >= 0, 'Invalid ball parameters')
    total = term = 1
    for j in range(1, min(length, radius) + 1):
        numerator = term * (length - j + 1) * (alphabet - 1)
        require(numerator % j == 0, 'Nonintegral ball recurrence')
        term = numerator // j
        total += term
    return total


def is_prime_power(q: int) -> bool:
    if q < 2:
        return False
    p = next((d for d in range(2, isqrt(q) + 1) if q % d == 0), q)
    z = q
    while z % p == 0:
        z //= p
    return z == 1


def expected_parameters(h: int) -> set[tuple[int, int, int, int]]:
    ans = set()
    for t in range(1, h + 1):
        for r in range(t, h + 1):
            if t <= 2 or r == t or t >= h - 4:
                continue
            if 3 * r >= 2 * h or 2 * r + 3 > h + t:
                continue
            if (t, r) in {(3, 4), (3, 5), (4, 5), (4, 6)}:
                continue
            for q in range(2, r):
                if is_prime_power(q):
                    ans.add((q, h, t, r))
    return ans


def griesmer(q: int, k: int, distance: int) -> int:
    denominator, index, answer = 1, 0, 0
    while index < k and denominator < distance:
        answer += -(-distance // denominator)
        denominator *= q
        index += 1
    return answer + k - index


def rejects_ordinary(q: int, n: int, h: int, d: int) -> bool:
    if ball(q, n, (d - 1) // 2) > q ** h:
        return True
    if q == 2 and d % 2 == 0:
        if h < 1 or ball(2, n - 1, (d - 2) // 2) > 2 ** (h - 1):
            return True
    return griesmer(q, n - h, d) > n


def verify_nodes(nodes: list[dict[str, Any]]) -> None:
    seen_parameters = set()
    for i, node in enumerate(nodes):
        p = tuple(node['parameters'])
        require(len(p) == 4 and p not in seen_parameters, f'Duplicate node {i}')
        seen_parameters.add(p)
        q, n, h, t = p
        require(is_prime_power(q) and h >= 0 and 1 <= t <= n - h,
                f'Invalid node parameters {p}')
        value, kind = node['value'], node['kind']
        require(t <= value <= h + t, f'Invalid claimed bound in node {i}')
        if kind == 'singleton':
            require(value == h + t, f'Singleton bound in node {i}')
        elif kind == 'ordinary':
            require(t == 1, f'Ordinary rule at higher order, node {i}')
            for d in range(value + 1, h + 2):
                require(rejects_ordinary(q, n, h, d),
                        f'Unexcluded ordinary distance {d} in node {i}')
        elif kind == 'collision':
            e = node['radius']
            require(e >= 0 and value >= (t + 1) * e, f'Collision size, node {i}')
            require(ball(q, n, e) > q ** (h + t - 1), f'Collision count, node {i}')
        elif kind == 'iteration':
            e, L, A = node['radius'], node['stages'], node['initial_order']
            require(e >= 1 and L >= 1 and A >= 1, f'Invalid iteration parameters, node {i}')
            S = 0
            for stage in range(L):
                S = min(t, S + A + (e - 1) * S + e * stage)
            B = e * (t + L)
            require(S >= t and value >= B and n > B, f'Iteration support size, node {i}')
            require(ball(q, n - B, e) > q ** (h + A - 1), f'Iteration count, node {i}')
        elif kind == 'residual':
            s = node['s']
            support_id = node['support_bound']
            require(1 <= s < t and 0 <= support_id < i, f'Residual order, node {i}')
            support = nodes[support_id]
            require(tuple(support['parameters']) == (q, n, h, s),
                    f'Wrong support-bound state, node {i}')
            jmax = support['value']
            children = node['children']
            js = [entry[0] for entry in children]
            require(bool(js) and len(js) == len(set(js)) and max(js) >= jmax
                    and set(js) == set(range(s, max(js) + 1)),
                    f'Incomplete support-size interval, node {i}')
            for j, child_id in children:
                require(0 <= child_id < i, f'Non-topological edge, node {i}')
                child = nodes[child_id]
                require(tuple(child['parameters']) == (q, n - j, h - j + s, t - s),
                        f'Wrong residual parameters, node {i}')
                require(value >= j + child['value'], f'Failed lifting bound, node {i}')
        else:
            raise ValueError(f'Unknown node type: {kind}')


def check_length_cap(q: int, h: int, t: int, r: int, claimed: int) -> None:
    m = 2 * r - t + 2
    require(2 <= m < h, 'Failure indices are not admissible')
    # Average over m-flats through x independent indexed columns.
    caps = []
    for x in range(m):
        numerator = (2 * r + 1 - x) * (q ** (h - x) - 1)
        denominator = q ** (m - x) - 1
        caps.append(x + numerator // denominator)
    # Independent-seed weighted projective-line bound, including zero images.
    dim = h - m + 2
    alpha = (t + 1) // (q + 1)
    beta = (t + 1) - alpha * (q + 1)
    points = sum(q ** j for j in range(dim))
    pencil = sum(q ** j for j in range(dim - 1))
    u = alpha * points
    if beta:
        u += 1 + (beta - 1) * pencil
    caps.append(m - 2 + u)
    require(claimed == min(caps), 'Incorrect length cap')


def verify_layer(path: Path, h: int) -> dict[str, int]:
    data = json.loads(path.read_text())
    require(data['format'] == 'packing-covering-proof-dag-v1', 'Unknown format')
    require(data['redundancy'] == h, 'Wrong redundancy label')
    nodes, roots = data['nodes'], data['roots']
    verify_nodes(nodes)
    actual = [tuple(root['parameters']) for root in roots]
    require(len(actual) == len(set(actual)), 'Duplicated parameter tuple')
    require(set(actual) == expected_parameters(h), 'Finite parameter enumeration is incomplete')
    count = dict(length=0, weight=0, nodes=len(nodes), tuples=len(roots))
    for root in roots:
        q, hh, t, r = root['parameters']
        require(hh == h, 'Wrong root redundancy')
        n0 = root['lower']
        base = max(r, h + t, h + 5 * t - 1, (5 * h) // 2 + 1)
        require(n0 >= base, 'Invalid shortening length')
        if n0 > base:
            require(ball(q ** t, n0 - 1, r) < q ** (h * t), 'Unproved length lower bound')
        check_length_cap(q, h, t, r, root['cap'])
        kind = root['kind']
        if kind == 'length':
            require(n0 > root['cap'], 'Noncontradictory length interval')
        elif kind == 'weight':
            index = root['node']
            require(0 <= index < len(nodes), 'Missing root proof')
            node = nodes[index]
            require(tuple(node['parameters']) == (q, n0, h, t), 'Wrong root weight state')
            require(node['value'] <= root['bound'] <= 2 * r + 2, 'Insufficient root weight bound')
        else:
            raise ValueError(f'Unresolved or unknown root kind {kind}')
        count[kind] += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--start', type=int, default=1)
    parser.add_argument('--end', type=int, default=50)
    parser.add_argument('--directory', type=Path, default=Path('certificates'))
    args = parser.parse_args()
    require(1 <= args.start <= args.end, 'Require 1 <= start <= end')
    total = dict(length=0, weight=0, nodes=0, tuples=0)
    layers = []
    for h in range(args.start, args.end + 1):
        count = verify_layer(args.directory / f'rho_{h:03d}.json', h)
        layers.append(dict(redundancy=h, **count))
        for key in total:
            total[key] += count[key]
        print(f'rho={h}: VERIFIED {count}', flush=True)
    print('ALL LAYERS VERIFIED', total, flush=True)
    (args.directory.parent / f'verified_{args.start}_{args.end}.json').write_text(
        json.dumps(dict(start=args.start, end=args.end, totals=total, layers=layers), indent=2) + '\n')

if __name__ == '__main__':
    main()
