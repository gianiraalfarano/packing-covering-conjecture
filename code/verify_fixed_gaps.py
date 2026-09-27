#!/usr/bin/env python3
"""Independent fixed-gap verifier. Imports the verifier, never the generator.

Re-enumerates the full finite range, checks all direct inequalities, and
validates exceptional length proofs and residual/iteration proof graphs.
"""
from __future__ import annotations
import argparse,json
from math import factorial
from pathlib import Path
from verify_range import require,ball,is_prime_power,check_length_cap,verify_nodes


def check(a: int, directory: Path) -> dict:
    data=json.loads((directory/f'gap_{a:02d}.json').read_text())
    require(data['format']=='packing-covering-fixed-gap-v1','Wrong format')
    require(data['gap']==a and data['min_redundancy']==51 and data['max_order']==131,
            'Wrong finite range')
    nodes=data['nodes'];verify_nodes(nodes)
    roots={tuple(z['parameters']):z for z in data['roots']}
    require(len(roots)==len(data['roots']),'Duplicate exceptional tuples')
    used=set();counts=dict(quick_length=0,quick_iteration=0,refined_length=0,weight=0)
    levels=[]
    for t in range(3,132):
        before=sum(counts.values());r=t+a;m=2*r-t+2
        e=(2*r+2)//(t+1)
        require(e*t>r and e*(t+1)<=2*r+2,'Invalid tail parameters')
        tail=(r*((e+1)*t+e-1)+e*t-r-1)//(e*t-r)
        L=a+1;B=2*(t+L)
        A=max(1,(t-2*(2**L-L-1)+2**L-2)//(2**L-1))
        S=0
        for stage in range(L):S=min(t,2*S+A+2*stage)
        require(S>=t,'Invalid iterated affine growth')
        for q in range(2,r):
            if not is_prime_power(q):continue
            powers=[1]
            for _ in range(tail+1): powers.append(powers[-1]*q)
            sums=[0]
            for z in powers: sums.append(sums[-1]+z)
            b=1
            while powers[b]<t:b+=1
            x=m-b
            require(0<=x<m,'Invalid seed rank')
            for h in range(max(51,m+1,3*r//2+1),tail):
                base=max(r,h+t,h+5*t-1,5*h//2+1)
                lower=max(base,(r*q**(t*(h-r)//r))//3+1)
                cap1=x+(2*r+1-x)*(q**(h-x)-1)//(q**(m-x)-1)
                dim=h-m+2;alpha,beta=divmod(t+1,q+1)
                u=alpha*sums[dim]
                if beta:u+=1+(beta-1)*sums[dim-1]
                cap=min(cap1,m-2+u)
                if lower>cap:
                    counts['quick_length']+=1;continue
                if lower>B and ball(q,lower-B,2)>q**(h+A-1):
                    counts['quick_iteration']+=1;continue
                key=(q,h,t,r)
                require(key in roots,f'Missing exceptional tuple {key}')
                root=roots[key];used.add(key);N=root['lower']
                require(N>=base,'Incorrect refined test length')
                if N>base:
                    require((N-1)**r < factorial(r)*q**(t*(h-r)),
                            f'Unproved factorial length lower bound {key}')
                check_length_cap(q,h,t,r,root['cap'])
                if root['kind']=='length':
                    require(N>root['cap'],'Noncontradictory refined length')
                    counts['refined_length']+=1
                elif root['kind']=='weight':
                    idx=root['node'];require(0<=idx<len(nodes),'Missing weight node')
                    nd=nodes[idx]
                    require(tuple(nd['parameters'])==(q,N,h,t),'Wrong weight state')
                    require(nd['value']<=root['bound']<=2*r+2,'Insufficient weight bound')
                    counts['weight']+=1
                else:raise ValueError('Unknown exceptional rule')
        levels.append(dict(order=t,tuples=sum(counts.values())-before))
    require(used==set(roots),'Unused exceptional records')
    require(counts==data['counts'] and levels==data['per_order'],'Wrong count summary')
    ans=dict(gap=a,**counts,tuples=sum(counts.values()),exceptions=len(roots),nodes=len(nodes))
    print('gap',a,'VERIFIED',ans,flush=True)
    return ans


def main()->None:
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--start',type=int,default=1);ap.add_argument('--end',type=int,default=5)
    ap.add_argument('--directory',type=Path,default=Path('fixed_gap_certificates'))
    args=ap.parse_args();require(1<=args.start<=args.end<=5,'Invalid gap range')
    rows=[check(a,args.directory) for a in range(args.start,args.end+1)]
    out=dict(start=args.start,end=args.end,layers=rows,
             tuples=sum(z['tuples'] for z in rows),nodes=sum(z['nodes'] for z in rows))
    (args.directory.parent/f'verified_gaps_{args.start}_{args.end}.json').write_text(json.dumps(out,indent=2)+'\n')
    print('ALL REQUESTED GAPS VERIFIED',out['tuples'],'tuples',out['nodes'],'nodes')

if __name__=='__main__':main()
