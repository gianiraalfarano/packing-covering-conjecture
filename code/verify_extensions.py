"""Checker for the 1,022 symbolic sequences of Proposition 4.7 and for the
range of Theorem 4.18 (r <= 2t or r-t <= 15, 3 <= t <= 31).
The search program that produced the records is not needed and not included. All parameter tuples are re-enumerated.
"""
import argparse,json
from math import factorial,isqrt
from pathlib import Path
from verify_range import require,ball,is_prime_power,check_length_cap,verify_nodes

def check_chain(q,n,h,t,B,route):
    W=S=0
    require(n>B and h+t>B,'Nontrivial chain domain')
    for e,s in route:
        require(isinstance(e,int) and isinstance(s,int) and e>=1 and 1<=s<=t-S,'Invalid chain step')
        require(h-W+S>=0 and n-h-S>=s,'Invalid residual dimensions')
        require(ball(q,n-W,e)>q**(h-W+S+s-1),'False strict volume inequality')
        W+=e*(s+1);S+=s
        require(W<=B,'Support budget exceeded')
    require(S==t,'Incomplete nullity')
    require(W-S<h,'Insufficient ambient rank for padding')

def symbolic(directory):
    rows=json.loads((directory/'symbolic.json').read_text())
    expected={(t,a) for a in range(6,10) for t in range(32,5*a*a)}
    require(len(rows)==len(expected),'Wrong symbolic row count')
    actual=set()
    for row in rows:
        t,a=row['t'],row['a'];actual.add((t,a));r=t+a;W=S=0
        for e,s in row['chain']:
            require(e>=1 and 1<=s<=t-S and e*t>=r and r>=12*e,'Invalid symbolic parameters')
            # Clear the denominator 5a. This checks the analytic exponent
            # directly, without a logarithm.
            need=s-1-W+S
            require(5*a*need <= 5*a*(e*t-2*r+e)+2*(e*t-r),'False symbolic exponent')
            W+=e*(s+1);S+=s
            require(W<=2*r+2,'Symbolic support budget')
        require(S==t,'Incomplete symbolic chain')
    require(actual==expected,'Missing or repeated symbolic pair')
    print('SYMBOLIC VERIFIED',len(rows),'paths',flush=True)
    return len(rows)

def layer(directory,t,gap):
    data=json.loads((directory/f'order_{t:03}.json').read_text())
    require(data['format']=='variable-affine-v1' and data['t']==t and data['gap']==gap,'Wrong layer metadata')
    verify_nodes(data['nodes'])
    roots={(z['q'],z['h'],z['r']):z for z in data['roots']}
    require(len(roots)==len(data['roots']),'Repeated exception')
    used=set();counts=dict(length=0,iteration=0,chain=0,weight=0)
    for r in range(t+1,max(2*t,t+gap)+1):
        a=r-t;m=2*r-t+2;B=2*r+2
        E=B//(t+1)
        require(E*t>r,'Invalid tail denominator')
        D=(r*((E+1)*t+E-1)+E*t-r-1)//(E*t-r)
        L=a+1;A=max(1,(t-2*(2**L-L-1)+2**L-2)//(2**L-1))
        S=0
        for i in range(L):S=min(t,2*S+A+2*i)
        require(S==t,'Invalid direct iteration growth')
        for q in range(2,r):
            if not is_prime_power(q):continue
            b=1
            while q**b<t:b+=1
            x=m-b
            require(0<=x<m,'Invalid incidence seed')
            for h in range(max(51,m+1,3*r//2+1),D):
                base=max(h+5*t-1,5*h//2+1,h+t,r)
                N=max(base,r*q**(t*(h-r)//r)//3+1)
                cap=x+(2*r+1-x)*(q**(h-x)-1)//(q**(m-x)-1)
                dim=h-m+2;alpha,beta=divmod(t+1,q+1)
                u=alpha*(q**dim-1)//(q-1)
                if beta:u+=1+(beta-1)*(q**(dim-1)-1)//(q-1)
                cap=min(cap,m-2+u)
                if N>cap:counts['length']+=1;continue
                if N>B and ball(q,N-B,2)>q**(h+A-1):counts['iteration']+=1;continue
                key=(q,h,r);require(key in roots,f'Missing tuple {(q,h,t,r)}')
                row=roots[key];used.add(key);n=row['n']
                require(row['t']==t and n>=base,'Invalid shortening length')
                if n>base:
                    if row.get('exact_length'):
                        require(ball(q**t,n-1,r)<q**(t*h),'Unproved exact covering lower bound')
                    else:
                        require((n-1)**r<factorial(r)*q**(t*(h-r)),'Unproved factorial lower bound')
                check_length_cap(q,h,t,r,row['cap'])
                kind=row['kind']
                if kind=='length':require(n>row['cap'],'Invalid length contradiction')
                elif kind=='chain':check_chain(q,n,h,t,B,row['chain'])
                elif kind=='weight':
                    i=row['node'];require(0<=i<len(data['nodes']),'Missing weight root')
                    nd=data['nodes'][i]
                    require(tuple(nd['parameters'])==(q,n,h,t) and nd['value']<=row['bound']<=B,'Insufficient weight proof')
                else:raise ValueError('Unknown rule')
                counts[kind]+=1
    require(used==set(roots),'Unused exceptions')
    require(counts==data['counts'],'Incorrect summary')
    result=dict(t=t,**counts,tuples=sum(counts.values()),exceptions=len(roots),nodes=len(data['nodes']))
    print('VERIFIED',result,flush=True)
    return result

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('--directory',type=Path,default=Path('certificates/extensions'))
    ap.add_argument('--max-t',type=int,default=31);ap.add_argument('--gap',type=int,default=15)
    args=ap.parse_args();require(args.max_t==31 and args.gap==15,'Final theorem range must be checked in full')
    n=symbolic(args.directory)
    rows=[layer(args.directory,t,args.gap) for t in range(3,32)]
    totals={k:sum(z[k] for z in rows) for k in ['length','iteration','chain','weight','tuples','exceptions','nodes']}
    output=dict(symbolic_paths=n,totals=totals,layers=rows)
    (args.directory/'verified.json').write_text(json.dumps(output,indent=2)+'\n')
    print('ALL EXTENSIONS VERIFIED',totals,flush=True)
