"""Checker for the large-radius bounds of Section 4.4 and the final finite
set of Theorem 4.19. The search program that produced the records is not
needed and is not part of this repository.
All proof comparisons and all affine-tail checks use exact arithmetic.
"""
import json
from pathlib import Path
from fractions import Fraction as Q
from math import factorial,isqrt
ROOT=Path(__file__).resolve().parents[1]

def need(ok,msg):
 if not ok:raise ValueError(msg)

def volume(q,n,e):
 need(n>=0 and e>=0,'ball domain')
 z=term=1
 for i in range(1,min(n,e)+1):
  term=term*(n-i+1)*(q-1)//i;z+=term
 return z

def primepower(q):
 p=next((d for d in range(2,isqrt(q)+1) if q%d==0),q)
 while q%p==0:q//=p
 return q==1

def verify_tails():
 rows=json.loads((ROOT/'certificates/tails.json').read_text())
 need([z['t'] for z in rows]==list(range(4,32)),'tail coverage')
 for z in rows:
  t,s,R=z['t'],z['s'],z['r0'];v=t+1;u=t-s
  need(1<=s<t and R>2*t,'tail domain')
  c1,c2=Q(*z['c1']),Q(*z['c2']);j1,j2=z['j1'],z['j2']
  for c,j in [(c1,j1),(c2,j2)]:need(c**j>=3*j and c>=Q(j+1,j),'factorial prefactor induction')
  need(Q(2*R+s-1,v)>=j1 and Q(R,t)>=j2,'minimum radii')
  # Compute the five affine expressions independently by evaluation at R,R+1.
  def margins(r):
   ehigh=Q(2*r+s+t-1,v);elow=Q(2*r+s-1,v)
   whigh=(s+1)*ehigh;wlow=(s+1)*elow
   fhigh=[Q(r+t-1,t),Q(2*r-wlow+2*t-1,v)]
   support1=2*r+2-whigh-(u+1)*fhigh[0]
   support2=2*r+2-Q(s,v)*whigh-Q((u+1)*(2*r+2*t-1),v)
   def ballmargin(w,e,c):return Q(1,2)*(r-Q(11*(w+e-1),4*2**t))-c*e
   return [support1,support2,ballmargin(0,ehigh,c1),*[ballmargin(whigh,f,c2) for f in fhigh]]
  base=margins(R);nxt=margins(R+1)
  need(all(a>=0 and b>=a for a,b in zip(base,nxt)),'negative affine tail margin')
 print('RATIONAL TAILS VERIFIED',len(rows),flush=True)
 # Order-three tails, including q>=4 as one uniform family.
 for q,R,c,j,cc,jj in [(3,48,Q(6,5),24,Q(3,2),8),(4,26,Q(4,3),13,Q(7,4),5)]:
  need(c**j>=3*j and c>=Q(j+1,j) and cc**jj>=3*jj and cc>=Q(jj+1,jj),'order3 prefactors')
  need(Q(R,2)>=j and Q(R,3)>=jj,'order3 radius lower bounds')
  def margins(r):
   first=Q(q-1,q)*(r-Q(11*(Q(r+1,2)-1),4*q**3))-c*Q(r+1,2)
   second=Q(q-1,q)*(r-Q(11*(r+Q(r+2,3)-1),4*q**3))-cc*Q(r+2,3)
   return [first,second]
  need(all(a>=0 and b>=a for a,b in zip(margins(R),margins(R+1))),'order3 affine inequalities')
 need(Q(1343,772)**2>3,'binary square root comparison')
 need(Q(343,256)**32>12*192 and Q(27,16)**16>12*192,'binary endpoints')
 need(Q(87,343*6)>Q(1,192) and Q(11,27*12)>Q(1,192),'binary monotonicity')
 print('ORDER THREE TAILS VERIFIED: q=2 r>=192; q=3 r>=48; q>=4 r>=26',flush=True)
 return {z['t']:z['r0'] for z in rows}

def symbolic(z):
 t,r,q=z['t'],z['r'],z['q0'];W=S=0
 need(q==2 and t>=4,'symbolic domain')
 for e,s in z['chain']:
  need(e>=1 and 1<=s<=t-S and e*t>=r,'symbolic orders')
  # Independently reconstruct the rational lower bound on (n-W-e+1)(q-1)/q.
  k=Q(q-1,q)*(Q(4*r,11)-Q(W+e-1,q**t))
  need(k>0 and k**e>=factorial(e),'symbolic factorial inequality')
  need(e*(t+1)-2*r>=s-1-W+S,'symbolic exponent')
  W+=e*(s+1);S+=s
  need(W<=2*r+2,'symbolic support')
 need(S==t,'symbolic completion')

def finite(tails):
 data=json.loads((ROOT/'certificates/new_finite.json').read_text())
 need(data['format']=='new-frontier-v1','format')
 sy={(z['t'],z['r']):z for z in data['symbolic']}
 need(len(sy)==len(data['symbolic']),'duplicate symbolic row')
 for z in sy.values():symbolic(z)
 roots={(z['q'],z['h'],z['t'],z['r']):z for z in data['roots']}
 need(len(roots)==len(data['roots']),'duplicate root')
 used=set();usy=set();counts=dict(symbolic_pairs=0,length=0,chain=0);per=[]
 for t in range(3,32):
  before=counts.copy();upper=192 if t==3 else tails[t]
  for r in range(max(2*t+1,t+16),upper):
   if (t,r) in sy:usy.add((t,r));counts['symbolic_pairs']+=1;continue
   e=(2*r+2)//(t+1);den=e*t-r
   need(den>0,'tail denominator')
   D=(r*((e+1)*t+e-1)+den-1)//den
   for q in range(2,r):
    if not primepower(q):continue
    bound=(192 if q==2 else 48 if q==3 else 26) if t==3 else tails[t]
    if r>=bound:continue
    h0=max(51,2*r+3+int(q==2)) if t==3 else max(51,2*r+1)
    for h in range(h0,D):
     key=q,h,t,r;need(key in roots,f'missing {key}');z=roots[key];used.add(key)
     n=z['n'];base=max(h+5*t-1,5*h//2+1)
     need(n>=base and n>=h+t,'test length')
     if n>base:
      # Exact covering-volume comparison (no root estimate is trusted).
      need(volume(q**t,n-1,r)<q**(t*h),'length not certified')
     if z['kind']=='length':
      m=2*r-t+2
      caps=[x+(2*r+1-x)*(q**(h-x)-1)//(q**(m-x)-1) for x in range(m)]
      alpha,beta=divmod(t+1,q+1);dim=h-m+2
      u=alpha*sum(q**i for i in range(dim))
      if beta:u+=1+(beta-1)*sum(q**i for i in range(dim-1))
      need(n>min(caps+[m-2+u]),'length comparison')
      counts['length']+=1
     elif z['kind']=='chain':
      W=S=0;B=2*r+2
      need(n>B and h+t>B,'chain domain')
      for e,s in z['chain']:
       need(isinstance(e,int) and isinstance(s,int) and e>=1 and 1<=s<=t-S,'step domain')
       need(n-h-S>=s and h-W+S>=0,'residual dimension')
       need(volume(q,n-W,e)>q**(h-W+S+s-1),'strict affine fiber inequality')
       S+=s;W+=e*(s+1)
       need(W<=B and W-S<h,'support and padding rank')
      need(S==t,'incomplete chain');counts['chain']+=1
     else:raise ValueError('unknown rule')
  row=dict(t=t,**{k:counts[k]-before[k] for k in counts});per.append(row);print('FINITE VERIFIED',row,flush=True)
 need(used==set(roots) and usy==set(sy),'unused records')
 need(counts==data['counts'] and per==data['per_order'],'incorrect counts')
 print('ALL NEW CERTIFICATES VERIFIED',counts,flush=True)
 (ROOT/'outputs').mkdir(exist_ok=True)
 (ROOT/'outputs/new_verified.json').write_text(json.dumps(dict(counts=counts,per_order=per,status='verified'),indent=2)+'\n')

if __name__=='__main__':finite(verify_tails())
