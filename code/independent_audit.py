"""Separately written checker, based directly on the conditions stated in the
paper. It imports none of the other programs. Standard-library exact arithmetic only.
"""
from pathlib import Path
from fractions import Fraction as F
from math import comb,factorial
from collections import Counter
import json,sys
ROOT=Path(__file__).resolve().parents[1]
BC=NC=ROOT/'certificates'

def must(p,message):
 if not p:raise RuntimeError(message)

def volume(q,n,e):
 must(n>=0 and e>=0,'ball domain')
 return sum(comb(n,i)*(q-1)**i for i in range(min(e,n)+1))

def fields(n):
 sieve=[True]*(n+1)
 if n>=0:sieve[0]=False
 if n>=1:sieve[1]=False
 out=set()
 for p in range(2,n+1):
  if sieve[p]:
   for j in range(p*p,n+1,p):sieve[j]=False
   v=p
   while v<=n:out.add(v);v*=p
 return sorted(out)

def chain(q,n,h,t,r,route):
 # Reconstruct the residual parameters after every actual support removal.
 remaining_n,remaining_h,remaining_k=n,h,n-h
 support=order=0
 must(n>2*r+2 and h+t>2*r+2,'nontrivial domain')
 for e,s in route:
  must(type(e)==int and type(s)==int and e>0 and 0<s<=t-order,'step domain')
  cost=e*(s+1)
  must(s<=remaining_k and cost-s<=remaining_h,'dimension or padding')
  must(volume(q,remaining_n,e)>q**(remaining_h+s-1),'strict fiber condition')
  remaining_n-=cost;remaining_k-=s;remaining_h-=cost-s
  order+=s;support+=cost
  must(remaining_n-remaining_k==remaining_h and support<=2*r+2,'state')
 must(order==t and support<=2*r+2,'endpoint')
 return len(route)

def cap(q,h,t,r):
 m=2*r-t+2
 # Check rational averages before taking floors, independently of original integer formulas.
 candidates=[x+int(F((2*r+1-x)*(q**(h-x)-1),q**(m-x)-1)) for x in range(m)]
 a,b=divmod(t+1,q+1);v=h-m+2
 P=lambda v:(q**v-1)//(q-1)
 candidates.append(m-2+a*P(v)+(1+(b-1)*P(v-1) if b else 0))
 return min(candidates)

def numeric(row):
 q,n,h,t,r=[row[k] for k in ['q','n','h','t','r']]
 simple=max(h+5*t-1,5*h//2+1)
 must(n>=simple,'test length domain')
 if n>simple:must(volume(q**t,n-1,r)<q**(t*h),'covering lower length')
 if row['kind']=='length':must(n>cap(q,h,t,r),'length contradiction');return 0
 must(row['kind']=='chain','unexpected certificate kind')
 return chain(q,n,h,t,r,row['chain'])

def symbolic(row,old=False):
 t=row['t'];r=t+row['a'] if old else row['r'];S=W=0
 for e,s in row['chain']:
  must(type(e)==int and type(s)==int and 0<s<=t-S and e*t>=r,'symbolic domain')
  if old:
   a=row['a'];must(12*e<=r,'old radius domain')
   allowance=F(e*(t+1)-2*r)+F(2*(e*t-r),5*a)
  else:
   must(row['q0']==2 and t>=4,'new symbolic domain')
   K=F(2*r,11)-F(W+e-1,2**(t+1))
   must(K>0 and K**e>=factorial(e),'rational factorial condition')
   allowance=e*(t+1)-2*r
  must(allowance>=S+s-1-W,'symbolic exponent')
  S+=s;W+=e*(s+1)
 must(S==t and W<=2*r+2,'symbolic endpoint')

def tail_margins(z,r):
 t,s=z['t'],z['s'];u=t-s;v=t+1
 lo=F(2*r+s-1,v);hi=lo+F(t,v)
 wl=(s+1)*lo;wh=(s+1)*hi
 fa=F(r+t-1,t);fb=F(2*r-wl+2*t-1,v)
 c1,c2=F(*z['c1']),F(*z['c2'])
 H=lambda W,e,c:F(r,2)-F(11*(W+e-1),2**(t+3))-c*e
 return [2*r+2-wh-(u+1)*fa,2*r+2-F(s,v)*wh-F((u+1)*(2*r+2*t-1),v),H(0,hi,c1),H(wh,fa,c2),H(wh,fb,c2)]

def stages(goal,e,A):
 S=L=0
 while S<goal:S=e*S+A+e*L;L+=1
 return L

def ceil_power(n,d,b):
 L=0
 while d*b**L<n:L+=1
 return L

def analytic_bases():
 states=0;mins=[]
 for t in range(32,64):
  slack=[]
  for a in range(t//3+1,t+1):
   b=max(0,(2*a-t+2)//3)
   val=a+2-b-4*stages(b,4,5)-3*stages(a-b,3,4)-2*stages(t-a,2,3)
   must(val>=0,'three-phase base');slack.append(val);states+=1
  mins.append([t,min(slack)])
 for t in range(64,128):
  f=3*ceil_power(t+4,4,3)+2*ceil_power(t+3,3,2)
  j=4*ceil_power(t+7,5,4)+3*ceil_power(t+2,2,3)+2*ceil_power(t+6,6,2)
  must(3*f<=t+6 and 6*j<=3*t+8,'dyadic base')
 for delta,r in [(0,3),(1,3),(2,4),(3,9)]:
  must(volume(8,2*r+2**(delta+3)-3,r)<8**(2*r+delta),'binary order3 endpoint')
 for q,delta in [(3,0),(3,1),(3,2),(4,0),(4,1),(4,2),(5,2),(7,2),(8,2),(9,2)]:
  r=q+1;h=delta+3;a,b=divmod(4,q+1)
  U=a*(q**h-1)//(q-1)+(1+(b-1)*(q**(h-1)-1)//(q-1) if b else 0)
  must(32*volume(q**3,2*r+U-3,r)<q**(3*(2*r+delta)),'nonbinary order3 endpoint')
 return {'three_phase_states':states,'dyadic_orders':64,'stage_minima':mins}

def order_three_tails():
 for q,R,c1,j1,c2,j2 in [(3,48,F(6,5),24,F(3,2),8),(4,26,F(4,3),13,F(7,4),5)]:
  for c,j in [(c1,j1),(c2,j2)]:must(c**j>=3*j and c*j>=j+1,'order-three prefactor')
  must(F(R,2)>=j1 and F(R,3)>=j2,'order-three minimum radii')
  def margins(r):
   e=F(r+1,2);f=F(r+2,3)
   H=lambda W,E,c:F(q-1,q)*(r-F(11*(W+E-1),4*q**3))-c*E
   return H(0,e,c1),H(r,f,c2)
  must(all(a>=0 and b>=a for a,b in zip(margins(R),margins(R+1))),'order-three affine tails')
 must(F(19,7)*(F(10*192,193)-1)>=24,'binary low endpoint')
 must(F(19,7)*(F(35*192,3*193)-1)>=28,'binary upper endpoint')
 must(F(3,2)**3>F(343,256),'binary concavity endpoint')
 must(F(1343,772)**2>3,'binary high prefactor')
 must(F(343,256)**32>2304 and F(27,16)**16>2304,'binary initial powers')
 must(F(87,6*343)>F(1,192) and F(11,12*27)>F(1,192),'binary tail propagation')
 return 'all three field regimes verified with rational inequalities'

def cutoff_boundaries():
 count=0
 for t in range(3,132):
  for r in range(t,max(193,t+133)):
   E=(2*r+2)//(t+1);den=E*t-r
   must(den>0 and E*(t+1)<=2*r+2 and E*t<=2*r,'cutoff radius')
   num=r*((E+1)*t+E-1);D=(num+den-1)//den
   must((D-1)*den<num<=D*den,'closed cutoff endpoint')
   count+=1
 return count

def populations():
 def D(t,r):
  e=(2*r+2)//(t+1)
  return -(-r*((e+1)*t+e-1)//(e*t-r))
 pp={r:len(fields(r-1)) for r in range(1,138)}
 rho=0
 for h in range(1,51):
  for t in range(3,h-4):
   for r in range(t+1,(h+t-3)//2+1):
    if 3*r<2*h and (t,r) not in {(3,4),(3,5),(4,5),(4,6)}:rho+=pp[r]
 gaps=0
 for a in range(1,6):
  for t in range(3,132):
   r=t+a
   gaps+=pp[r]*max(0,D(t,r)-max(51,2*r-t+3,3*r//2+1))
 boundary=0
 for t in range(3,32):
  for r in range(t+1,max(2*t,t+15)+1):
   boundary+=pp[r]*max(0,D(t,r)-max(51,2*r-t+3,3*r//2+1))
 must((rho,gaps,boundary)==(59086,2396360,619126),'finite populations')
 return {'redundancy':rho,'fixed_gaps':gaps,'boundary':boundary}

def main():
 out={'independent_population_counts':populations(),'analytic_bases':analytic_bases(),'order_three_tails':order_three_tails(),'cutoff_boundary_pairs':cutoff_boundaries()}
 tails=json.loads((NC/'tails.json').read_text())
 must([z['t'] for z in tails]==list(range(4,32)),'tail coverage')
 for z in tails:
  t,s,R=z['t'],z['s'],z['r0']
  must(1<=s<t and R>2*t,'tail domain')
  for cj,jj in [('c1','j1'),('c2','j2')]:
   c=F(*z[cj]);j=z[jj]
   must(type(j)==int and j>=1 and c>1 and c**j>=3*j and c*j>=j+1,'prefactor')
  must(F(2*R+s-1,t+1)>=z['j1'] and F(R,t)>=z['j2'],'radius floor')
  at0=tail_margins(z,0);at1=tail_margins(z,1)
  slopes=[b-a for a,b in zip(at0,at1)]
  must(all(a>=0 for a in slopes) and all(a+R*b>=0 for a,b in zip(at0,slopes)),'infinite affine domain')
 out['tail_rows']=len(tails)
 syold=json.loads((BC/'extensions/symbolic.json').read_text())
 must(Counter((z['t'],z['a']) for z in syold)==Counter({(t,a):1 for a in range(6,10) for t in range(32,5*a*a)}),'old symbolic coverage')
 for z in syold:symbolic(z,True)
 data=json.loads((NC/'new_finite.json').read_text())
 sy={(z['t'],z['r']):z for z in data['symbolic']}
 must(len(sy)==len(data['symbolic']),'new symbolic duplicates')
 for z in sy.values():symbolic(z)
 roots={(z['q'],z['h'],z['t'],z['r']):z for z in data['roots']}
 must(len(roots)==len(data['roots']),'new numerical duplicates')
 expected=set();expected_sy=set();cuts={z['t']:z['r0'] for z in tails}
 for t in range(3,32):
  bound=192 if t==3 else cuts[t]
  for r in range(max(t+16,2*t+1),bound):
   if (t,r) in sy:expected_sy.add((t,r));continue
   E=(2*r+2)//(t+1);D=-(-r*((E+1)*t+E-1)//(E*t-r))
   for q in fields(r-1):
    if t==3 and r>={2:192,3:48}.get(q,26):continue
    h0=max(51,2*r+(4 if q==2 else 3)) if t==3 else max(51,2*r+1)
    expected.update((q,h,t,r) for h in range(h0,D))
 must(expected==set(roots) and expected_sy==set(sy),'finite exhaustion')
 counts=Counter();steps=0
 for z in roots.values():steps+=numeric(z);counts[z['kind']]+=1
 out['new_numeric_counts']=dict(counts);out['new_chain_steps']=steps
 ex=exsteps=0
 for f in sorted((BC/'extensions').glob('order_*.json')):
  d=json.loads(f.read_text())
  for z in d['roots']:
   if z['kind']=='chain':
    exsteps+=chain(z['q'],z['n'],z['h'],z['t'],z['r'],z['chain']);ex+=1
 out['extension_chains']=ex;out['extension_chain_steps']=exsteps
 out['old_symbolic']=len(syold);out['new_symbolic']=len(sy)
 nodes=files_count=0
 for folder in ['redundancy','fixed_gaps','extensions']:
  for f in (BC/folder).glob('*.json'):
   d=json.loads(f.read_text())
   if not isinstance(d,dict) or 'nodes' not in d:continue
   ns=d['nodes'];seen=set();todo=[z['node'] for z in d['roots'] if z['kind']=='weight']
   while todo:
    i=todo.pop()
    if i in seen:continue
    must(0<=i<len(ns),'node id');seen.add(i);n=ns[i]
    if n['kind']=='residual':todo+=[n['support_bound']]+[j for _,j in n['children']]
   must(seen==set(range(len(ns))),'unused DAG nodes');nodes+=len(ns);files_count+=1
 out['reachable_nodes']=nodes;out['DAG_files']=files_count
 out['status']='PASS'
 print(json.dumps(out,indent=2))

if __name__=='__main__':main()
