import itertools
import json
from fractions import Fraction as Q
import numpy as np
from scipy.integrate import quad

def fs(s):
    if s <= 0: return 0.
    if s <= 1: return s*s/2
    if s < 2: return 1-(2-s)**2/2
    return 1.

def h(s,n):
    if s <= 0: return 0.
    if s >= 2: return 1.
    lo,hi=max(0.,s-1),min(1.,s)
    def integrand(y):
        top=np.clip((y-lo)/(hi-lo),0,1)
        if s <= 1:
            trunc=1. if y>=s else (s*y-y*y/2)/fs(s)
        else:
            trunc=(y if y<=s-1 else s*y-y*y/2-(s-1)**2/2)/fs(s)
        return 1-top*trunc**(n-2)
    return quad(integrand,0,1,points=[lo,hi],epsabs=1e-10)[0]

rng=np.random.default_rng(20260919)
out={}
grid=np.linspace(0,2,1001)
for n in [2,3,5]:
    hv=np.array([h(s,n) for s in grid])
    bids=grid+1.2-hv
    assert np.min(np.diff(bids))>0
    x=rng.random((500000,n));y=rng.random((500000,n));s=x+y
    winner=np.argmax(s,axis=1)
    second=np.partition(s,-2,axis=1)[:,-2]
    pa=np.interp(second,grid,bids)
    y[np.arange(len(y)),winner]=-np.inf
    rev=pa+.3+np.max(y,axis=1)
    target=1.5+quad(lambda z:1-fs(z)**n-n*(1-fs(z))*fs(z)**(n-1),0,2,points=[1])[0]
    se=rev.std()/len(rev)**.5
    assert abs(rev.mean()-target)<6*se+1e-5
    out[f'strong_memory_n{n}']={'expected':target,'MC':float(rev.mean()),'SE':float(se),'min_bid_grid_increment':float(np.diff(bids).min())}

out['exact_n2_second_sum']=str(Q(23,30))
for d in [0.,.2,.7,1.,1.5]:
    numeric=quad(lambda z:(1-z)*(1-max(0,z-d)),0,1,points=[min(1,d)])[0]
    exact=1/3+d/2-d*d/2+d**3/6 if d<=1 else .5
    assert abs(numeric-exact)<1e-10
    out[f'delayed_shift_{d}']=exact

K=5
M=rng.integers(0,5,(K,K)).astype(float)
np.fill_diagonal(M,0)
delta=1.;alpha_l=1.;H=1+M.sum()
def solve(pi,mu,t=0):
    if t==K:return np.zeros(2),0.
    j=pi[t]
    nxt=[]
    for winner in [0,1]:
        updated=mu+alpha_l*M[j]
        updated=updated.copy();updated[winner]+=delta*M[j]
        nxt.append(solve(pi,updated,t+1))
    v=np.array([H,0.])+mu[:,j]
    b=np.array([v[0]+nxt[0][0][0]-nxt[1][0][0],v[1]+nxt[1][0][1]-nxt[0][0][1]])
    expected=np.array([v[0]+2*delta*sum(M[j,k] for k in pi[t+1:]),v[1]])
    assert np.allclose(b,expected)
    assert b[0]>b[1]>=0
    winner=int(np.argmax(b));price=b[1-winner]
    util=nxt[winner][0].copy();util[winner]+=v[winner]-price
    return util,price+nxt[winner][1]

brute=[]
for pi in itertools.permutations(range(K)):
    _,r=solve(pi,np.zeros((2,K)))
    target=sum(M[pi[s],pi[t]] for s in range(K) for t in range(s+1,K))
    assert abs(r-target)<1e-8
    brute.append(target)
dp={0:0.}
for mask in range(1,1<<K):
    dp[mask]=max(dp[mask^(1<<k)]+sum(M[j,k] for j in range(K) if j!=k and mask>>j&1) for k in range(K) if mask>>k&1)
assert dp[(1<<K)-1]==max(brute)
out['generalK_dominance']={'permutations_verified':len(brute),'optimum':max(brute),'DP':dp[(1<<K)-1]}
print(json.dumps(out,indent=2))
