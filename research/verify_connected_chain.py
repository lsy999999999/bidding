"""Connected A->B (+), B->C (-) memory auctions.

Continuous n=2: simulate proven equilibria for all six orders.
Finite n=3: compute a symmetric logit fixed point of the auxiliary auction,
then explicitly report maximum type-conditional best-response regret.
The latter is a finite-grid approximate BNE, NOT a continuous-game solution.
Run with Python + numpy + scipy.
"""
import argparse
from fractions import Fraction as F
from itertools import permutations
import json
import numpy as np
from scipy.special import softmax
from math import comb, factorial


def q_exact(k):
    """E min of two independent sums of k U[0,1] variables."""
    answer=F(0)
    for left in range(k):
        c=[F(0)]*(k+1)
        for j in range(left+1):
            for r in range(k+1):
                c[r]+=F((-1)**j*comb(k,j)*comb(k,r)*(-j)**(k-r),factorial(k))
        surv=[(F(1) if i==0 else F(0))-c[i] for i in range(k+1)]
        sq=[F(0)]*(2*k+1)
        for i,x in enumerate(surv):
            for j,y in enumerate(surv):sq[i+j]+=x*y
        answer+=sum(co*F((left+1)**(i+1)-left**(i+1),i+1)
                    for i,co in enumerate(sq))
    return answer


def signed_tree_check(N=10000):
    K=5; base=60.; al=1.; aw=2.; delta=1.
    # Rooted out-tree; D entering an internal node exceeds 1+sum_desc(D+1).
    edges=[(0,1,5.),(0,2,-5.),(1,3,2.),(2,4,-2.)]
    M=np.zeros((K,K))
    for u,v,m in edges:M[u,v]=m
    rng=np.random.default_rng(20260922); x=base+rng.random((N,2,K)); rows=np.arange(N)
    Q={k:float(q_exact(k)) for k in range(1,K+1)}
    max_identity_error=0.; max_z=0.; records=[]
    for order in permutations(range(K)):
        pos={v:t for t,v in enumerate(order)}
        active=[(u,v,m) for u,v,m in edges if pos[u]<pos[v]]
        parent={v:(u,m) for u,v,m in active}
        comp={}; parity={}
        for v in order:
            if v not in parent:comp[v]=v; parity[v]=1
            else:
                u,m=parent[v]; comp[v]=comp[u]; parity[v]=parity[u]*(1 if m>0 else -1)
        groups={r:[v for v in range(K) if comp[v]==r] for r in set(comp.values())}
        shifts={r:sum(parity[v]*delta*abs(m) for u,v,m in active if comp[v]==r)
                for r in groups}
        theory=K*base+sum(Q[len(g)] for g in groups.values())
        for u,v,m in active:
            theory+=((al+parity[v]*delta)*m if m>0 else (aw-parity[v]*delta)*m)
        mu=np.zeros_like(x); rev=np.zeros(N); util=np.zeros(N); bw=np.zeros(N)
        cut=np.zeros(N); hist=[]; root_winners={}
        for v in order:
            value=x[:,:,v]+mu[:,:,v]
            if v in groups:
                group=groups[v]; neg=sum(parity[j]<0 for j in group)
                score=sum(parity[j]*(x[:,:,j]-base) for j in group)
                bids=base+shifts[v]+(score+neg)/len(group)
            else:bids=value.copy()
            winner=np.argmax(bids,axis=1); price=np.min(bids,axis=1)
            if v in groups:root_winners[v]=winner
            expected=root_winners[comp[v]] if parity[v]>0 else 1-root_winners[comp[v]]
            assert np.all(winner==expected)
            rev+=price; util+=value[rows,winner]-price; bw+=x[rows,winner,v]
            for u,ow in hist:cut+=M[u,v]*(ow!=winner)
            hist.append((v,winner))
            mu+=al*M[v,:][None,None,:]; mu[rows,winner,:]+=delta*M[v,:][None,:]
        W=sum(m for u,v,m in active)
        max_identity_error=max(max_identity_error,float(np.max(np.abs(rev-(aw*W-delta*cut+bw-util)))))
        se=rev.std(ddof=1)/np.sqrt(N); z=abs(rev.mean()-theory)/se
        max_z=max(max_z,float(z)); assert z<6
        records.append((theory,all((pos[u]<pos[v])==(m>0) for u,v,m in edges)))
    best=max(a for a,b in records)
    assert all(abs(a-best)<1e-10 for a,b in records if b)
    assert all(a<best-1e-8 for a,b in records if not b)
    return {'goods':K,'samples':N,'permutations':factorial(K),
            'Q_exact':{k:str(q_exact(k)) for k in range(1,K+1)},
            'all_winner_patterns_verified':True,
            'positive_forward_negative_backward_exactly_optimal':True,
            'optimal_theoretical_revenue':best,
            'max_monte_carlo_standard_errors':max_z,
            'identity_max_error':max_identity_error}


def continuous_two(N, seed):
    rng=np.random.default_rng(seed)
    base=5.; al=1.; aw=2.; p=5.; m=1.5
    delta=aw-al; P=delta*p; D=delta*m
    assert D>=1 and P>D+2
    x=base+rng.random((N,2,3)); rows=np.arange(N)
    M=np.zeros((3,3)); M[0,1]=p; M[1,2]=-m
    Q=float(F(1021,840)); results=[]
    for order in permutations(range(3)):
        mu=np.zeros_like(x); revenue=np.zeros(N); utilities=np.zeros(N)
        bw=np.zeros(N); cut=np.zeros(N); histories=[]
        # Which graph edges are activated by this order?
        pos={k:r for r,k in enumerate(order)}
        on_ab=pos[0]<pos[1]; on_bc=pos[1]<pos[2]
        for t,k in enumerate(order):
            val=x[:,:,k]+mu[:,:,k]
            if on_ab and on_bc:  # ABC only
                if t==0:
                    s=(x[:,:,0]-base)+(x[:,:,1]-base)-(x[:,:,2]-base)
                    bids=base+P-D+(s+1)/3
                else:bids=val.copy()
            elif on_ab and k==0:
                bids=base+P+((x[:,:,0]-base)+(x[:,:,1]-base))/2
            elif on_bc and k==1:
                bids=base-D+(1+(x[:,:,1]-base)-(x[:,:,2]-base))/2
            else:bids=val.copy()
            assert min(val.min(),bids.min())>=0
            w=np.argmax(bids,axis=1); price=np.min(bids,axis=1)
            revenue+=price; utilities+=val[rows,w]-price; bw+=x[rows,w,k]
            for j,oldw in histories:cut+=M[j,k]*(oldw!=w)
            histories.append((k,w))
            mu+=al*M[k,:][None,None,:]
            mu[rows,w,:]+=delta*M[k,:][None,:]
        W=sum(M[j,k] for r,j in enumerate(order) for k in order[r+1:])
        err=np.max(np.abs(revenue-(aw*W-delta*cut+bw-utilities)))
        assert err<1e-10
        if on_ab and on_bc:
            truth=3*base+Q+aw*p-(aw+delta)*m
            wins=dict(histories)
            assert np.all(wins[0]==wins[1]) and np.all(wins[1]!=wins[2])
            expected_A=3*base+Q-2*D
            assert abs(np.mean(bw-utilities)-expected_A)<.015
        elif on_ab:truth=3*base+1.1+aw*p
        elif on_bc:truth=3*base+1.1-(aw+delta)*m
        else:truth=3*base+1.
        se=float(revenue.std(ddof=1)/np.sqrt(N)); mean=float(revenue.mean())
        assert abs(mean-truth)<6*se
        results.append({'order':''.join('ABC'[k] for k in order),
                        'theory':truth,'simulation':mean,'standard_error':se,
                        'identity_max_error':float(err)})
    return results


def finite_three(step=.1):
    # Private X,Y,Z independent Bernoulli(1/2). Own type compresses to (T=X+Y,Z).
    # Opponent states are (Y,Z,q); integrate their unobserved X into their mass.
    bids=np.linspace(-2,2,int(round(4/step))+1); J=len(bids)
    types=[(t,z) for t in range(3) for z in range(2)]
    index={v:i for i,v in enumerate(types)}
    sy=np.repeat([0,0,1,1],J)
    sz=np.repeat([0,1,0,1],J)
    sq=np.tile(bids,4)
    q1=sq[:,None]; q2=sq[None,:]
    ymax=np.maximum(sy[:,None],sy[None,:])
    high=np.maximum(q1,q2)
    own=bids[:,None,None]
    eligible=own>=high
    tiecount=1+(q1==own).astype(int)+(q2==own).astype(int)
    winprob=np.where(eligible,1/tiecount,0.)
    win_cost=winprob*(-ymax-high)
    # Opponent j can win iff tied for the top including the focal bidder.
    mx=np.maximum(own,high)
    denom=(own==mx).astype(int)+(q1==mx).astype(int)+(q2==mx).astype(int)
    opp1=(q1==mx)/denom; opp2=(q2==mx)/denom
    # When own Z=1, continuation surplus equals 1 if the OTHER loser has Z=0.
    loser_surplus=opp1*(1-sz[None,None,:])+opp2*(1-sz[None,:,None])
    cost_flat=win_cost.reshape(J,-1)
    win_flat=winprob.reshape(J,-1)
    lose_flat=loser_surplus.reshape(J,-1)
    def mass_of(strategy):
        mass=np.zeros((2,2,J))
        for xx in range(2):
            for yy in range(2):
                for zz in range(2):
                    mass[yy,zz,:]+=.125*strategy[index[(xx+yy,zz)],:]
        return mass.ravel()
    def payoffs(strategy):
        mass=mass_of(strategy); joint=np.outer(mass,mass).ravel()
        pw=win_flat@joint; cw=cost_flat@joint; cl=lose_flat@joint
        return np.array([t*pw+cw+z*cl for t,z in types])
    strategy=np.ones((6,J))/J
    logs=[]
    # Fixed-point convergence is checked numerically; no general convergence claim.
    for temp in [.5,.25,.15,.1,.07,.05,.035,.025]:
        residual=float('inf')
        for it in range(12000):
            u=payoffs(strategy); target=softmax(u/temp,axis=1)
            residual=float(np.max(np.abs(target-strategy)))
            if residual<2e-10:break
            strategy=.94*strategy+.06*target
        u=payoffs(strategy)
        regret=u.max(axis=1)-(u*strategy).sum(axis=1)
        logs.append({'temperature':temp,'iterations':it+1,'fixed_point_residual':residual,
                     'max_type_regret':float(regret.max())})
        print('Finite-game check:',logs[-1],flush=True)
        if residual>1e-6:break
    # Exact finite sum for expected seller constant; no sampling of opponent actions.
    mass=mass_of(strategy)
    # For three iid states, enumerate the first winner, averaging uniform ties.
    constant=0.; cut_checks=0
    S=len(mass)
    for i in range(S):
        if mass[i]<1e-14:continue
        qa=sq[i]; qb=sq[:,None]; qc=sq[None,:]
        top=np.maximum(qa,np.maximum(qb,qc))
        nt=(qa==top).astype(int)+(qb==top).astype(int)+(qc==top).astype(int)
        second=qa+qb+qc-top-np.minimum(qa,np.minimum(qb,qc))
        # B price base part: max Y among A losers. C price: min Z among those losers.
        if_awins=np.maximum(sy[:,None],sy[None,:])+np.minimum(sz[:,None],sz[None,:])
        if_bwins=np.maximum(sy[i],sy[None,:])+np.minimum(sz[i],sz[None,:])
        if_cwins=np.maximum(sy[i],sy[:,None])+np.minimum(sz[i],sz[:,None])
        rest=((qa==top)*if_awins+(qb==top)*if_bwins+(qc==top)*if_cwins)/nt
        constant+=mass[i]*float(np.sum(np.outer(mass,mass)*(second+rest)))
    return {'type_distribution':'independent Bernoulli(1/2) for X,Y,Z',
            'bid_grid':bids.tolist(),'tie_rule':'uniform among highest bids',
            'types':types,'strategy':strategy.tolist(),'iterations':logs,
            'max_type_regret':logs[-1]['max_type_regret'],
            'auxiliary_revenue_constant':constant,
            'ABC_revenue_at_a_b_c_5_p5_m1_5_al1_aw2':15+10-1.5+constant,
            'scope':'finite-grid approximate BNE; continuous general-n equilibrium not solved'}


def main():
    p=argparse.ArgumentParser(); p.add_argument('--samples',type=int,default=200000)
    p.add_argument('--output',default='connected_chain_results.json')
    p.add_argument('--skip-finite',action='store_true')
    a=p.parse_args()
    result={'samples':a.samples,'seed':20260921,
            'continuous_two_bidder':continuous_two(a.samples,20260921)}
    print(json.dumps(result,indent=2),flush=True)
    if not a.skip_finite:result['finite_three_bidder']=finite_three()
    result['signed_tree']=signed_tree_check()
    with open(a.output,'w') as f:json.dump(result,f,indent=2)

if __name__=='__main__':main()
