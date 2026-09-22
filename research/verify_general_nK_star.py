"""Connected signed star: arbitrary-n,K theorem illustration at n=3,K=5.

Finite types: all private noises independent Bernoulli(1/2).
Finite root bids: grid 0.1, after normalization q in [-2,3].
Two positive leaves, two negative leaves; negative leaves sold before root.
Computes a symmetric approximate BNE of the root auxiliary auction and
checks its type-wise best-response regret. Not a continuous-game solver.
Dependencies: numpy, scipy. Run: python verify_general_nK_star.py
"""
import json
import argparse
from math import comb
from fractions import Fraction as F
import numpy as np
from scipy.special import softmax


def order_stats_sum(n,d):
    """Exact top and second-top expectations for Binomial(d,1/2) signals."""
    top=F(0); second=F(0); cdf=F(0)
    for s in range(d):
        cdf+=F(comb(d,s),2**d)
        top+=1-cdf**n
        second+=1-cdf**n-n*(1-cdf)*cdf**(n-1)
    return top,second


def solve_root():
    n=3; r=2; qnegative=2; K=5
    bids=np.linspace(-2,3,51); J=len(bids)
    leaf_types=np.array([[0,0],[0,1],[1,0],[1,1]])
    leaf=np.repeat(leaf_types,J,axis=0); sq=np.tile(bids,4)
    q1=sq[:,None]; q2=sq[None,:]; own=bids[:,None,None]
    high=np.maximum(q1,q2)
    tiecount=1+(q1==own).astype(int)+(q2==own).astype(int)
    winprob=np.where(own>=high,1/tiecount,0.)
    leaf_cost=np.maximum(leaf[:,None,:],leaf[None,:,:]).sum(axis=2)
    win_cost=winprob*(-leaf_cost-high)
    pflat=winprob.reshape(J,-1); cflat=win_cost.reshape(J,-1)
    scores=np.arange(4)
    allowed=(bids[None,:]>=scores[:,None]-r-1e-9)&(bids[None,:]<=scores[:,None]+1e-9)
    def mass_of(strategy):
        mass=np.zeros((4,J))
        for state,y in enumerate(leaf_types):
            for x in [0,1]:mass[state]+=.125*strategy[int(x+y.sum())]
        return mass.ravel()
    def payoffs(strategy):
        mass=mass_of(strategy); joint=np.outer(mass,mass).ravel()
        return scores[:,None]*(pflat@joint)[None,:]+(cflat@joint)[None,:]
    strategy=allowed/allowed.sum(axis=1,keepdims=True)
    logs=[]
    for temp in [.5,.25,.15,.1,.07,.05,.035,.025]:
        for it in range(15000):
            u=payoffs(strategy)
            target=softmax(np.where(allowed,u/temp,-np.inf),axis=1)
            residual=float(np.max(np.abs(target-strategy)))
            if residual<2e-10:break
            strategy=.94*strategy+.06*target
        u=payoffs(strategy)
        # Include all grid actions, also the weakly dominated actions excluded above.
        regret=u.max(axis=1)-(u*strategy).sum(axis=1)
        item={'temperature':temp,'iterations':it+1,'residual':residual,
              'maximum_type_regret':float(regret.max())}
        logs.append(item); print(item,flush=True)
        if residual>1e-6:break
    mass=mass_of(strategy); C=0.
    for i in range(len(mass)):
        if mass[i]<1e-15:continue
        qa=sq[i]; qb=sq[:,None]; qc=sq[None,:]
        top=np.maximum(qa,np.maximum(qb,qc))
        nt=(qa==top).astype(int)+(qb==top).astype(int)+(qc==top).astype(int)
        second=qa+qb+qc-top-np.minimum(qa,np.minimum(qb,qc))
        cost_a=leaf_cost
        cost_b=np.maximum(leaf[i,:],leaf[None,:,:]).sum(axis=2)
        cost_c=np.maximum(leaf[i,:],leaf[:,None,:]).sum(axis=2)
        future=((qa==top)*cost_a+(qb==top)*cost_b+(qc==top)*cost_c)/nt
        C+=mass[i]*float(np.sum(np.outer(mass,mass)*(second+future)))
    t1,s1=order_stats_sum(n,1); ts,ss=order_stats_sum(n,r+1)
    error=ts-ss+r*(1-s1)+qnegative*(t1-s1)
    lower=50+10+ss-r+(K-1)*s1
    upper=50+10+ts+qnegative*t1
    revenue=50+10+C+qnegative*float(s1)
    assert revenue>=float(lower)-1e-8
    assert float(upper-lower)==float(error)
    # Every type is checked; the full multi-round game has exact best-response
    # continuations, so only this root-grid regret remains.
    return {'n':n,'K':K,'positive_leaf_weights':[2,3],
            'negative_leaf_weights':[-1.5,-2],'alpha_L':1,'alpha_W':2,
            'base_values':[10]*K,'type_distribution':'iid Bernoulli(1/2) noise',
            'tie_rule':'uniform among highest bids','root_bid_grid':bids.tolist(),
            'root_strategies_by_private_sum':strategy.tolist(),'iterations':logs,
            'maximum_grid_type_regret':logs[-1]['maximum_type_regret'],
            'auxiliary_revenue_constant':C,'computed_policy_revenue':revenue,
            'guaranteed_policy_revenue_lower':str(lower),
            'all_order_equilibrium_revenue_upper':str(upper),
            'additive_ordering_guarantee':str(error),
            'scope':'finite-game approximate BNE; theorem existence is for finite games; no continuous-equilibrium claim'}


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',default='general_nK_star_results.json')
    args=parser.parse_args()
    result=solve_root()
    with open(args.output,'w') as f:json.dump(result,f,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k not in
                     ['root_strategies_by_private_sum','root_bid_grid','iterations']},indent=2))
