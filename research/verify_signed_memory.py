"""Signed-memory auctions: analytic equilibrium checks and Monte Carlo revenue.

Run: python verify_signed_memory.py --samples 200000 --output signed_results.json
Requires numpy. This checks explicit proven strategies; it is NOT a general
equilibrium solver. See signed_memory_derivations.md for assumptions/proofs.
"""
import argparse
from fractions import Fraction as F
import itertools
import json
import numpy as np


def exact_deviation_checks():
    """Integrate (s-z) f(z) exactly, across both triangular score laws."""
    count = 0
    for lo, hi, pieces in [
        (F(0), F(2), [(F(0), F(1), F(0), F(1)),
                      (F(1), F(2), F(2), F(-1))]),
        (F(-1), F(1), [(F(-1), F(0), F(1), F(1)),
                       (F(0), F(1), F(1), F(-1))]),
    ]:
        grid = [lo + (hi-lo)*F(i, 40) for i in range(41)]
        for s in grid:
            for t in grid:
                left, right = min(s,t), max(s,t)
                integral = F(0)
                for l, r, a, b in pieces:
                    l, r = max(l,left), min(r,right)
                    if l >= r:
                        continue
                    def primitive(z):
                        return s*a*z + (s*b-a)*z*z/2 - b*z*z*z/3
                    integral += primitive(r)-primitive(l)
                gain = integral if t >= s else -integral
                assert gain <= 0, (s,t,gain)
                count += 1
    return count


def simulate_order(x, M, order, alpha_l=1., alpha_w=2., base=5.):
    N, n, K = x.shape
    assert n == 2
    delta = alpha_w-alpha_l
    mu = np.zeros_like(x)
    revenue = np.zeros(N)
    total_utility = np.zeros(N)
    base_welfare = np.zeros(N)
    signed_cut = np.zeros(N)
    histories = []
    seen = set()
    partner = {0:1,1:0,2:3,3:2}
    W = sum(M[j,k] for r,j in enumerate(order) for k in order[r+1:])
    rows = np.arange(N)
    for j in order:
        k = partner[j]
        value = x[:,:,j]+mu[:,:,j]
        if k not in seen:
            m = M[j,k]
            assert abs(delta*m) >= 1
            u, v = x[:,:,j]-base, x[:,:,k]-base
            if m > 0:
                bid = base + delta*m + (u+v)/2
            else:
                bid = base + delta*m + (1+u-v)/2
        else:
            bid = value.copy()
        assert np.min(bid) >= 0 and np.min(value) >= 0
        winner = np.argmax(bid,axis=1)
        payment = np.min(bid,axis=1)
        revenue += payment
        total_utility += value[rows,winner]-payment
        base_welfare += x[rows,winner,j]
        for oldj, oldwinner in histories:
            signed_cut += M[oldj,j]*(oldwinner != winner)
        histories.append((j,winner))
        mu += alpha_l*M[j,:][None,None,:]
        mu[rows,winner,:] += delta*M[j,:][None,:]
        seen.add(j)
    rhs = alpha_w*W - delta*signed_cut + base_welfare-total_utility
    error = float(np.max(np.abs(rhs-revenue)))
    assert error < 1e-10
    return revenue, histories, error


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--samples',type=int,default=200000)
    parser.add_argument('--output',default='signed_results.json')
    args = parser.parse_args()
    N=args.samples
    rng=np.random.default_rng(20260920)
    x=5.+rng.random((N,2,4))
    M=np.zeros((4,4))
    M[0,1],M[1,0]=1.,-1.
    M[2,3],M[3,2]=-1.5,-1.
    records=[]
    paths={}
    max_error=0.
    for order in itertools.permutations(range(4)):
        rev,hist,error=simulate_order(x,M,order)
        max_error=max(max_error,error)
        position={j:r for r,j in enumerate(order)}
        ms=[M[j,k] if position[j]<position[k] else M[k,j]
            for j,k in [(0,1),(2,3)]]
        exact=20.+2.*23./30.+sum(2*m if m>0 else 3*m for m in ms)
        se=float(rev.std(ddof=1)/np.sqrt(N))
        mean=float(rev.mean())
        assert abs(mean-exact) <= 6*se+1e-12
        orientation=tuple(ms)
        if orientation in paths:
            assert np.max(np.abs(rev-paths[orientation]))<1e-10
        else:
            paths[orientation]=rev
        wins=dict(hist)
        for j,k in [(0,1),(2,3)]:
            m=M[j,k] if position[j]<position[k] else M[k,j]
            assert np.all((wins[j]==wins[k]) == (m>0))
        records.append({'order':''.join('ABCD'[j] for j in order),
                        'theoretical_revenue':exact,'sample_mean':mean,
                        'standard_error':se})
    records.sort(key=lambda z:-z['theoretical_revenue'])
    result={'samples':N,'seed':20260920,
            'exact_threshold_deviation_checks':exact_deviation_checks(),
            'pathwise_revenue_identity_max_error':max_error,
            'all_24_orders_checked':True,
            'interleaving_invariance_checked':True,
            'winner_pattern_checked':True,'orders':records}
    with open(args.output,'w',encoding='utf-8') as f:
        json.dump(result,f,ensure_ascii=False,indent=2)
    print(json.dumps({k:v for k,v in result.items() if k!='orders'},indent=2))
    print('Best example:',records[0])
    print('Worst example:',records[-1])


if __name__=='__main__':
    main()
