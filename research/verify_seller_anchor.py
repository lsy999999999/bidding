from fractions import Fraction as F
from itertools import permutations
import json
from random import Random

K, n = 5, 3
alpha_l, alpha_w = F(1), F(2)
delta = alpha_w-alpha_l
theta = [F(1,10), F(2,5), F(9,10)]
a = [F(0),F(100),F(0),F(0),F(0)]
lam = [F(1)]*K
M = [[F(0) for _ in range(K)] for _ in range(K)]
M[2][3],M[3][4],M[4][2]=F(3),F(2),F(1)
ordinary_total=sum(sum(row) for row in M)
B=ordinary_total+K+1
for k in range(1,K):M[0][k]=B

def W(pi):
    return sum(M[pi[s]][pi[t]] for s in range(K) for t in range(s+1,K))

def update(mu,j,w):
    return [[mu[i][k]+(alpha_w if i==w else alpha_l)*M[j][k] for k in range(K)] for i in range(n)]

def values(mu,j):
    return [a[j]+lam[j]*theta[i]+mu[i][j] for i in range(n)]

def truthful_tail(pi,t,mu):
    util=[F(0)]*n
    prices=[]
    winners=[]
    for j in pi[t:]:
        v=values(mu,j)
        w=max(range(n),key=lambda i:v[i])
        price=sorted(v)[-2]
        util[w]+=v[w]-price
        prices.append(price);winners.append(w)
        mu=update(mu,j,w)
    return util,prices,winners

state_checks=0
def verify_all_histories(pi,t,mu,leader):
    global state_checks
    if t==K:return
    j=pi[t];v=values(mu,j)
    assert all(v[leader]>v[i] for i in range(n) if i!=leader)
    util,_,winners=truthful_tail(pi,t,mu)
    assert all(w==leader for w in winners)
    for deviator in range(n):
        if deviator!=leader:
            after,_,_=truthful_tail(pi,t+1,update(mu,j,deviator))
            dev_util=v[deviator]-v[leader]+after[deviator]
            assert dev_util<util[deviator]
        else:
            for other in range(n):
                if other==leader:continue
                after,_,_=truthful_tail(pi,t+1,update(mu,j,other))
                assert after[leader]<util[leader]
    state_checks+=1
    for w in range(n):verify_all_histories(pi,t+1,update(mu,j,w),leader)

all_scores=[]
for suffix in permutations(range(1,K)):
    pi=(0,)+suffix
    total=W(pi)
    for leader in range(n):
        mu=[[F(0)]*K for _ in range(n)]
        mu=update(mu,0,leader)
        verify_all_histories(pi,1,mu,leader)
        tail_u,tail_p,tail_w=truthful_tail(pi,1,mu)
        T=max(theta[j] for j in range(n) if j!=leader)
        first_price=a[0]+lam[0]*T+delta*total
        winning_u=a[0]+lam[0]*theta[leader]-first_price+tail_u[leader]
        assert winning_u==sum(lam)*(theta[leader]-T)
    leader=max(range(n),key=lambda i:theta[i])
    T=sorted(theta)[-2]
    _,tail_p,_=truthful_tail(pi,1,update([[F(0)]*K for _ in range(n)],0,leader))
    first_price=a[0]+lam[0]*T+delta*total
    assert first_price+sum(tail_p)==sum(a)+sum(lam)*T+alpha_w*total
    all_scores.append((pi,total))

pi=(0,1,2,3,4)
rho_theta=F(n-1,n+1)
revenue=sum(a)+sum(lam)*rho_theta+alpha_w*W(pi)
other_revenue=revenue-alpha_w*B
first_anchor=a[0]+rho_theta+delta*W(pi)
first_terminal=a[1]+rho_theta
assert revenue==F(417,2)
assert other_revenue==F(369,2)

dp={0:F(0)}
items=list(range(1,K))
for mask in range(1,1<<len(items)):
    dp[mask]=max(dp[mask^(1<<q)]+sum(M[items[p]][items[q]] for p in range(len(items)) if p!=q and mask>>p&1) for q in range(len(items)) if mask>>q&1)
assert dp[(1<<len(items))-1]==max(score for _,score in all_scores)-(K-1)*B
best=max(all_scores,key=lambda z:z[1])
nonanchor_upper=sum(a)+sum(lam)+alpha_w*((K-2)*B+ordinary_total)
anchor_lower=sum(a)+sum(lam)*rho_theta+alpha_w*(K-1)*B
assert anchor_lower>nonanchor_upper

# This checks a pathwise identity for independent multidimensional types;
# arbitrary residual bids here are NOT claimed to constitute an equilibrium.
rng=Random(20260920)
general_checks=0
for _ in range(30):
    base=[[F(rng.randrange(101),100) for _ in range(K)] for _ in range(n)]
    residual=[F(rng.randrange(-400,501),100) for _ in range(n)]
    w=max(range(n),key=lambda i:residual[i])
    q_second=sorted(residual)[-2]
    competitor_cost=sum(max(base[i][k] for i in range(n) if i!=w) for k in range(1,K))
    phi=sum(base[w])-competitor_cost
    constant=q_second+competitor_cost
    for suffix in permutations(range(1,K)):
        order=(0,)+suffix
        first=q_second+delta*W(order)
        memory=update([[F(0)]*K for _ in range(n)],0,w)
        total_price=first
        utility=base[w][0]-first
        for j in suffix:
            v=[base[i][j]+memory[i][j] for i in range(n)]
            assert max(range(n),key=lambda i:v[i])==w
            price=sorted(v)[-2]
            total_price+=price
            utility+=v[w]-price
            memory=update(memory,j,w)
        assert total_price==constant+alpha_w*W(order)
        assert utility==phi-q_second
        general_checks+=1

print(json.dumps({
    'off_path_states_checked':state_checks,
    'suffix_permutations_checked':len(all_scores),
    'cue_strength':str(B),
    'best_suffix_gain':str(dp[(1<<len(items))-1]),
    'best_order':best[0],
    'best_W':str(best[1]),
    'anchor_first_expected_revenue':str(revenue),
    'terminal_first_expected_revenue':str(other_revenue),
    'anchor_first_round_expected_payment':str(first_anchor),
    'terminal_first_round_expected_payment':str(first_terminal),
    'anchor_any_order_revenue_lower_bound':str(anchor_lower),
    'nonanchor_any_equilibrium_revenue_upper_bound':str(nonanchor_upper),
    'independent_type_pathwise_identities_checked':general_checks,
},indent=2))
