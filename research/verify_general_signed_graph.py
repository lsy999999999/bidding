"""Exact rational checks for arbitrary signed graphs; NOT an exact-equilibrium solver.

Opponents use the common-shift benchmark. Compute the deviator's optimal
history-dependent response for independent Bernoulli private item values.
No discretization of bids is used: the positive-part Bellman recursion solves
the optimal threshold bid against the two possible rival maximum bids.
"""
from fractions import Fraction as F
from itertools import permutations, product
from functools import lru_cache
from math import comb
from pathlib import Path
import json


def subset_order_dp(weights):
    k = len(weights)
    score = {0: F(0)}
    order = {0: ()}
    for mask in range(1, 1 << k):
        options = []
        for last in range(k):
            if mask >> last & 1:
                prev = mask ^ (1 << last)
                incoming = sum((weights[j][last] for j in range(k)
                                if prev >> j & 1), F(0))
                options.append((score[prev] + incoming, order[prev] + (last,)))
        score[mask], order[mask] = max(options)
    return score[(1 << k) - 1], order[(1 << k) - 1]


def forward_sum(matrix, order):
    return sum((matrix[j][k] for t, k in enumerate(order)
                for j in order[:t]), F(0))


def pair_utility(matrix, order, owned_mask):
    return sum((matrix[j][k] for t, k in enumerate(order)
                for j in order[:t]
                if (owned_mask >> j & 1) and (owned_mask >> k & 1)), F(0))


def response_values(n, matrix, order, own_type, delta):
    """Interim utility of optimal response and of benchmark, exact fractions.

    Common item offsets cancel from value minus rival price. Opponent bids
    ignore history; independent items make earlier observations uninformative
    about later rival types. Thus (time, own holdings) is sufficient here.
    This sufficiency is NOT asserted for opponents' equilibrium strategies.
    """
    k = len(order)
    p0 = F(1, 2 ** (n - 1))
    max_probs = ((0, p0), (1, 1 - p0))
    tie_one = sum((F(comb(n - 1, j), 2 ** (n - 1)) / (j + 1)
                   for j in range(1, n)), F(0)) / (1 - p0)

    def win_probability(value, rival_max):
        if value > rival_max:
            return F(1)
        if value < rival_max:
            return F(0)
        return F(1, n) if value == 0 else tie_one

    @lru_cache(None)
    def best(t, owned):
        if t == k:
            return F(0)
        item = order[t]
        lose = best(t + 1, owned)
        win = best(t + 1, owned | (1 << item))
        shift = delta * sum((matrix[j][item] for j in range(k)
                             if owned >> j & 1), F(0))
        threshold = own_type[item] + shift + win - lose
        return lose + sum((prob * max(F(0), threshold - rival_max)
                           for rival_max, prob in max_probs), F(0))

    @lru_cache(None)
    def benchmark(t, owned):
        if t == k:
            return F(0)
        item = order[t]
        lose = benchmark(t + 1, owned)
        win = benchmark(t + 1, owned | (1 << item))
        shift = delta * sum((matrix[j][item] for j in range(k)
                             if owned >> j & 1), F(0))
        out = F(0)
        for rival_max, prob in max_probs:
            wp = win_probability(own_type[item], rival_max)
            out += prob * (wp * (own_type[item] - rival_max + shift + win)
                           + (1 - wp) * lose)
        return out

    return best(0, 0), benchmark(0, 0)


def as_json(value):
    if isinstance(value, F):
        return {"exact": str(value), "decimal": float(value)}
    if isinstance(value, dict):
        return {str(k): as_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [as_json(v) for v in value]
    return value


def run():
    n, k, alpha_l, delta, base = 3, 4, F(1), F(2, 5), F(3)
    # A->B positive, B->C negative, C->A positive closes a cycle.
    # C->D positive, D->B negative; A->B->D->A is an all-positive cycle.
    # A->D, B->A, D->C are exactly zero. Every node has an outgoing edge.
    m = [list(map(F, row)) for row in (
        (0, F(3, 10), F(-1, 10), 0),
        (0, 0, F(-2, 5), F(1, 10)),
        (F(1, 5), F(1, 20), 0, F(1, 4)),
        (F(3, 10), F(-3, 20), 0, 0),
    )]
    alpha_w = alpha_l + delta
    for item in range(k):
        assert base > alpha_w * sum((max(F(0), -m[j][item])
                                     for j in range(k)), F(0))
    common = [[alpha_l * w for w in row] for row in m]
    upper = [[alpha_l * w + delta * max(w, F(0)) for w in row] for row in m]
    g_opt, pi_opt = subset_order_dp(common)
    upper_opt, _ = subset_order_dp(upper)
    types = tuple(product((0, 1), repeat=k))
    mean_top = 1 - F(1, 2 ** n)
    mean_second = 1 - F(1 + n, 2 ** n)
    r0 = k * (base + mean_second)
    rent0 = k * (mean_top - mean_second)
    records = []
    for pi in permutations(range(k)):
        w = forward_sum(m, pi)
        potentials = [pair_utility(m, pi, mask) for mask in range(1 << k)]
        eps_range = delta * (max(potentials) - min(potentials))
        eps_abs = delta * forward_sum([[abs(v) for v in row] for row in m], pi)
        response = [response_values(n, m, pi, typ, delta) for typ in types]
        regrets = [br - bench for br, bench in response]
        assert min(regrets) >= 0
        assert max(regrets) <= eps_range <= eps_abs
        total_utility = n * sum((bench for _, bench in response), F(0)) / len(types)
        assert total_utility == rent0 + delta * w / n
        a_correction = k * (base + mean_top) - total_utility
        l_correction = (1 - F(1, n)) * w
        revenue = r0 + alpha_l * w
        assert revenue == alpha_w * w - delta * l_correction + a_correction
        records.append({
            "order": ''.join('ABCD'[j] for j in pi),
            "W": w, "benchmark_revenue": revenue,
            "max_interim_dynamic_regret": max(regrets),
            "mean_dynamic_regret": sum(regrets, F(0)) / len(regrets),
            "epsilon_subset_range": eps_range, "epsilon_absolute_edges": eps_abs,
            "total_buyer_utility": total_utility,
            "A_correction": a_correction, "expected_L": l_correction,
        })
    assert g_opt == max(forward_sum(common, p) for p in permutations(range(k)))
    assert upper_opt == max(forward_sum(upper, p) for p in permutations(range(k)))
    chosen = ''.join('ABCD'[j] for j in pi_opt)
    out = {
        "scope": "Exact rational benchmark and best responses, not a computed Nash equilibrium",
        "n": n, "K": k, "alpha_L": alpha_l, "alpha_W": alpha_w,
        "base_value": base, "noise_distribution": "independent Bernoulli(1/2)",
        "M": m, "best_benchmark_order": chosen,
        "chosen_result": next(r for r in records if r['order'] == chosen),
        "equilibrium_revenue_upper_bound_all_orders": k * (base + mean_top) + upper_opt,
        "upper_minus_chosen_benchmark_revenue": rent0 + upper_opt - g_opt,
        "all_orders": records,
    }
    # δ=0 sanity check: the same Bellman deviations must give zero gain.
    for typ in types:
        br, bench = response_values(n, m, pi_opt, typ, F(0))
        assert br == bench
    return out


if __name__ == '__main__':
    result = as_json(run())
    path = Path(__file__).with_name('general_signed_graph_results.json')
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps({k: v for k, v in result.items() if k not in ('all_orders', 'M')},
                     ensure_ascii=False, indent=2))
