"""气象过程块配对 bootstrap 与分层汇总。

重抽单位是调用方预先划定的整块：每个副本有放回抽取 B 个块（B=块总数），
块重复则块内全部配对随块一起重复；同一副本的块索引序列对所有系统、
reference、lead 与评分族完全一致，绝不按系统各自重抽。lead 分层内逐配对
等权：先在每个块内对该 lead 的配对取等权均值，副本再按块多重性对块均值
平均，因此配对数不同的块不会因为块内日期更多而获得更大权重。
skill=1-系统均分/参考均分在每个副本内部重算；参考均分为零或为负的
副本保留诊断、不参与该 skill 分位，且明确标注区间并非全部副本有效。
分位采用公开的 type-7 经验分位（与 numpy 默认 linear 同一定义）。
"""
import random

REASON_REF_ZERO = 'reference_mean_zero'
REASON_REF_NEGATIVE = 'reference_mean_negative'
REASON_NO_BLOCKS = 'no_blocks_drawn_for_lead'
REASON_BELOW_MIN_BLOCKS = 'effective_blocks_below_minimum'
REASON_NO_PAIRS = 'no_pairs_in_family'

METRICS = ('crps', 'brier')
FAMILIES = ('empirical', 'fair')


def percentile_type7(values, q):
    """type-7 经验分位：h=(n-1)q，在升序值间线性插值。"""
    if not values:
        raise ValueError('空样本无法计算分位')
    if not 0.0 <= q <= 1.0:
        raise ValueError('分位概率必须落在 [0,1]')
    xs = sorted(values)
    n = len(xs)
    if n == 1:
        return xs[0]
    h = (n - 1) * q
    lo = int(h)
    frac = h - lo
    if lo + 1 >= n:
        return xs[-1]
    return xs[lo] * (1.0 - frac) + xs[lo + 1] * frac


def _block_means(pairs, system, family, metric):
    """块内配对等权均值：block_id -> 该块该 lead 的等权均分。"""
    totals = {}
    counts = {}
    for p in pairs:
        entry = p['records'][system]
        value = entry[family][metric] if family == 'empirical' else entry['fair'][metric]
        totals[p['block_id']] = totals.get(p['block_id'], 0.0) + value
        counts[p['block_id']] = counts.get(p['block_id'], 0) + 1
    return {bid: totals[bid] / counts[bid] for bid in totals}


def _weighted_mean(pairs, weights, system, family, metric):
    """两级聚合：块内配对等权，副本内按块多重性对块均值平均。"""
    block_mean = _block_means(pairs, system, family, metric)
    total = 0.0
    weight_sum = 0
    for bid, mean in block_mean.items():
        w = weights.get(bid, 0)
        if not w:
            continue
        total += w * mean
        weight_sum += w
    if weight_sum == 0:
        return None, 0
    return total / weight_sum, weight_sum


def _comparisons(systems):
    out = []
    for i in range(len(systems)):
        for j in range(i + 1, len(systems)):
            out.append((f'{systems[j]}_minus_{systems[i]}', systems[i], systems[j]))
    return out


def _family_aggregate(pairs, weights, all_systems, family):
    out = {}
    for metric in METRICS:
        means = {}
        wsum_seen = 0
        for s in all_systems:
            mean, wsum = _weighted_mean(pairs, weights, s, family, metric)
            means[s] = mean
            wsum_seen = max(wsum_seen, wsum)
        out[metric] = None if not wsum_seen else {'means': means}
    return out


def run_bootstrap(pairing, *, n_draws, seed, percentiles, min_blocks):
    if type(n_draws) is not int or n_draws < 1:
        raise ValueError('bootstrap_draws 必须为正整数')
    if type(min_blocks) is not int or min_blocks < 1:
        raise ValueError('min_blocks 必须为 >=1 的整数')
    if not isinstance(percentiles, list) or not percentiles:
        raise ValueError('percentiles 必须为非空列表')
    qs = [float(q) for q in percentiles]
    for q in qs:
        if not 0.0 <= q <= 1.0:
            raise ValueError('分位概率必须落在 [0,1]')
    if type(seed) is not int:
        raise ValueError('seed 必须为整数')

    systems = pairing['systems']
    reference = pairing['reference']
    all_systems = pairing['all_systems']
    comparisons = _comparisons(systems)
    block_universe = [b['id'] for b in pairing['blocks']]
    n_blocks = len(block_universe)
    index_of = {bid: i for i, bid in enumerate(block_universe)}

    rng = random.Random(seed)
    draws = []
    for _ in range(n_draws):
        picked = [rng.randrange(n_blocks) for _ in range(n_blocks)] if n_blocks else []
        multiplicity = {}
        for idx in picked:
            bid = block_universe[idx]
            multiplicity[bid] = multiplicity.get(bid, 0) + 1
        draws.append(multiplicity)

    leads = sorted({p['lead'] for p in pairing['empirical_pairs']})
    result = {'leads': {}}
    for lead in leads:
        lead_key = str(lead)
        result['leads'][lead_key] = {}
        for family in FAMILIES:
            source = pairing['empirical_pairs'] if family == 'empirical' else pairing['fair_pairs']
            pairs = [p for p in source if p['lead'] == lead]
            section = _family_section(
                pairs, family, draws, all_systems, systems, reference,
                comparisons, block_universe, index_of, qs, min_blocks, n_draws)
            result['leads'][lead_key][family] = section

    result['bootstrap'] = {
        'n_draws': n_draws,
        'seed': seed,
        'percentiles': qs,
        'min_blocks': min_blocks,
        'block_universe': block_universe,
        'n_blocks': n_blocks,
        'draws': [{'replicate': k, 'multiplicity': m} for k, m in enumerate(draws)],
        'resampling_unit': 'whole weather-process block, with replacement, paired across all systems',
    }
    return result


def _pair_equal_means(pairs, all_systems, family):
    """点估计口径：lead 分层内全部配对逐配对等权，不经过块大小调整。"""
    out = {}
    for metric in METRICS:
        means = {}
        for s in all_systems:
            vals = []
            for p in pairs:
                entry = p['records'][s]
                vals.append(entry[family][metric]
                            if family == 'empirical' else entry['fair'][metric])
            means[s] = sum(vals) / len(vals)
        out[metric] = {'means': means}
    return out


def _family_section(pairs, family, draws, all_systems, systems, reference,
                    comparisons, block_universe, index_of, qs, min_blocks, n_draws):
    if not pairs:
        return {'available': False, 'reason': REASON_NO_PAIRS,
                'detail': '该 lead 没有满足 m>=2 的共同配对，fair 评分不可估'
                if family == 'fair' else '该 lead 无共同配对'}
    present_blocks = sorted({p['block_id'] for p in pairs},
                            key=lambda bid: index_of[bid])
    n_effective_blocks = len(present_blocks)

    # 点估计：分层内逐配对等权（与块大小无关）
    point_agg = _pair_equal_means(pairs, all_systems, family)
    point_means = {m: point_agg[m]['means'] for m in METRICS}
    point = {'means': point_means, 'deltas': {}, 'skills': {}}
    for metric in METRICS:
        means = point_means[metric]
        point['deltas'][metric] = {
            label: means[j] - means[i] for label, i, j in comparisons}
        point['skills'][metric] = {}
        for s in systems:
            ref = means[reference]
            if ref == 0:
                point['skills'][metric][s] = {'value': None, 'valid': False,
                                              'reason': REASON_REF_ZERO}
            elif ref < 0:
                point['skills'][metric][s] = {'value': None, 'valid': False,
                                              'reason': REASON_REF_NEGATIVE}
            else:
                point['skills'][metric][s] = {'value': 1.0 - means[s] / ref,
                                              'valid': True}

    section = {
        'available': True,
        'n_pairs': len(pairs),
        'n_effective_blocks': len(present_blocks),
        'effective_blocks': present_blocks,
        'point_weighting': 'lead 分层内逐配对等权（不按块大小调整）',
        'replicate_weighting': '块内配对等权求块均值，副本按块多重性对块均值平均',
        'point': point,
    }
    if n_effective_blocks < min_blocks:
        section['bootstrap'] = {
            'status': 'point_only',
            'reason': REASON_BELOW_MIN_BLOCKS,
            'detail': f'有效块 {n_effective_blocks} < min_blocks {min_blocks}，仅报点估计',
        }
        return section

    # 副本收集器
    rep_means = {m: {s: [] for s in all_systems} for m in METRICS}
    rep_deltas = {m: {label: [] for label, _, _ in comparisons} for m in METRICS}
    rep_skills = {m: {s: [] for s in systems} for m in METRICS}
    skill_invalid = {m: {s: {} for s in systems} for m in METRICS}
    stratum_invalid = {}
    for rep_id, multiplicity in enumerate(draws):
        agg = _family_aggregate(pairs, multiplicity, all_systems, family)
        if any(agg[m] is None for m in METRICS):
            stratum_invalid[rep_id] = REASON_NO_BLOCKS
            continue
        for metric in METRICS:
            means = agg[metric]['means']
            for s in all_systems:
                rep_means[metric][s].append(means[s])
            for label, i, j in comparisons:
                rep_deltas[metric][label].append(means[j] - means[i])
            for s in systems:
                ref = means[reference]
                if ref == 0 or ref < 0:
                    reason = REASON_REF_ZERO if ref == 0 else REASON_REF_NEGATIVE
                    skill_invalid[metric][s][reason] = \
                        skill_invalid[metric][s].get(reason, 0) + 1
                else:
                    rep_skills[metric][s].append(1.0 - means[s] / ref)

    valid_replicates = n_draws - len(stratum_invalid)
    intervals = {'means': {}, 'deltas': {}, 'skills': {}}
    for metric in METRICS:
        intervals['means'][metric] = {
            s: {f'p{q:g}': percentile_type7(rep_means[metric][s], q) for q in qs}
            for s in all_systems}
        intervals['deltas'][metric] = {
            label: {f'p{q:g}': percentile_type7(rep_deltas[metric][label], q)
                    for q in qs}
            for label, _, _ in comparisons}
        intervals['skills'][metric] = {}
        for s in systems:
            values = rep_skills[metric][s]
            if not values:
                no_block_invalid = ({REASON_NO_BLOCKS: len(stratum_invalid)}
                                    if stratum_invalid else {})
                intervals['skills'][metric][s] = {
                    'percentiles': None, 'interval_valid': False,
                    'replicate_counts': {'total': n_draws, 'valid': 0,
                                         'invalid': n_draws},
                    'invalid_reasons': {**skill_invalid[metric][s],
                                        **no_block_invalid}}
                continue
            counts = {'total': n_draws, 'valid': len(values),
                      'invalid': n_draws - len(values)}
            no_block_invalid = ({REASON_NO_BLOCKS: len(stratum_invalid)}
                                if stratum_invalid else {})
            intervals['skills'][metric][s] = {
                'percentiles': {f'p{q:g}': percentile_type7(values, q) for q in qs},
                'interval_valid': counts['invalid'] == 0,
                'replicate_counts': counts,
                'invalid_reasons': {**skill_invalid[metric][s],
                                    **no_block_invalid}}
    section['bootstrap'] = {
        'status': 'ok',
        'replicate_counts': {'total': n_draws, 'valid': valid_replicates,
                             'invalid': len(stratum_invalid)},
        'invalid_reasons': ({REASON_NO_BLOCKS: len(stratum_invalid)}
                            if stratum_invalid else {}),
        'intervals': intervals,
    }
    return section
