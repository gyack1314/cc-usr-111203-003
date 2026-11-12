"""配对检验门面：读取资料 -> 严格配对 -> 分布评分 -> 过程块配对 bootstrap。

报告先声明比较的是哪一种分布评分（经验 CRPS / fair CRPS / 经验 Brier /
fair Brier，均为越低越好的误差型评分），再给出配对排除原因、经验与 fair
各自的共同样本、原始分数，以及按气象过程块的配对 bootstrap 区间。
块由调用方给定并在看到任何评分差异之前冻结，本模块不根据差异挑选块。
"""
from .io import integer, load
from .pairing import build_pairs
from .bootstrap import run_bootstrap

SCORE_DEFINITIONS = {
    'empirical_crps': {
        'family': 'empirical',
        'metric': 'crps',
        'direction': 'lower_is_better',
        'formula': 'mean_i|x_i-y| - (1/(2 m^2)) * sum_{i,j}|x_i-x_j|，等权成员经验分布 CRPS',
    },
    'fair_crps': {
        'family': 'fair',
        'metric': 'crps',
        'direction': 'lower_is_better',
        'formula': 'mean_i|x_i-y| - (1/(2 m(m-1))) * sum_{i,j}|x_i-x_j|，仅 m>=2 可估',
        'assumption': '成员为同一潜在预报分布的独立同分布抽样；确定性或强相关成员不自动满足',
    },
    'empirical_brier': {
        'family': 'empirical',
        'metric': 'brier',
        'direction': 'lower_is_better',
        'event': '观测严格大于锁定阈值 threshold，即 y > threshold',
        'formula': '(p-o)^2，p=成员经验事件频率，o=事件示性 (0/1)',
    },
    'fair_brier': {
        'family': 'fair',
        'metric': 'brier',
        'direction': 'lower_is_better',
        'event': '观测严格大于锁定阈值 threshold，即 y > threshold',
        'formula': '(p-o)^2 - p(1-p)/(m-1)，仅 m>=2 可估，允许为负、不裁零',
        'assumption': '成员为同一潜在预报分布的独立同分布抽样；确定性或强相关成员不自动满足',
    },
}

COMPARISON_CONVENTION = {
    'delta': 'new_minus_old 方向按 systems 列表顺序定义：systems[j] 均分 - systems[i] 均分；负值表示后者更优',
    'skill': 'skill = 1 - 系统均分 / reference 均分，在每个 bootstrap 副本内重算',
    'pairing': '配对身份为 (valid_day, lead)；任一系统、reference 或观测不完整则整对排除，禁止按系统各自算完再相减',
    'bootstrap': '整块有放回重抽，块重复则块内全部配对重复；所有系统与 reference 使用同一份块副本',
}


def run_verification(path):
    """从资料路径运行完整检验，返回可直接 JSON 序列化的报告字典。"""
    payload, digest = load(path)
    pairing = build_pairs(payload)
    bootstrap = run_bootstrap(
        pairing,
        n_draws=integer(payload['bootstrap_draws'], 1),
        seed=integer(payload['seed']),
        percentiles=list(payload['percentiles']),
        min_blocks=integer(payload['min_blocks'], 1),
    )
    report = {
        'source': {'path': str(path), 'sha256': digest},
        'score_definitions': SCORE_DEFINITIONS,
        'comparison_convention': COMPARISON_CONVENTION,
        'event': {'threshold': pairing['threshold'],
                  'rule': 'observation > threshold (strictly greater)'},
        'as_of_lag_days': pairing['as_of_lag_days'],
        'systems': pairing['systems'],
        'reference': pairing['reference'],
        'blocks': pairing['blocks'],
        'pairing': {
            'counts': pairing['counts'],
            'reason_legend': pairing['reason_legend'],
            'exclusions': pairing['exclusions'],
            'empirical_common_sample': [
                {'valid_day': p['valid_day'], 'lead': p['lead'],
                 'block_id': p['block_id'], 'observation': p['observation']}
                for p in pairing['empirical_pairs']],
            'fair_common_sample': [
                {'valid_day': p['valid_day'], 'lead': p['lead'],
                 'block_id': p['block_id'], 'observation': p['observation']}
                for p in pairing['fair_pairs']],
        },
        'raw_scores': _raw_scores(pairing),
        'results': bootstrap,
    }
    return report


def _raw_scores(pairing):
    out = []
    for p in pairing['empirical_pairs']:
        row = {'valid_day': p['valid_day'], 'lead': p['lead'],
               'block_id': p['block_id'], 'observation': p['observation'],
               'fair_eligible': p['fair_eligible'], 'systems': {}}
        for s, entry in p['records'].items():
            row['systems'][s] = {
                'record_id': entry['record_id'], 'issued_day': entry['issued_day'],
                'recorded_day': entry['recorded_day'], 'revision': entry['revision'],
                'm': entry['m'],
                'empirical': entry['empirical'], 'fair': entry['fair']}
        out.append(row)
    return out
