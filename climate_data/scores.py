"""等权集合成员的分布评分原语。

事件口径：观测或成员值 *严格大于* 锁定阈值 threshold（y > threshold）。

CRPS（连续排位概率评分，越低越好），m 为成员数：
- 经验 CRPS = mean_i |x_i - y| - 有序成员对距离和 / (2 m^2)
  其中有序对距离和等于 2 * 无序对距离和，故第二项实现为 sum_{i<j}|x_i-x_j| / m^2。
- fair CRPS 仅把第二项分母改为 2 m (m-1)，即 sum_{i<j}|x_i-x_j| / (m(m-1))，
  仅 m>=2 可估；m=1 只有经验分，禁止把经验分混入 fair 均值。

Brier（阈值事件概率误差平方，越低越好），p=经验事件频率，o=事件示性：
- 经验 Brier = (p-o)^2
- fair Brier = (p-o)^2 - p(1-p)/(m-1)，仅 m>=2 可估，允许为负、不裁零。
fair 修正假设成员是同一潜在预报分布的独立同分布抽样；确定性或强相关集合不自动满足。
"""
from .io import number


class FairNotEstimable(ValueError):
    """成员数 m<2 时 fair 评分不可估计。"""


def _clean(members):
    if not isinstance(members, list) or not members or any(x is None for x in members):
        raise ValueError('成员列表必须为非空列表且不含 null')
    return [number(x) for x in members]


def _pair_distance_sum(members):
    """无序成员对 i<j 的绝对距离和；等于全部有序对距离和的一半。"""
    return sum(abs(xi - xj) for i, xi in enumerate(members) for xj in members[i + 1:])


def empirical_crps(members, observation):
    members = _clean(members)
    y = number(observation)
    m = len(members)
    mae = sum(abs(x - y) for x in members) / m
    spread = _pair_distance_sum(members) / (m * m)
    return mae - spread


def fair_crps(members, observation):
    members = _clean(members)
    y = number(observation)
    m = len(members)
    if m < 2:
        raise FairNotEstimable('fair CRPS 要求 m>=2，单成员不可估')
    mae = sum(abs(x - y) for x in members) / m
    spread = _pair_distance_sum(members) / (m * (m - 1))
    return mae - spread


def event_indicator(value, threshold):
    """严格大于阈值记为 1。"""
    return 1 if number(value) > number(threshold) else 0


def empirical_probability(members, threshold):
    members = _clean(members)
    return sum(event_indicator(x, threshold) for x in members) / len(members)


def empirical_brier(members, observation, threshold):
    """返回 (经验 Brier, 经验概率 p, 事件示性 o)。"""
    members = _clean(members)
    p = empirical_probability(members, threshold)
    o = event_indicator(observation, threshold)
    return (p - o) ** 2, p, o


def fair_brier_value(p, o, m):
    """fair Brier 公式原形，不做裁零；仅用于口径检验时 p 可与计数频率不一致。"""
    if type(m) is not int or m < 2:
        raise FairNotEstimable('fair Brier 要求 m>=2，单成员不可估')
    p = number(p)
    o = number(o)
    return (p - o) ** 2 - p * (1.0 - p) / (m - 1)


def fair_brier(members, observation, threshold):
    """返回 (fair Brier, 经验概率 p, 事件示性 o)；m<2 抛 FairNotEstimable。"""
    members = _clean(members)
    m = len(members)
    if m < 2:
        raise FairNotEstimable('fair Brier 要求 m>=2，单成员不可估')
    p = empirical_probability(members, threshold)
    o = event_indicator(observation, threshold)
    return fair_brier_value(p, o, m), p, o
