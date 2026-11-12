"""按冻结 as_of 规则选择各系统最新版本并构造严格配对。

配对身份为 (valid_day, lead)。每个系统在配对内只选 issued_day 与 recorded_day
均不晚于 valid_day - as_of_lag_days 的记录，先取最大 issued_day，再取该签发日
下的最大 revision；成员不得跨版本拼接。任一比较系统、reference 或观测不完整，
整个配对排除并记录原因；经验与 fair 各自有明确的共同样本，不隐式混用。
"""
from collections import defaultdict
from .io import date, integer, number, unique
from .scores import empirical_brier, empirical_crps, fair_brier, fair_crps

# 排除原因代码（英文稳定标识，报告同时给中文说明）
REASON_MISSING_OBSERVATION = 'missing_observation'
REASON_NO_VISIBLE_RECORD = 'no_visible_record'
REASON_LATEST_INCOMPLETE = 'latest_record_incomplete'
REASON_PAYLOAD_CONFLICT = 'same_issue_revision_payload_conflict'
REASON_STRUCTURAL = 'structural_violation'

_REASON_TEXT = {
    REASON_MISSING_OBSERVATION: '该验证日无观测或观测值为 null',
    REASON_NO_VISIBLE_RECORD: '该系统在 as_of 冻结点内无可见记录',
    REASON_LATEST_INCOMPLETE: '该系统最新可见记录成员为空或含 null',
    REASON_PAYLOAD_CONFLICT: '同签发日同 revision 存在载荷不同的记录，拒绝拼接',
    REASON_STRUCTURAL: '资料结构违例（见报错信息），无法配对',
}


class StructuralError(ValueError):
    """资料违反配对契约，属整体拒收而非单配对排除。"""


def reason_text(code):
    return _REASON_TEXT[code]


def _validate_forecast(r):
    issued = date(r['issued_day'])
    verified = date(r['valid_day'])
    lead = integer(r['lead'], 1)
    if (verified - issued).days != lead:
        raise StructuralError('valid_day-issued_day 必须等于 lead')
    if not isinstance(r['system'], str) or not r['system']:
        raise StructuralError('系统身份缺失')
    integer(r['revision'], 1)
    date(r['recorded_day'])
    if not isinstance(r['members'], list):
        raise StructuralError('成员应为列表')
    for x in r['members']:
        if x is not None:
            number(x)
    return issued, verified, lead


def _record_complete(r):
    return bool(r['members']) and all(x is not None for x in r['members'])


def _validate_blocks(raw_blocks, common_days):
    blocks = []
    unique(raw_blocks, 'id')
    for b in raw_blocks:
        start = date(b['start'])
        end = date(b['end'])
        if end < start:
            raise StructuralError('块结束日早于开始日')
        blocks.append({'id': b['id'], 'start': start, 'end': end})
    blocks.sort(key=lambda b: b['start'])
    for prev, nxt in zip(blocks, blocks[1:]):
        if nxt['start'] <= prev['end']:
            raise StructuralError('气象过程块彼此重叠')
    if not common_days:
        return blocks, {}
    day_to_block = {}
    for day in common_days:
        hits = [b['id'] for b in blocks if b['start'] <= day <= b['end']]
        if not hits:
            raise StructuralError(f'共同验证日 {day.isoformat()} 未被任何块覆盖')
        day_to_block[day] = hits[0]
    unused = [b['id'] for b in blocks
              if not any(b['start'] <= day <= b['end'] for day in common_days)]
    if unused:
        raise StructuralError(f'块未覆盖任何共同验证日期: {unused}')
    return blocks, day_to_block


def build_pairs(payload):
    """返回配对结果字典；结构违例抛 StructuralError。"""
    try:
        systems = payload['systems']
        reference = payload['reference']
        lag = integer(payload['as_of_lag_days'], 0)
        threshold = number(payload['event_threshold'])
        raw_blocks = payload['blocks']
    except (KeyError, TypeError) as e:
        raise StructuralError('顶层字段缺失或类型错误') from e
    if not isinstance(systems, list) or not systems or any(
            not isinstance(s, str) or not s for s in systems):
        raise StructuralError('systems 应为非空非空字符串列表')
    if len(set(systems)) != len(systems):
        raise StructuralError('systems 存在重复系统')
    if not isinstance(reference, str) or not reference:
        raise StructuralError('reference 身份缺失')
    if reference in systems:
        raise StructuralError('reference 不得同时列入比较系统')
    all_systems = list(systems) + [reference]

    unique(payload['forecasts'], 'record_id')
    obs = {}
    obs_days = set()
    for o in payload['observations']:
        day = date(o['day'])
        if day in obs_days:
            raise StructuralError('验证观测日期重复')
        obs_days.add(day)
        obs[day] = o['value']
        if o['value'] is not None:
            number(o['value'])

    # key=(valid_day, lead) -> system -> 可见记录列表
    visible = defaultdict(lambda: defaultdict(list))
    forecast_systems = set()
    for r in payload['forecasts']:
        issued, verified, lead_value = _validate_forecast(r)
        recorded = date(r['recorded_day'])
        forecast_systems.add(r['system'])
        cutoff = verified.fromordinal(verified.toordinal() - lag)
        if issued <= cutoff and recorded <= cutoff:
            visible[(verified, lead_value)][r['system']].append(r)
    unknown = [s for s in all_systems if s not in forecast_systems]
    if unknown:
        raise StructuralError(f'系统在预报记录中完全不存在: {unknown}')

    candidate_keys = sorted(visible.keys(), key=lambda k: (k[0], k[1]))
    exclusions = []
    selected_rows = []
    for key in candidate_keys:
        verified, lead_value = key
        pair_reasons = []
        if verified not in obs or obs[verified] is None:
            pair_reasons.append({'system': None,
                                 'reason': REASON_MISSING_OBSERVATION,
                                 'detail': reason_text(REASON_MISSING_OBSERVATION)})
        chosen = {}
        for s in all_systems:
            rows = visible[key].get(s, [])
            if not rows:
                pair_reasons.append({'system': s, 'reason': REASON_NO_VISIBLE_RECORD,
                                     'detail': reason_text(REASON_NO_VISIBLE_RECORD)})
                continue
            latest_issue = max(date(r['issued_day']) for r in rows)
            latest_rows = [r for r in rows if date(r['issued_day']) == latest_issue]
            latest_revision = max(integer(r['revision']) for r in latest_rows)
            finalists = [r for r in latest_rows
                         if integer(r['revision']) == latest_revision]
            payloads = {tuple(r['members']) for r in finalists}
            if len(payloads) > 1:
                pair_reasons.append({'system': s, 'reason': REASON_PAYLOAD_CONFLICT,
                                     'detail': reason_text(REASON_PAYLOAD_CONFLICT)})
                continue
            pick = finalists[0]
            if not _record_complete(pick):
                pair_reasons.append({'system': s, 'reason': REASON_LATEST_INCOMPLETE,
                                     'detail': reason_text(REASON_LATEST_INCOMPLETE)})
                continue
            chosen[s] = pick
        if pair_reasons:
            exclusions.append({
                'valid_day': verified.isoformat(),
                'lead': lead_value,
                'reasons': pair_reasons,
            })
            continue
        observation = obs[verified]
        row = {'valid_day': verified.isoformat(), 'lead': lead_value,
               'observation': observation, 'threshold': threshold,
               'records': {}}
        fair_eligible = True
        for s in all_systems:
            pick = chosen[s]
            members = list(pick['members'])
            crps_e = empirical_crps(members, observation)
            brier_e, p, o_event = empirical_brier(members, observation, threshold)
            entry = {'record_id': pick['record_id'], 'issued_day': pick['issued_day'],
                     'recorded_day': pick['recorded_day'], 'revision': integer(pick['revision']),
                     'm': len(members), 'members': members,
                     'empirical': {'crps': crps_e, 'brier': brier_e,
                                   'event_probability': p}}
            if len(members) < 2:
                fair_eligible = False
                entry['fair'] = None
            else:
                crps_f = fair_crps(members, observation)
                brier_f, _, _ = fair_brier(members, observation, threshold)
                entry['fair'] = {'crps': crps_f, 'brier': brier_f}
            row['records'][s] = entry
        row['fair_eligible'] = fair_eligible
        selected_rows.append(row)

    common_days = sorted({date(r['valid_day']) for r in selected_rows})
    blocks, day_to_block = _validate_blocks(raw_blocks, common_days)
    for r in selected_rows:
        r['block_id'] = day_to_block[date(r['valid_day'])]

    empirical_pairs = selected_rows
    fair_pairs = [r for r in selected_rows if r['fair_eligible']]
    return {
        'systems': systems,
        'reference': reference,
        'all_systems': all_systems,
        'threshold': threshold,
        'as_of_lag_days': lag,
        'blocks': [{'id': b['id'], 'start': b['start'].isoformat(),
                    'end': b['end'].isoformat()} for b in blocks],
        'exclusions': exclusions,
        'reason_legend': dict(_REASON_TEXT),
        'empirical_pairs': empirical_pairs,
        'fair_pairs': fair_pairs,
        'counts': {
            'candidate_pairs': len(candidate_keys),
            'excluded_pairs': len(exclusions),
            'empirical_common_pairs': len(empirical_pairs),
            'fair_common_pairs': len(fair_pairs),
        },
    }
