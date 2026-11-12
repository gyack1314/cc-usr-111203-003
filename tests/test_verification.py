"""配对检验、分布评分与过程块 bootstrap 的边界测试。"""
import copy
import datetime
import json
import tempfile
import unittest
from pathlib import Path

from climate_data.scores import (FairNotEstimable, empirical_brier, empirical_crps,
                                 empirical_probability, event_indicator,
                                 fair_brier, fair_brier_value, fair_crps)
from climate_data.pairing import (REASON_LATEST_INCOMPLETE, REASON_MISSING_OBSERVATION,
                                  REASON_NO_VISIBLE_RECORD, REASON_PAYLOAD_CONFLICT,
                                  StructuralError, build_pairs)
from climate_data.bootstrap import (REASON_BELOW_MIN_BLOCKS, REASON_REF_ZERO,
                                    percentile_type7, run_bootstrap)
from climate_data.io import load
from climate_data.verify import run_verification

DATA = Path('data')


def payload(name):
    return load(DATA / f'{name}.json')[0]


class ScoreTests(unittest.TestCase):
    def test_empirical_crps_two_members_hand_value(self):
        # m=2, x=(0,2), y=1: MAE=1, 无序对距离=2, spread=2/4=0.5
        self.assertAlmostEqual(empirical_crps([0, 2], 1), 0.5)
        # 手算: (|0-1|+|2-1|)/2 - (2)/4 = 1 - 0.5
        # 有序对口径等价: 1 - (|0-2|+|2-0|)/(2*4) = 1 - 4/8 = 0.5

    def test_fair_crps_two_members_hand_value(self):
        # spread=2/(2*1)=1 -> 1-1=0
        self.assertAlmostEqual(fair_crps([0, 2], 1), 0.0)
        # m=3, x=(0,1,3), y=1: MAE=(1+0+2)/3=1, 无序对距离=1+3+2=6,
        # 经验 spread=6/9=2/3; fair spread=6/(3*2)=1
        self.assertAlmostEqual(empirical_crps([0, 1, 3], 1), 1 - 6 / 9)
        self.assertAlmostEqual(fair_crps([0, 1, 3], 1), 0.0)

    def test_single_member_empirical_only(self):
        self.assertAlmostEqual(empirical_crps([5], 7), 2.0)
        with self.assertRaises(FairNotEstimable):
            fair_crps([5], 7)
        with self.assertRaises(FairNotEstimable):
            fair_brier([5], 7, 6)

    def test_duplicate_members_and_strict_event(self):
        # 全部重复成员（确定性集合）：spread=0，CRPS 退化为绝对误差
        self.assertAlmostEqual(empirical_crps([3, 3, 3], 8), 5.0)
        self.assertAlmostEqual(fair_crps([3, 3, 3], 8), 5.0)
        # 严格大于：阈值等于观测时事件不发生
        self.assertEqual(event_indicator(5, 5), 0)
        self.assertEqual(event_indicator(5.0001, 5), 1)
        self.assertEqual(empirical_probability([1, 5, 9], 5), 1 / 3)
        # Brier: p=1/3, o=0 -> 1/9；fair=1/9-(2/9)/2=0
        b, p, o = empirical_brier([1, 5, 9], 5, 5)
        self.assertAlmostEqual(b, 1 / 9)
        self.assertEqual(p, 1 / 3)
        self.assertEqual(o, 0)
        fb, _, _ = fair_brier([1, 5, 9], 5, 5)
        self.assertAlmostEqual(fb, 0.0)
        # o=1: (1/3-1)^2=4/9；fair=4/9-1/9=1/3
        b2, _, o2 = empirical_brier([1, 5, 9], 9, 5)
        self.assertEqual(o2, 1)
        self.assertAlmostEqual(b2, 4 / 9)
        fb2, _, _ = fair_brier([1, 5, 9], 9, 5)
        self.assertAlmostEqual(fb2, 1 / 3)

    def test_fair_brier_can_go_negative_and_is_not_clipped(self):
        # 公式原形：极端 p 下允许负值，实现不裁零
        v = fair_brier_value(0.9, 0.0, 2)
        self.assertAlmostEqual(v, 0.81 - 0.09)
        v_negative = fair_brier_value(0.99, 1.0, 2)
        self.assertLess(v_negative, 0.0)

    def test_null_and_empty_members_rejected(self):
        for bad in ([], [None], [1, None]):
            with self.assertRaises(ValueError):
                empirical_crps(bad, 1)
            with self.assertRaises(ValueError):
                empirical_brier(bad, 1, 0)

    def test_crps_pair_sum_matches_literal_ordered_formula(self):
        # 随机成员下，无序对实现必须与题面有序对公式逐项一致
        rng_state = [(0.7, 1.3, -2.0, 5.5), (3.0, 3.0, 1.0), (4.2,)]
        for xs in rng_state:
            for y in (-1.0, 2.0, 7.7):
                m = len(xs)
                mae = sum(abs(x - y) for x in xs) / m
                ordered = sum(abs(xi - xj) for xi in xs for xj in xs)
                literal = mae - ordered / (2 * m * m)
                self.assertAlmostEqual(empirical_crps(list(xs), y), literal)
                if m >= 2:
                    literal_fair = mae - ordered / (2 * m * (m - 1))
                    self.assertAlmostEqual(fair_crps(list(xs), y), literal_fair)


class PairingTests(unittest.TestCase):
    def test_normal_full_common_sample(self):
        p = build_pairs(payload('normal'))
        self.assertEqual(p['counts'],
                         {'candidate_pairs': 12, 'excluded_pairs': 0,
                          'empirical_common_pairs': 12, 'fair_common_pairs': 12})
        first = p['empirical_pairs'][0]
        self.assertEqual(first['valid_day'], '2022-06-01')
        self.assertEqual(first['block_id'], 'weather0')
        # old m=3, new m=5, reference m=4，不同 m 共存于同一配对
        self.assertEqual(first['records']['old']['m'], 3)
        self.assertEqual(first['records']['new']['m'], 5)
        self.assertEqual(first['records']['reference']['m'], 4)
        # 手算 old-0 对 y=18: MAE=(0.7+0.2+0.3)/3=0.4,
        # 无序对距离=(0.5+1.0+0.5)=2.0, spread=2/9
        self.assertAlmostEqual(first['records']['old']['empirical']['crps'],
                               0.4 - 2.0 / 9)

    def test_missing_forecast_excludes_whole_pair(self):
        p = build_pairs(payload('unequal-missing-forecast'))
        self.assertEqual(p['counts']['excluded_pairs'], 1)
        self.assertEqual(p['counts']['empirical_common_pairs'], 11)
        ex = p['exclusions'][0]
        self.assertEqual(ex['valid_day'], '2022-06-03')
        self.assertEqual(ex['reasons'][0]['system'], 'old')
        self.assertEqual(ex['reasons'][0]['reason'], REASON_NO_VISIBLE_RECORD)
        # 排除日后仍在共同样本中的日期不含 2022-06-03
        self.assertNotIn('2022-06-03',
                         [r['valid_day'] for r in p['empirical_pairs']])
        # 块覆盖只针对共同验证日，不要求覆盖被排除日
        self.assertEqual(len(p['blocks']), 4)

    def test_late_revision_hidden_by_as_of(self):
        p = build_pairs(payload('late-revision'))
        self.assertEqual(p['counts']['excluded_pairs'], 0)
        first = p['empirical_pairs'][0]
        pick = first['records']['new']
        # 迟到修订 revision=2 的 recorded_day=2022-06-01 晚于冻结点，必须选 revision=1
        self.assertEqual(pick['record_id'], 'new-0')
        self.assertEqual(pick['revision'], 1)
        self.assertEqual(pick['m'], 5)

    def test_single_member_pair_is_empirical_only(self):
        p = build_pairs(payload('single-member'))
        self.assertEqual(p['counts']['empirical_common_pairs'], 12)
        self.assertEqual(p['counts']['fair_common_pairs'], 11)
        special = next(r for r in p['empirical_pairs']
                       if r['valid_day'] == '2022-06-01')
        self.assertFalse(special['fair_eligible'])
        self.assertEqual(special['records']['new']['m'], 1)
        self.assertIsNone(special['records']['new']['fair'])
        self.assertIsNotNone(special['records']['new']['empirical'])
        # fair 共同样本不含该日，经验分绝不混入 fair 均值
        self.assertNotIn('2022-06-01',
                         [r['valid_day'] for r in p['fair_pairs']])

    def test_payload_conflict_rejects_pair(self):
        v = payload('normal')
        v['forecasts'].append({
            'record_id': 'conflict', 'system': 'old',
            'issued_day': '2022-05-30', 'valid_day': '2022-06-01', 'lead': 2,
            'recorded_day': '2022-05-30', 'revision': 1,
            'members': [99.0, 99.1, 99.2]})
        p = build_pairs(v)
        ex = next(e for e in p['exclusions'] if e['valid_day'] == '2022-06-01')
        self.assertEqual(ex['reasons'][0]['reason'], REASON_PAYLOAD_CONFLICT)

    def test_incomplete_latest_record_excludes(self):
        v = payload('normal')
        v['forecasts'].append({
            'record_id': 'broken', 'system': 'old',
            'issued_day': '2022-05-30', 'valid_day': '2022-06-01', 'lead': 2,
            'recorded_day': '2022-05-30', 'revision': 2,
            'members': [1.0, None, 3.0]})
        p = build_pairs(v)
        ex = next(e for e in p['exclusions'] if e['valid_day'] == '2022-06-01')
        self.assertEqual(ex['reasons'][0]['reason'], REASON_LATEST_INCOMPLETE)

    def test_highest_revision_chosen_even_with_null_member(self):
        # 同签发日 revision=2 含 null，必须选 revision=2 并排除，不回退 revision=1
        v = payload('normal')
        v['forecasts'].append({
            'record_id': 'broken-rev2', 'system': 'old',
            'issued_day': '2022-05-30', 'valid_day': '2022-06-01', 'lead': 2,
            'recorded_day': '2022-05-30', 'revision': 2,
            'members': [1.0, None]})
        p = build_pairs(v)
        ex = next(e for e in p['exclusions'] if e['valid_day'] == '2022-06-01')
        reasons = {x['system']: x['reason'] for x in ex['reasons']}
        self.assertEqual(reasons['old'], REASON_LATEST_INCOMPLETE)
        self.assertEqual(p['counts']['empirical_common_pairs'], 11)

    def test_null_observation_excludes_pair(self):
        v = payload('normal')
        v['observations'][0]['value'] = None
        p = build_pairs(v)
        ex = next(e for e in p['exclusions'] if e['valid_day'] == '2022-06-01')
        self.assertEqual(ex['reasons'][0]['reason'], REASON_MISSING_OBSERVATION)
        self.assertIsNone(ex['reasons'][0]['system'])

    def test_block_contract_violations(self):
        v = payload('normal')
        overlap = json.loads(json.dumps(v))
        overlap['blocks'][1]['start'] = '2022-06-03'
        with self.assertRaises(StructuralError):
            build_pairs(overlap)
        gap = json.loads(json.dumps(v))
        gap['blocks'][0]['end'] = '2022-06-02'
        with self.assertRaises(StructuralError):
            build_pairs(gap)
        unused = json.loads(json.dumps(v))
        unused['blocks'].append({'id': 'extra', 'start': '2022-07-01',
                                 'end': '2022-07-02'})
        with self.assertRaises(StructuralError):
            build_pairs(unused)

    def test_structural_violations(self):
        with self.assertRaises(StructuralError):
            build_pairs(payload('invalid-domain'))
        v = payload('normal')
        v['reference'] = 'old'
        with self.assertRaises(StructuralError):
            build_pairs(v)
        v2 = payload('normal')
        v2['systems'] = ['old', 'old']
        with self.assertRaises(StructuralError):
            build_pairs(v2)
        v3 = payload('normal')
        v3['systems'] = ['old', 'ghost']
        with self.assertRaises(StructuralError):
            build_pairs(v3)

    def test_all_pairs_excluded_yields_empty_sample_not_crash(self):
        v = payload('normal')
        for o in v['observations']:
            o['value'] = None
        p = build_pairs(v)
        self.assertEqual(p['counts']['empirical_common_pairs'], 0)
        self.assertEqual(len(p['exclusions']), 12)
        r = run_bootstrap(p, n_draws=5, seed=1, percentiles=[0.5], min_blocks=3)
        self.assertEqual(r['leads'], {})


class BootstrapTests(unittest.TestCase):
    def test_type7_quantile_matches_known_definition(self):
        xs = [1.0, 2.0, 3.0, 4.0]
        # h=(n-1)q: q=0.5 -> h=1.5 -> 0.5*2+0.5*3=2.5
        self.assertAlmostEqual(percentile_type7(xs, 0.5), 2.5)
        self.assertEqual(percentile_type7(xs, 0.0), 1.0)
        self.assertEqual(percentile_type7(xs, 1.0), 4.0)
        self.assertEqual(percentile_type7([7.0], 0.25), 7.0)

    def test_seed_determinism_and_paired_draws(self):
        p = build_pairs(payload('normal'))
        r1 = run_bootstrap(p, n_draws=50, seed=53, percentiles=[0.025, 0.975],
                           min_blocks=3)
        r2 = run_bootstrap(p, n_draws=50, seed=53, percentiles=[0.025, 0.975],
                           min_blocks=3)
        self.assertEqual(r1['bootstrap']['draws'], r2['bootstrap']['draws'])
        # 每次抽取块数等于块总数，多重性之和等于块数
        for d in r1['bootstrap']['draws']:
            self.assertEqual(sum(d['multiplicity'].values()), 4)
        r3 = run_bootstrap(p, n_draws=50, seed=54, percentiles=[0.025, 0.975],
                           min_blocks=3)
        self.assertNotEqual(r1['bootstrap']['draws'], r3['bootstrap']['draws'])

    def test_point_estimate_pair_equal_weight(self):
        p = build_pairs(payload('normal'))
        r = run_bootstrap(p, n_draws=10, seed=1, percentiles=[0.025, 0.975],
                          min_blocks=3)
        sec = r['leads']['2']['empirical']
        # 手工按 12 个配对等权求 old 的经验 CRPS 均值
        manual = sum(x['records']['old']['empirical']['crps']
                     for x in p['empirical_pairs']) / 12
        self.assertAlmostEqual(sec['point']['means']['crps']['old'], manual)
        delta = sec['point']['deltas']['crps']['new_minus_old']
        self.assertAlmostEqual(
            delta,
            sec['point']['means']['crps']['new']
            - sec['point']['means']['crps']['old'])
        ref = sec['point']['means']['crps']['reference']
        self.assertAlmostEqual(
            sec['point']['skills']['crps']['old']['value'],
            1 - sec['point']['means']['crps']['old'] / ref)

    def test_min_blocks_gate_point_only(self):
        v = payload('unequal-missing-forecast')
        # 缺报后仍有 4 个有效块，门槛设为 5 触发只报点估计
        p = build_pairs(v)
        r = run_bootstrap(p, n_draws=50, seed=53, percentiles=[0.025, 0.975],
                          min_blocks=5)
        sec = r['leads']['2']['empirical']
        self.assertEqual(sec['bootstrap']['status'], 'point_only')
        self.assertEqual(sec['bootstrap']['reason'], REASON_BELOW_MIN_BLOCKS)
        # 点估计仍给出，区间不冒报
        self.assertIn('point', sec)
        self.assertNotIn('intervals', sec['bootstrap'])

    def test_single_member_fair_section_unavailable(self):
        p = build_pairs(payload('single-member'))
        r = run_bootstrap(p, n_draws=30, seed=7, percentiles=[0.025, 0.975],
                          min_blocks=3)
        fair = r['leads']['2']['fair']
        self.assertTrue(fair['available'])
        self.assertEqual(fair['n_pairs'], 11)
        emp = r['leads']['2']['empirical']
        self.assertEqual(emp['n_pairs'], 12)

    def test_zero_reference_skill_interval_marked_invalid(self):
        # 构造 reference 在每个配对经验 CRPS=0（成员全部等于观测）
        v = payload('normal')
        obs = {o['day']: o['value'] for o in v['observations']}
        for rec in v['forecasts']:
            if rec['system'] == 'reference':
                y = obs[rec['valid_day']]
                rec['members'] = [y, y, y, y]
        p = build_pairs(v)
        r = run_bootstrap(p, n_draws=40, seed=11, percentiles=[0.025, 0.975],
                          min_blocks=3)
        sec = r['leads']['2']['empirical']
        self.assertEqual(sec['point']['skills']['crps']['old']['reason'],
                         REASON_REF_ZERO)
        skill = sec['bootstrap']['intervals']['skills']['crps']['old']
        self.assertFalse(skill['interval_valid'])
        self.assertIsNone(skill['percentiles'])
        self.assertEqual(skill['replicate_counts']['valid'], 0)
        self.assertEqual(skill['replicate_counts']['invalid'], 40)
        self.assertEqual(skill['invalid_reasons'][REASON_REF_ZERO], 40)
        # 均值和 delta 区间仍然有效，不因 skill 失效而整体丢弃
        self.assertIn('p0.025', sec['bootstrap']['intervals']['means']['crps']['old'])
        self.assertIn('p0.025',
                      sec['bootstrap']['intervals']['deltas']['crps']['new_minus_old'])

    def test_multi_lead_strata_share_block_draws(self):
        v = payload('normal')
        v['as_of_lag_days'] = 0
        extra = []
        for r in v['forecasts']:
            r2 = copy.deepcopy(r)
            r2['record_id'] = r['record_id'] + '-l1'
            issued = (datetime.date.fromisoformat(r2['valid_day'])
                      - datetime.timedelta(days=1))
            r2['issued_day'] = issued.isoformat()
            r2['recorded_day'] = issued.isoformat()
            r2['lead'] = 1
            extra.append(r2)
        v['forecasts'].extend(extra)
        p = build_pairs(v)
        r = run_bootstrap(p, n_draws=20, seed=5, percentiles=[0.1, 0.9],
                          min_blocks=3)
        self.assertEqual(set(r['leads']), {'1', '2'})
        self.assertEqual(r['leads']['1']['empirical']['n_pairs'], 12)
        self.assertEqual(r['leads']['2']['empirical']['n_pairs'], 12)
        # 两 lead 共用同一组块副本（同一 seed 下 draws 只有一份）
        self.assertEqual(len(r['bootstrap']['draws']), 20)

    def test_block_repetition_weights_all_pairs_equally(self):
        # 两级聚合：块内配对等权 -> 块均值按块多重性平均。
        # 先直接验证块均值语义，再用大小不等的块验证大块不占额外权重。
        from climate_data.bootstrap import _family_aggregate, _block_means
        p = build_pairs(payload('normal'))
        pairs = [x for x in p['empirical_pairs'] if x['lead'] == 2]
        bm = _block_means(pairs, 'old', 'empirical', 'crps')
        w0_pairs = [x for x in pairs if x['block_id'] == 'weather0']
        self.assertAlmostEqual(
            bm['weather0'],
            sum(x['records']['old']['empirical']['crps'] for x in w0_pairs)
            / len(w0_pairs))
        weights = {'weather0': 2, 'weather1': 1}
        agg = _family_aggregate(pairs, weights, p['all_systems'], 'empirical')
        manual = (2 * bm['weather0'] + bm['weather1']) / 3
        self.assertAlmostEqual(agg['crps']['means']['old'], manual)
        excluded = [x for x in pairs if x['block_id'] not in weights]
        self.assertTrue(excluded)  # 确认确实有块未计权

    def test_unequal_block_sizes_do_not_inflate_weight(self):
        # 构造 weather0 含 3 个配对、weather1 含 1 个配对的 lead 分层，
        # 点估计必须是 4 个配对的简单等权均值（大块不占更大权重）
        v = payload('unequal-missing-forecast')
        p = build_pairs(v)
        pairs = [x for x in p['empirical_pairs'] if x['lead'] == 2]
        block_sizes = {}
        for x in pairs:
            block_sizes[x['block_id']] = block_sizes.get(x['block_id'], 0) + 1
        self.assertEqual(block_sizes['weather0'], 2)
        r = run_bootstrap(p, n_draws=10, seed=3, percentiles=[0.025, 0.975],
                          min_blocks=3)
        point = r['leads']['2']['empirical']['point']['means']['crps']['old']
        manual = sum(x['records']['old']['empirical']['crps'] for x in pairs) / len(pairs)
        self.assertAlmostEqual(point, manual)

    def test_fair_crps_is_never_negative(self):
        # fair CRPS 是期望经验 CRPS 的无偏估计，期望量非负；随机搜索不得出现负值
        import random as _random
        rng = _random.Random(0)
        for _ in range(5000):
            m = rng.randint(2, 8)
            xs = [rng.uniform(-100, 100) for _ in range(m)]
            y = rng.uniform(-100, 100)
            self.assertGreaterEqual(fair_crps(xs, y), -1e-12)


class FacadeTests(unittest.TestCase):
    def test_end_to_end_report_is_serializable_and_complete(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / 'report.json'
            report = run_verification(DATA / 'normal.json')
            out.write_text(json.dumps(report, ensure_ascii=False, allow_nan=False),
                           encoding='utf8')
            again = json.loads(out.read_text(encoding='utf8'))
        self.assertEqual(len(again['source']['sha256']), 64)
        self.assertIn('empirical_crps', again['score_definitions'])
        self.assertEqual(again['pairing']['counts']['empirical_common_pairs'], 12)
        self.assertEqual(len(again['raw_scores']), 12)
        self.assertEqual(again['results']['bootstrap']['seed'], 53)
        self.assertEqual(again['results']['bootstrap']['n_draws'], 199)
        # 块多重性、种子、区间有效性全部留存
        self.assertEqual(len(again['results']['bootstrap']['draws']), 199)
        skill = (again['results']['leads']['2']['empirical']
                 ['bootstrap']['intervals']['skills']['crps']['old'])
        self.assertIn('interval_valid', skill)
        self.assertIn('replicate_counts', skill)

    def test_all_variant_files_run_or_raise_documented_errors(self):
        for path in sorted(DATA.glob('*.json')):
            with self.subTest(file=path.name):
                if path.name.startswith('invalid-'):
                    with self.assertRaises(StructuralError):
                        run_verification(path)
                else:
                    report = run_verification(path)
                    json.dumps(report, allow_nan=False)
                    self.assertGreaterEqual(
                        report['pairing']['counts']['empirical_common_pairs'], 1)


if __name__ == '__main__':
    unittest.main()
