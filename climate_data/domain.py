"""清点预报修订和成员数，不跨版本拼接。"""
from collections import Counter
from .io import date,integer,number,unique

def summarize(v):
    unique(v['forecasts'],'record_id');systems=Counter();sizes=[]
    for r in v['forecasts']:
        issued=date(r['issued_day']);verified=date(r['valid_day']);lead=integer(r['lead'])
        if (verified-issued).days!=lead:raise ValueError('提前期与日期不符')
        integer(r['revision'],1);date(r['recorded_day'])
        if not isinstance(r['system'],str) or not r['system']:raise ValueError('系统身份缺失')
        if not isinstance(r['members'],list):raise ValueError('成员应为列表')
        for x in r['members']:
            if x is not None:number(x)
        systems[r['system']]+=1;sizes.append(len(r['members']))
    days=set()
    for r in v['observations']:
        day=date(r['day'])
        if day in days:raise ValueError('验证观测日期重复')
        days.add(day)
        if r['value'] is not None:number(r['value'])
    return {'forecast_records':len(v['forecasts']),'system_records':dict(sorted(systems.items())),'member_counts':sorted(sizes),'observation_days':len(days)}
