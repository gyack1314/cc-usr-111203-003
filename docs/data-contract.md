# 集合预报对照

valid_day-issued_day必须等于lead。验证配对以(valid_day,lead)标识，as_of_lag_days=a表示最晚可知日为valid_day-a，issued_day与recorded_day都不得晚于该日。对每系统该配对先取最大issued_day，再取已知最大revision，记录日同样用于可见性；同业务键、同修订而载荷不同应拒绝，不从版本拼接成员。members中null是缺成员，空列表或任何null使该记录不完整。

blocks是闭区间公历日期，彼此无重叠且覆盖所有共同验证日期。同块可含多个lead，lead分层的原始均值按配对等权；每次先抽整块，再在副本内按lead汇总。fair评分要求全部比较系统和reference的该配对m>=2，经验与fair使用各自明确的共同样本，不能不披露地混用。bootstrap区间采用公开percentiles的type-7经验分位；任何无法定义的skill副本意味着该lead的skill区间不可正式给出。

## 配对检验口径（climate_data.scores / pairing / bootstrap / verify）

比较的分布评分均为越低越好的误差型评分，事件统一按观测严格大于锁定阈值（y > threshold）定义：

- 经验CRPS = mean_i|x_i-y| − 有序成员对距离和/(2m²)，等权成员经验分布，任意m>=1可估。
- fair CRPS把第二项分母改为2m(m−1)，仅m>=2可估；单成员只给经验分，禁止混入fair均值。
- 经验Brier = (p−o)²，p为成员经验事件频率、o为事件示性。
- fair Brier = (p−o)² − p(1−p)/(m−1)，仅m>=2可估，允许为负且不裁零。
- fair修正假设成员来自同一潜在预报分布的独立同分布抽样；确定性或强相关集合不自动满足，报告同时保留经验与fair两套口径。

配对身份是(valid_day, lead)。每系统只取issued_day与recorded_day均不晚于valid_day−as_of_lag_days的记录，先最大issued_day再最大revision；同签发日同revision载荷不同则整对拒绝，成员不跨版本拼接。任一比较系统、reference或观测不完整时整对排除，报告给出每对的排除原因代码与中文说明，以及经验族与fair族各自的最终共同样本；禁止按系统各自算完再相减。

点估计在lead分层内逐配对等权，与块大小无关。bootstrap以调用方预先冻结的整块为重抽单位：每副本有放回抽取与块总数等量的块，块重复则块内全部配对随块重复，所有系统、reference、lead与评分族共用同一份块副本；副本内先在块内对配对取等权均值，再按块多重性对块均值平均，块内日期数不同不会获得额外权重。delta与skill=1−系统均分/reference均分均在副本内重算；reference均分为零（或负）的副本保留失效诊断、不参与该skill分位，并标注interval_valid=false，绝不丢弃后冒报完整区间。有效块数低于min_blocks时该lead只报点估计。原始分数、每副本块多重性、种子与区间有效性随报告落盘。块只能由调用方给定，不得依据评分差异挑选。
