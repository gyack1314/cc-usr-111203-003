# 不同成员数集合预报的配对概率检验

这是独立的气候资料Python项目。`data/`保存合成资料，`climate_data.io`负责严格读取、公历日期和来源字节摘要，`climate_data.domain`提供现有描述统计，`climate_data.cli`将摘要与配对检验报告输出为JSON。所有观测值、站点、人群、情景与金额都是虚构。字段和研究约定见[资料约定](docs/data-contract.md)。

当前Python 3.11及以上可直接使用，没有第三方Python依赖。项目根目录运行：

```sh
python -m unittest discover -s tests -v
python -m compileall -q climate_data tests
python -m climate_data.cli inspect data/normal.json
python -m climate_data.cli verify data/normal.json --out out/verification-normal.json
```

`inspect`显示当前原始摘要，不覆盖资料。`verify`运行配对检验并落盘完整报告：评分定义与适用口径（经验/fair CRPS与Brier，均越低越好）、每项配对排除原因、经验族与fair族各自的共同样本、逐配对原始分数、气象过程块配对bootstrap的块多重性、种子、type-7分位区间与区间有效性。块由资料调用方预先冻结，不根据评分差异挑选；`invalid-`前缀文件为结构违例资料，检验会明确拒收而非静默排除。

评分口径要点：事件为观测严格大于锁定阈值；fair CRPS与fair Brier仅在该配对全部系统与reference的成员数m>=2时可估，单成员只保留经验分；fair Brier允许为负且不裁零，其修正依赖成员独立同分布假设，确定性或强相关集合需谨慎引用。

容器使用固定基础镜像，可运行`docker build -t climate-material .`，随后`docker run --rm climate-material`检查资料。`docker run --rm climate-material python -m unittest discover -s tests -v`执行已有检查。运行不依赖远端账号或实时服务。
