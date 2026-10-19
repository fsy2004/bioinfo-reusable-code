# 604 · Prescription exposure audit

检索词：prescription exposure, glucocorticoid, drug dictionary, route, missingness, 药物暴露、糖皮质激素、用药核对。

用途：把研究者提供的药物和给药途径字典应用到处方记录，保留每个原始字段，输出逐记录标记和受试者层面的记录摘要。代码为原创，依赖 Python 标准库。**有用药记录不等于滥用；没有映射到的记录不等于未暴露。** 本模块不采集数据，也不提供糖皮质激素的临床判定阈值。

```powershell
python 604_prescription_exposure_audit.py --input prescriptions.csv --drug-map drug_map.csv --route-map route_map.csv --output-dir results_new
```

输入 `prescriptions.csv` 每行一条记录，必需字段为 `person_id,drug_name,route,duration_days`，其他字段原样保留；不得使用保留前缀 `audit_`。`drug_map.csv` 必须有 `alias,ingredient,drug_class`，要标记为糖皮质激素时 `drug_class` 使用 `glucocorticoid`。`route_map.csv` 必须有 `alias,route_type`，`route_type` 取 `systemic,inhaled,topical,local_other,unknown`。字典应由研究团队结合所用数据库版本核定并版本化。仅进行去首尾空格、空白折叠与大小写标准化后的精确匹配；不做模糊药名猜测或仅凭通用名推断给药途径。

输出 `medication_audit.csv`、`person_record_summary.csv` 和带输入 SHA256 的 `audit.json`。保留重复记录并标记；时长仅校验缺失、非数字、负值、非有限值。受试者摘要报告最大单记录时长，**不累加重叠处方、不等同累计疗程、不转换泼尼松等效剂量**。若需要剂量、起止日期、依从性、用药适应证或滥用定义，应建立单独经过审查的映射规则。不会覆盖已有结果。

验证：标准库行为检查覆盖原值保留、未映射药物、途径缺失、非法时长、重复记录和字典冲突；使用人工合成记录，不能证明任何数据库的暴露可测性或科学效应。来源固定记录在 `module.json` 的 `sources`。
