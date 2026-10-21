# AI 按需检索入口

全库使用稳定三位模块 ID；领域与子类使用 `NN_english_slug`，模块使用 `NNN_method_purpose`。现有单文件模板保留 ID 与后缀，可由索引直接定位；不为统一外观改写科学代码。

读取顺序：根 README → `catalog/domains.json` 或关键词查询 → 命中模块 README / module.json → 实际入口。避免把整个代码库或原始数据一次装入上下文。

```powershell
python scripts/query_modules.py NHANES
python scripts/query_modules.py 癌症 --domain 17
python scripts/query_modules.py 032 --json
python scripts/query_modules.py colocalization
```

`catalog/modules.jsonl` 一行一个模块，适合 `rg` 与批量检索；`modules.json` 是同源结构化版本。字段包括稳定 ID、路径、领域/子类中英文、用途、输入输出、依赖、入口、来源、验证状态与旧路径。JSONL 的检索结果不等于已读过代码或验证运行。

权威元数据：旧模块用途来自 `modules/_framework/module_metadata.json`（2026-10-19 从当前目录冻结）与 `new_modules.json`；本次新模块优先读取本目录 `module.json`。生成文件由脚本维护，手动修改会在下次生成时被替换。

`legacy_status` 保留历史目录符号；`verification` 独立标记本次证据。`historical_metadata_not_revalidated` 表示仅定位用途，不能推断当前依赖已安装或真实数据验证成功。只有实际验证记录可以提升状态。来源未在本次核对的旧模块，sources 可为空并保留原 citation，不补造 commit/许可。

`catalog/path_aliases.json` 记录旧模块路径到新模块路径；模块 ID 不变。旧路径的后代文件可用最长前缀替换来解析；旧目录 README 提供可点击迁移导引。旧路径不是运行脚本兼容层。

仅索引 Git 可见的编号模块。不扫描被忽略的上游源码、研究数据或本地草稿；本地保留材料的内容不写入公共索引。编号配套资源可以没有主入口，标记 `reference_or_companion_resource`。

MR与共定位有第二层方法索引：`catalog/mr_coloc_methods.jsonl` / `mr_coloc_sources.jsonl`，由模块611的 `methods.json` 生成。先 `python scripts/mr_coloc.py query <关键词>`，也可 `--category mr_estimator` 或 `colocalization`。`record_kind=family_overview` 是方法族导引，不能算独立估计器；实现与参数变体也不等于独立算法。源码位于忽略的 `_upstream_sources/mr_coloc`，按来源ID和commit读取相关文件，勿整库扫描。具体条件、许可和复核流程见 [方法库说明](methods/mr_coloc/README.md)。

方法元数据变更后先运行 `python scripts/mr_coloc.py build` / `build --check`，再更新下面的全库模块索引。历史获取快照是有日期的证据；换机需重新获取并核验本地文件。

更新后执行：

```powershell
python scripts/build_ai_index.py
python modules/_framework/gen_catalog.py
python scripts/update_readme.py
python scripts/build_ai_index.py --check
python scripts/update_readme.py --check
python -m unittest tests.test_ai_index tests.test_readme_update -v
```
