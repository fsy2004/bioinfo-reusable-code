# Repository instructions

整理与管理以方便 AI 索引、检索和按需读取为目标。保持稳定的三位模块 ID、语义化英文路径和单一权威元数据；重命名同步引用与旧新路径映射。

1. 先读 `README.md`；按问题运行 `python scripts/query_modules.py <关键词>`，或检索 `catalog/modules.jsonl`。只读取命中模块及其依赖。
2. 模块用途与历史状态在 `modules/_framework/module_metadata.json`；新模块契约在各自 `module.json`。变更后执行 `python scripts/build_ai_index.py`；随后运行 `python modules/_framework/gen_catalog.py` 和 `python scripts/update_readme.py`。`--check` 检查索引漂移。
3. `catalog/path_aliases.json` 保存旧路径。模块 ID 用作引用锚点；路径中的领域/子类可以迁移。`docs/AI_INDEX.md` 解释字段、状态和检索。
4. 为本课题准备方法时读 `docs/research/2026_glucocorticoids_cancer/README.md`。论文、方法与 GitHub 来源分别登记；未实现方法保留为资源记录，不能声称已有可运行代码。
5. 原始数据、用户草稿与忽略的上游源码保留在原位；不索引其正文或把它们意外提交。第三方方法记录许可与版本，不批量复制源码。
6. 合成接口检查、依赖探测和真实数据科学验证分别报告。现有目录状态是历史记录，不能作为本机已安装/本次已验证的证据。
7. MR/共定位先运行 `python scripts/mr_coloc.py query <关键词>`，按 `catalog/mr_coloc_methods.jsonl` 的稳定方法ID按需读取。`methods.json` 是方法元数据权威；`select` 只查声明的输入条件。新增方法核原论文/官方源码/许可，保留算法、参数变体、设计/诊断、纯fine-mapping的分类；取得源码、原生接口验证和科学验证分开记，更新方法索引后再更新全库索引。
