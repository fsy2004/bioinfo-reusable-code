# README 写作与持续更新规范

本规范服务于科研代码库和 Agent Skill 仓库。主页要让新读者在一分钟内回答四个问题：这是什么、能解决什么、如何跑第一个例子、哪些结论不能从当前验证中推出。

## 信息顺序

1. **名称 + 一句话价值。** 直接写用户得到什么，不写口号。
2. **Quick start。** 给最短可复制命令，并说明预期输出。
3. **能力/领域地图。** 用表格导航，不把所有文件平铺到主页。
4. **可信边界。** 区分“接口跑通”“方法已执行”“科学结论已验证”。
5. **工作方式。** 展示输入、处理、输出和人工检查点。
6. **开发与维护。** 列出真实验证命令、贡献入口、引用和许可。

段落以主题句开头；使用主动语态、具体名词和并列结构；删除重复的中英文长篇镜像。正文默认中文，保留英文名称、命令、文件名和 DOI。

## 自动生成与人工维护边界

- `<!-- readme-metrics:* -->` 与 `<!-- readme-domains:* -->` 之间由 `scripts/update_readme.py` 生成，禁止手改。
- 价值主张、Quick start、科学边界和导航由维护者编辑。
- 模块数、领域数和目录表来自实际文件树；新增或移动模块后运行生成器。
- 任何性能、测试或复现数字必须链接到当前仓库中的报告；结果失效时删除数字，不能只改措辞。

## 更新流程

```powershell
python modules\_framework\gen_catalog.py
python scripts\update_readme.py
python scripts\update_readme.py --check
python -m unittest tests.test_readme_update tests.test_2026_method_modules -v
git diff --check
```

每次领域/模块/安装方式/验证边界变化时同步更新 README。CI 在 push、pull request 和每周计划任务中运行 `--check`；发现派生区漂移就失败，提示维护者本地重生成。

## 结构参照

本规范借鉴这些仓库的信息组织，不复制其文案：

- nf-core/rnaseq: https://github.com/nf-core/rnaseq/blob/master/README.md
- Scanpy: https://github.com/scverse/scanpy/blob/main/README.md
- Snakemake: https://github.com/snakemake/snakemake/blob/main/README.md
- OpenAI Skills: https://github.com/openai/skills/blob/main/README.md
- metafor: https://github.com/wviechtb/metafor/blob/master/README.md
- revtools: https://github.com/mjwestgate/revtools/blob/master/README.md
