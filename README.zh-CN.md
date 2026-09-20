# target-gene-scrna-stagecraft

[English README](README.md) · 中文说明 · v2.1.0

[![CI](https://github.com/Guoming-Lynn/target-gene-scrna-stagecraft/actions/workflows/ci.yml/badge.svg)](https://github.com/Guoming-Lynn/target-gene-scrna-stagecraft/actions/workflows/ci.yml)

## 这是什么

这是一个围绕**一个预先指定目标基因**设计的多研究单细胞 RNA 测序分析工作流。它要回答的问题很具体：在不同研究、不同 library 的数据放在一起时，我们能不能在 donor-unit 层面，可靠地判断这个基因的信号与其他转录特征是否相关？

项目把质量控制、细胞状态整理、供体级统计和来源敏感性分析连成一条可追溯的流程。它不会把细胞条码数当成生物学重复，也不会用细胞级 DEG 代替 donor-level 分析。

更完整的安装、挂载路径和验证说明以 [英文 README](README.md) 为准。中文版覆盖同一条主路径。

## 现在能跑什么

| 形式 | 现状 |
|---|---|
| **AI skill** | `SKILL.md` 与 `references/` 中的 Part 1–6 完整规范 |
| **CLI 助手** | 审阅表、图目录、Part 5/6 构造与 verdict |
| **未捆绑** | 即开即用的 Part 1 Scanpy 跑数程序（QC、Harmony、Leiden）。请用 `scripts/part1_init.py` 加 Part 1 规范 |
| **Part 6** | 已规定、未做成 turnkey；需要用户自备、已授权的 Geneformer 环境 |

先看 [阶段路线](references/start-here.md) 和 [文档目录](references/README.md)。

**科学校准尚未完成。** smoke test 通过或 `FROZEN_PASS` 都不等于经验 FDR 控制、功效或独立复制。见 [calibration-and-provenance.md](references/calibration-and-provenance.md)。1.9.8 的 null pilot 显示 repeated-donor 模式下 BH 假发现升高，因此 `joint_common_slope` 只能作为探索性分析。

## 先试试看

需要 Python 3.10 或 3.11：

```bash
python -m pip install -r requirements.txt
python quickstart.py --out quickstart_output
```

这条命令只检查 helper 和 toy 演示路径，不会给出生物学结论，也不会执行 Geneformer。

Part 1/3 的聚类额外依赖：

```bash
python -m pip install -r requirements-part1.txt
```

## 怎么安装

```bash
python -m pip install -r requirements.txt
python -m pip install -e .
```

测试与 lint 另装 `requirements-dev.txt`。Part 5 还需要 R，以及 `Matrix`、`limma`、`edgeR`、`fgsea`、`statmod`、`jsonlite`、`yaml` 和 `digest`。Part 6 需要另建 Geneformer 环境。详细说明见 [INSTALL.md](INSTALL.md) 和 [R 依赖](references/r-requirements.md)。

解释器选择与 Windows Conda 注意项也在 INSTALL.md。维护者基线见 [validated-environment.md](references/validated-environment.md)；CI 同时覆盖 Linux 与 macOS。

本项目也可以作为 AI skill 使用。请保持目录名为 `target-gene-scrna-stagecraft`，并将包含 `SKILL.md` 的目录挂载到 Codex、Cursor 或 Claude Code 对应的 skills 目录。路径表见英文 README。

## 六个阶段

1. **Part 1：建立 atlas。** 每个 library 单独做 QC，再取严格共同基因，运行 Leiden，并由人确认分辨率和细胞类型标签。仓库提供 `part1_init.py` 和审阅表，不提供完整 QC 跑数器。
2. **Part 2：看目标基因在哪里。** 在锁定的 atlas 上描述检测率、表达强度和 donor-unit 汇总，不做 DEG。
3. **Part 3：整理一个 compartment。** 从 counts 重新计算 HVG、PCA、Harmony 和 Leiden；发现污染或低质量 cluster 时，逐轮删除并重新计算。
4. **Part 4：查看 subtype。** 在已经命名的 subtype 上重复目标基因描述，并判断哪些 donor/source 组合有足够的暴露范围。
5. **Part 5：做 donor-level 关联。** 构建排除目标基因的 pseudobulk，运行 limma/edgeR、CAMERA/fgsea、donor LOO 和 source-block LODO。
6. **Part 6：虚拟 KO/OE。** 在单独的 Geneformer 环境中比较 embedding-axis shift。它不是表达预测，也不是因果证明；当前状态是 **specified / not turnkey**。

Part 1–4 都有必须由人完成的停点。agent 可以整理候选表和图，但不能替人选择 Leiden 分辨率、KEEP/DELETE cluster 或命名 cluster。

输入对象与层： [objects-and-layers.md](references/objects-and-layers.md)。阶段目录： [stage-layout.md](references/stage-layout.md)。反模式案例： [examples.md](examples.md)。术语： [glossary.md](references/glossary.md)。

## 运行检查

与 CI 一致：

```bash
python -m pytest -q tests
python quickstart.py --out quickstart_output
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

`scripts/simulation_contract.py` 只写 manifest，不会真正跑 1000 次校准。`tests/run_r_integration.py` 是上述两个 R 脚本的包装器，不是语法检查。

## 发布 skill

```bash
python scripts/package_skill.py --out-dir release
```

发布包严格按照 `release-files.txt` 生成。发布前请核对同目录下的 `.sha256.txt` 文件。

## 引用、维护者、许可证

请使用 [CITATION.cff](CITATION.cff)。

GuomingLin · [GitHub](https://github.com/Guoming-Lynn) · [guoming.lin.med@gmail.com](mailto:guoming.lin.med@gmail.com)

[CHANGELOG](CHANGELOG.md) · [CONTRIBUTING](CONTRIBUTING.md) · [SECURITY](SECURITY.md) · [CODE OF CONDUCT](CODE_OF_CONDUCT.md)

MIT
