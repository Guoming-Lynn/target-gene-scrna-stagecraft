# target-gene-scrna-stagecraft

[English README](README.md) · 中文说明

## 这是什么

这是一个围绕**一个预先指定目标基因**设计的多研究单细胞 RNA 测序分析工作流。它要回答的问题很具体：在不同研究、不同 library 的数据放在一起时，我们能不能在 donor-unit 层面，可靠地判断这个基因的信号与其他转录特征是否相关？

项目把质量控制、细胞状态整理、供体级统计和来源敏感性分析连成一条可追溯的流程。它不会把细胞条码数当成生物学重复，也不会用细胞级 DEG 代替 donor-level 分析。

## 先试试看

准备好 Python 3.10 或 3.11 后，先运行：

```bash
python quickstart.py --out quickstart_output
python scripts/generate_toy_data.py --out toy.h5ad
```

这两条命令只用于熟悉目录结构和检查 helper 是否能运行，不会给出任何生物学结论，也不会执行 Geneformer。

## 六个阶段

1. **Part 1：建立 atlas。** 每个 library 单独做 QC，再取严格共同基因，运行 Leiden，并由人确认分辨率和细胞类型标签。
2. **Part 2：看目标基因在哪里。** 在锁定的 atlas 上描述检测率、表达强度和 donor-unit 汇总，不做 DEG。
3. **Part 3：整理一个 compartment。** 从 counts 重新计算 HVG、PCA、Harmony 和 Leiden；发现污染或低质量 cluster 时，逐轮删除并重新计算。
4. **Part 4：查看 subtype。** 在已经命名的 subtype 上重复目标基因描述，并判断哪些 donor/source 组合有足够的暴露范围。
5. **Part 5：做 donor-level 关联。** 构建排除目标基因的 pseudobulk，运行 limma/edgeR、CAMERA/fgsea、donor LOO 和 source-block LODO。
6. **Part 6：虚拟 KO/OE。** 在单独的 Geneformer 环境中比较 embedding-axis shift。它不是表达预测，也不是因果证明；当前状态是 **specified / not turnkey**。

Part 1–4 都有必须由人完成的停点。agent 可以整理候选表和图，但不能替人选择 Leiden 分辨率、KEEP/DELETE cluster 或命名 cluster。没有 `resolution_choice.yaml` 和填完整的 KEEP/DELETE CSV，就不能继续写锁定对象。

## 怎么安装

先装 Python 依赖：

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-dev.txt
```

Part 5 还需要 R，以及 `Matrix`、`limma`、`edgeR`、`fgsea`、`statmod`、`jsonlite`、`yaml` 和 `digest`。Part 6 需要另建 Geneformer 环境、固定版本的模型和用户自备的模型权重。详细说明见 [INSTALL.md](INSTALL.md)。

本项目也可以作为 AI skill 使用。请保持目录名为 `target-gene-scrna-stagecraft`，并将包含 `SKILL.md` 的目录挂载到 Codex、Cursor 或 Claude Code 对应的 skills 目录。

## 运行检查

```bash
python scripts/_part5_smoke.py
python scripts/_part6_smoke.py
python -m unittest discover -s tests -p "test_*.py"
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

通过 smoke test 只说明辅助代码和环境检查通过，不代表统计校准、功效、独立复制或 Geneformer 模型 parity 已经成立。科学校准的限制见 [`references/calibration-and-provenance.md`](references/calibration-and-provenance.md)。

## 发布 skill

```bash
python scripts/package_skill.py --out-dir release
```

发布包严格按照 `release-files.txt` 生成，不包含数据、密钥或模型权重。发布前请核对同目录下的 `.sha256.txt` 文件。

## 维护者

GuomingLin · [GitHub](https://github.com/Guoming-Lynn) · [guoming.lin.med@gmail.com](mailto:guoming.lin.med@gmail.com)

## 许可证

MIT
