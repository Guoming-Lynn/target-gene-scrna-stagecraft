# target-gene-scrna-stagecraft

[English README](README.md) · 中文说明

## 这个项目回答什么问题

对于一个预先指定的目标基因，它在多研究单细胞数据中的信号能否被诚实测量，并在 **donor-unit** 层面进行关联分析，同时保留 QC、来源敏感性和明确的结论边界？

这不是通用 Scanpy 教程：QC 按 library 独立执行，统计单位是 `dataset × donor_id`，禁止把细胞数当样本量，也禁止直接做细胞级 DEG。

## 先跑起来

```bash
python quickstart.py --out quickstart_output
python scripts/generate_toy_data.py --out toy.h5ad
```

quickstart 只运行环境和 helper smoke 检查，不产生生物学结论，也不执行 Geneformer。

一个面向**预先指定目标基因**的多研究单细胞 RNA 测序分析技能与 CLI 工作流。

## 项目做什么

项目把分析拆成六个阶段，强调可追溯、可复现和不过度声明：

1. **Part 1：** 逐 library QC、Scrublet、MAD、严格共同基因、Leiden atlas 和人工锁定标签。
2. **Part 2：** 在锁定 atlas 上描述目标基因的定位、检测率和强度，不做全基因组 DEG。
3. **Part 3：** 对选定 compartment 重新执行 HVG → PCA → Harmony → Leiden，并逐轮删除低质量 cluster。
4. **Part 4：** 在锁定 subtype 上复用目标基因目录，并评估 donor-level 可识别性。
5. **Part 5：** 以 donor-unit 为统计单位执行 target-excluded pseudobulk、limma/edgeR、CAMERA/fgsea、donor LOO 和 source-block LODO。
6. **Part 6：** 使用固定版本的 Geneformer 做虚拟 KO/OE，报告 embedding-axis shift，不把它解释为表达预测或因果证明。

科学校准尚未完成。通过 smoke test 或 `FROZEN_PASS` 不等于 FDR、统计功效或独立复制已经得到经验验证；详见
[`references/calibration-and-provenance.md`](references/calibration-and-provenance.md)。

## 安装与发布

建议使用 Python 3.10 或 3.11：

```bash
python -m pip install -r requirements.txt
python -m pip install -r requirements-dev.txt
```

Part 5 还需要 R 及 `Matrix`、`limma`、`edgeR`、`fgsea`、`statmod`、`jsonlite`、`yaml`、`digest`。
若 `Rscript` 不在 PATH，可设置 `STAGECRAFT_RSCRIPT`。

本项目也可以作为 AI skill 安装到 Codex、Cursor 或 Claude Code。请统一使用目录名 `target-gene-scrna-stagecraft`，不要改成其他名称；详细路径见 [INSTALL.md](INSTALL.md)。

## 快速验证

```bash
python quickstart.py --out quickstart_output
python scripts/_part5_smoke.py
python scripts/_part6_smoke.py
python -m unittest discover -s tests -p "test_*.py"
Rscript --vanilla tests/audit_boundaries.R
Rscript --vanilla tests/scientific_regression.R
```

正式分析前请阅读 [`SKILL.md`](SKILL.md)、[`references/start-here.md`](references/start-here.md) 和对应 Part 的规范。Part 1–4 含必须由人确认的停点：agent 只能生成候选表和图，不能自行选择 Leiden 分辨率、KEEP/DELETE 或 cluster 名称；没有人工返回的 worksheet/YAML，不得写入锁定 h5ad。图表属于分析合同的一部分，使用 `scripts/validate_figure_manifest.py` 检查统计 sidecar。

Part 1–5 是主分析路径。Part 6 是独立的 Geneformer 分支，不是开箱即用的默认安装内容。

## 发布

```bash
python scripts/package_skill.py --out-dir release
```

发布包只包含 `release-files.txt` 中明确列出的文件，不包含数据、密钥或模型权重。

## 维护者

GuomingLin · [GitHub](https://github.com/Guoming-Lynn) · [guoming.lin.med@gmail.com](mailto:guoming.lin.med@gmail.com)

## 许可证

MIT
