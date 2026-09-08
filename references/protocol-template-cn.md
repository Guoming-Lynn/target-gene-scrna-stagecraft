# 阶段协议模板（先冻结，再写代码）

复制到 `00_protocol_manifest/PROTOCOL_CN.md`。每一节都要填。空着等于没冻结。

协议版本：`x.y.z`
冻结日期：YYYY-MM-DD
证据上限：`descriptive` | `formal` | `exploratory` | `orthogonal_external` | `geometry_only`
`can_only_downgrade`：true | false
执行者：低成本模型。**本文件是唯一权威规格。** 与本文件冲突时停下报告，不得自行改设计。

Part 5 按 [part5-donor-association.md](part5-donor-association.md)。来源块 LODO 写在本协议内。`NOT_ESTIMABLE` 就是结果。

建模前冻结估计对象（`subtype_specific` 或 `joint_common_slope`）、同一供体
重复亚型的处理、source-block 的解释，以及 VIF、相关系数和条件数的预拟合审查阈值。
Part 6 按 [part6-virtual-knockout.md](part6-virtual-knockout.md)。这是 embedding shift，不是表达预测。KO 与 OE 细胞集不同时不得画成配对机制。

---

## 0. 本阶段回答什么

一句话：人群、对比或斜率、统计单位。

声明类型：观察性关联 / 快照几何顺序 / 体外药理正交 / 仅审计措辞。

正文可以写（最多 4 句）。
正文不能写（最多 6 句）。这是规格，不是文风建议。

## 1. 硬边界

- 上游目录只读，列出路径。
- 是否重聚类、改标签、重算 embedding：默认否。
- 是否下载新数据：默认否。
- 本阶段不允许升级的既有判定。

## 2. 输入（只读）

表：用途、路径、期望 shape / 必备列、冻结后 SHA-256。

每个输入必须写明它来自哪个上游锁定对象、由哪个脚本读取、缺失时停止
还是进入预先写好的降级状态。不得用当前目录里“看起来相似”的文件替代，
不得在运行后修改 hash 或把输出文件重新当成输入。记录 counts、metadata、
标签、source_block、GMT、模型字典和任何外部验证文件的完整路径。

缺列则停。不换代理基因，不换标注字段。

## 3. 细胞宇宙与标签

用哪一个冻结 annotation key。哪些标签进入、哪些先验排除（排除理由必须写在看图之前）。

不属于这条流形的谱系，要在重算 HVG/PCA **之前**去掉，然后整图重算。删点留旧图不等于重算。

## 4. 统计单位与入选

单位键。每臂/每单位最少细胞。formal / exploratory / insufficient 门槛（单位数、数据集数、残差 df、暴露变异）。

深度协变量。`duplicateCorrelation` 的 block。

## 5. 冻结清单

基因、负对照、GMT 及 SHA-256。发现库与特写面板分开。

暴露基因从结局矩阵和 leading edge 中排除；除非本阶段把它当作身份对比里的**被检验基因**（必须写明）。

## 6. 模型

主模型：设计式、软件、是否 robust eBayes。
支持模型：同一设计、另一引擎。
展示层：细胞图、UMAP、GAM 曲线。细胞上不打发现星号。

一个 FDR 家族，写名字。敏感性行不进该家族。

## 7. 停止解释的门控

身份回收、PAGA 阈值、根规则、上游系数复现锚点、共线性、模块稀疏。

失败 → 停止生物学语言；仍写覆盖和失败的门。

## 8. 稳健性

若存在多个研究：按 **source_block** 做 LODO；供体 LOO；细胞数阈值；几何层的邻域/PC。

写明一个 fold 必须保住什么（符号，不是 p 值）。

## 9. 判定表

按顺序，第一个命中生效。必须包含 `NOT_ESTIMABLE` 与 `INCONCLUSIVE`。

每一行必须写 token、触发条件、证据文件和允许的措辞后果。条件应来自
实际的单位数、rank、残差 df、暴露范围、来源留出、模型一致性和失败审计，
不能只靠“显著/不显著”。没有命中的结果必须保持 `INCONCLUSIVE`，不能由
执行者临时放宽阈值或改写统计单位来获得通过。

措辞后果：现有句子 → 必须改成的句子。

## 10. 产出

表、图（PNG/PDF/SVG + source_data + sidecar）、`verdict.json` 必备键、报告路径。

## 11. 明确不做

条列。"要是有就好了"的项目停在这里，本阶段不执行。

## 12. 执行顺序

编号。需要人眼看的检查点（例如原身份在新流形上的 UMAP）必须停，等人通过后再往下。

