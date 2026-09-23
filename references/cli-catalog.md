# CLI catalog

Examples for the helper commands. These are not a turnkey analysis and do not
change the frozen protocol. Figure specifications:

- [part1-figures.md](part1-figures.md)
- [part2-figures.md](part2-figures.md)
- [part3-figures.md](part3-figures.md)
- [part4-figures.md](part4-figures.md)
- [part5-figures.md](part5-figures.md)
- [part6-figures.md](part6-figures.md)

```bash
python scripts/part1_init.py --out analysis/01_qc_global_atlas --gene TARGET_GENE
python scripts/demo_toy_run.py --out toy_demo_output
python scripts/check_protocol.py path/to/PROTOCOL.md
python scripts/check_stage_layout.py path/to/stage_dir
python scripts/hash_inputs.py path/to/analysis_config.yaml --verify
python scripts/cluster_review_tables.py clustered.h5ad \
    --leiden-key leiden_r0_5 --out tables/Leiden_markers \
    --lineage-genes ACTA2,PECAM1,CD3D,CD79A,LYZ
python scripts/part2_figures.py locked_atlas.h5ad \
    --gene TARGET_GENE --group cell_type \
    --colors config/celltype_colors.yaml --out figures/part2
python scripts/cluster_review_tables.py round_leiden.h5ad \
    --leiden-key leiden_r0_5 --out tables/round_00_initial/markers \
    --strict-positive --exclude-genes TARGET_GENE
python scripts/part3_decision_template.py \
    tables/round_00_initial/markers/leiden_r0_5_cluster_qc.csv \
    --round round_00_initial --leiden-key leiden_r0_5 \
    --out tables/round_00_initial/manual_decision_template.csv
python scripts/part3_prepare_removal.py \
    --parent-raw objects/round_00_initial/raw_counts.h5ad \
    --parent-clustered objects/round_00_initial/leiden.h5ad \
    --leiden-key leiden_r0_5 \
    --decision tables/round_00_initial/manual_decision.csv \
    --removed-out objects/round_00_initial/removed_clusters.h5ad \
    --child-raw objects/round_01_after_removal/raw_counts.h5ad \
    --tables-out tables/round_00_initial
python scripts/part3_figures.py stop1 clustered.h5ad --out figures/round_00_initial
python scripts/part4_figures.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype \
    --colors config/subtype_colors.yaml --out figures/part4
python scripts/part5_pseudobulk.py locked_subtype.h5ad \
    --gene TARGET_GENE --group subtype --out 03_pseudobulk/
python scripts/part5_source_blocks.py 03_pseudobulk/metadata.csv \
    --map 00_protocol_manifest/source_block_map.yaml \
    --out 02_tables/
python scripts/part5_eligibility.py 02_tables/metadata_with_source_block.csv \
    --out 02_tables/eligibility.csv
Rscript scripts/part5_run_models.R 00_protocol_manifest/analysis_config.yaml
Rscript scripts/part5_run_pathways.R 00_protocol_manifest/analysis_config.yaml
python scripts/part5_verdict.py 05_logs/model_audit.json \
    --table 00_protocol_manifest/verdict_table.yaml \
    --out 05_logs/verdict.json
python scripts/part5_cell_exploratory.py cells.csv \
    --group subtype --value log1p_expr --group-a A --group-b B \
    --authorization I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK \
    --out 02_tables/exploratory/cell_level_result.json
python scripts/simulation_contract.py --seed 20260906 \
    --declared-replicates 1000 --out 00_protocol_manifest/simulation_manifest.json
python scripts/part5_figures.py holdout 02_tables/gene_effects.csv \
    --gene TARGET_GENE --out 03_figures
python scripts/part6_endpoints.py --sets endpoints.yaml \
    --model-genes model_visible_genes.txt --target TARGET_GENE \
    --out 00_input_audit/endpoint_coverage.csv
python scripts/part6_token_audit.py 02_tables/token_ledger.csv \
    --out 02_tables/token_audit.csv
python scripts/part6_sign_tests.py 02_tables/donor_effects_eligible.csv \
    --target TARGET_GENE --family-size 2 --out 02_tables/sign_tests.csv
python scripts/part6_verdict.py 05_logs/model_audit.json \
    --table 00_protocol_manifest/verdict_table.yaml \
    --out 05_logs/verdict.json
python scripts/part6_figures.py donors 02_tables/donor_effects.csv \
    --out 03_figures
python scripts/_part6_smoke.py
```

Figure QA helper:

```bash
python scripts/validate_figure_manifest.py figures/F05_08.parameters.statistics.json
```

For a formal figure, pass the complete statistics mapping to
`plotting_style.save_figure(..., statistics=...)`; it writes
`<stem>.parameters.statistics.json`, which is the input to this validator.
Rendering-only `parameters` are provenance metadata and are not a statistics
manifest.
