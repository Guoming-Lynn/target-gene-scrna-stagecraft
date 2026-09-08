# Figure statistics contract

Every quantitative panel must be auditable without guessing.

Required fields in a figure statistics manifest:

```yaml
estimand: ""
biological_unit: ""
n_definition: ""
n: null
model_or_test: ""
effect_scale: ""
error_bar: "none | SD | SEM | IQR | 95% CI | other"
multiple_testing_family: "none | declared family"
adjustment: ""
claim_ceiling: ""
```

For Part 5, `n` must distinguish units, donors, datasets, source blocks and
residual df. For Part 6, donor eligibility, family size, nominal/shared-axis
status and perturbation population must be stated. A figure cannot upgrade a
`NOT_ESTIMABLE`, `SINGLE_SOURCE_DEPENDENT`, or nominal result.

Report effect sizes and intervals where available. Significance symbols are
optional and discouraged; if used, define them, show the exact test and retain
the underlying source data. Cell measurements are subsamples unless the frozen
estimand explicitly authorizes a cell-level exploratory fallback.

