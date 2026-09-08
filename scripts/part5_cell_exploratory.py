#!/usr/bin/env python3
"""Explicitly authorized cell-level exploratory fallback."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd
from scipy.stats import mannwhitneyu

AUTH = "I_ACCEPT_CELL_LEVEL_FALSE_POSITIVE_RISK"

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("table",type=Path); p.add_argument("--group",required=True)
    p.add_argument("--value",required=True); p.add_argument("--group-a",required=True)
    p.add_argument("--group-b",required=True); p.add_argument("--authorization",required=True)
    p.add_argument("--out",type=Path,required=True); a=p.parse_args(argv)
    if a.authorization != AUTH: raise SystemExit(f"Exact authorization required: {AUTH}")
    f=pd.read_csv(a.table); req={a.group,a.value}
    if not req.issubset(f.columns): raise SystemExit(f"missing columns: {sorted(req-set(f.columns))}")
    x=pd.to_numeric(f.loc[f[a.group].astype(str).eq(a.group_a),a.value],errors="coerce").dropna()
    y=pd.to_numeric(f.loc[f[a.group].astype(str).eq(a.group_b),a.value],errors="coerce").dropna()
    if len(x)<2 or len(y)<2: raise SystemExit("both groups need at least two finite cells")
    out={"analysis":"CELL_LEVEL_EXPLORATORY_FALLBACK","status":"EXPLORATORY_ONLY","false_positive_risk_accepted":True,"authorization":AUTH,"group":a.group,"value":a.value,"group_a":a.group_a,"group_b":a.group_b,"n_cells_a":len(x),"n_cells_b":len(y),"test":"Mann-Whitney U; cells are dependent subsamples","p_value":float(mannwhitneyu(x,y,alternative="two-sided").pvalue),"warning":"Not a donor-level finding; not replication; not a Part 5 verdict."}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if a.out.exists(): raise SystemExit(f"Refusing to overwrite: {a.out}")
    a.out.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(f"wrote exploratory cell-level result: {a.out}")
if __name__ == "__main__": raise SystemExit(main())

