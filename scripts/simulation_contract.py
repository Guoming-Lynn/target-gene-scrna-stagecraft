#!/usr/bin/env python3
"""Create a frozen calibration manifest; execution remains stage-specific."""
from __future__ import annotations
import argparse,json
from pathlib import Path
SCENARIOS=("NULL","DONOR_EFFECT","SOURCE_CONFOUNDED","DEPTH_COLLINEAR","ONE_DONOR_DRIVEN","QC_ATTRITION","FAILED_HOLDOUT")
def main(argv=None):
    p=argparse.ArgumentParser(); p.add_argument("--out",type=Path,required=True); p.add_argument("--seed",type=int,required=True); p.add_argument("--replicates",type=int,default=1000); p.add_argument("--scenario",choices=SCENARIOS,action="append"); a=p.parse_args(argv)
    if a.replicates<100: raise SystemExit("Use at least 100 replicates")
    out={"protocol":"target-gene-scrna-stagecraft simulation calibration","seed":a.seed,"replicates":a.replicates,"scenarios":a.scenario or list(SCENARIOS),"required_outputs":["false_positive_rate","FDR","power","CI_coverage","NOT_ESTIMABLE_rate","LODO_failure_rate"],"status":"MANIFEST_ONLY","note":"Run the frozen pipeline per scenario; generic OLS is not a substitute."}
    a.out.parent.mkdir(parents=True,exist_ok=True)
    if a.out.exists(): raise SystemExit(f"Refusing to overwrite: {a.out}")
    a.out.write_text(json.dumps(out,indent=2)+"\n",encoding="utf-8")
    print(f"wrote simulation manifest: {a.out}")
if __name__ == "__main__": raise SystemExit(main())

