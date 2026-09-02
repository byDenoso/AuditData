#!/usr/bin/env python3
from __future__ import annotations
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def main():
    result_file=ROOT/'results/r1/G12.json'
    if not result_file.is_file(): print('REFUSE: G12 verified R1 benchmark receipt is missing; historical speedups will not be plotted as fresh evidence.'); return 2
    data=json.loads(result_file.read_text(encoding='utf-8'))
    if data.get('status')!='verified': print('REFUSE: G12 is not verified'); return 2
    import matplotlib.pyplot as plt
    bench=data.get('observed',{}); batches=bench.get('batch_sizes',[]); speeds=bench.get('speedups',[])
    if not batches or len(batches)!=len(speeds): raise SystemExit('invalid G12 receipt')
    out=ROOT/'results/generated/workload_speedup.png'; out.parent.mkdir(parents=True,exist_ok=True); fig,ax=plt.subplots(); ax.plot(batches,speeds,marker='o'); ax.set_xlabel('Batch size'); ax.set_ylabel('Measured speedup'); ax.set_title('R1 verified nuisance-heavy workload'); fig.tight_layout(); fig.savefig(out,dpi=180); plt.close(fig); print(f'PASS: {out}'); return 0
if __name__=='__main__': raise SystemExit(main())
