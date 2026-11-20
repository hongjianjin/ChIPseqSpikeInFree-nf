#!/usr/bin/env python3
"""Step 1: bin the 33 Dataset 1 bigWigs of GSE108344 straight from GEO (HTTP range requests), mm9, 1 kb bins, pseudo-depth 20e6.
Run from the folder that holds dataset1_samples_and_bigwig_urls.csv, csif_bw/bin/csif_bigwig_bin.py (the tool in bin/ of the pipeline) and bench/."""
import subprocess, sys, os, time
from concurrent.futures import ThreadPoolExecutor
import pandas as pd
os.makedirs('bench/counts', exist_ok=True)
t1 = pd.read_csv('dataset1_table1_bw.csv')   # = dataset1_samples_and_bigwig_urls.csv
def work(row):
    gsm = row.gsm; outp = f'bench/counts/{gsm}'
    if os.path.exists(outp + '_bigwig_qc.json'): return gsm, 0, 'cached', 0.0
    t = time.time()
    r = subprocess.run([sys.executable, 'csif_bw/bin/csif_bigwig_bin.py', '--bigwig', row.bw, '--id', f'{gsm}.bw', '--genome', 'mm9',
                        '--prefix', outp, '--bin-size', '1000', '--pseudo-depth', '20000000'], capture_output=True, text=True)
    return gsm, r.returncode, (r.stdout.strip() or r.stderr.strip())[-300:], time.time() - t
with ThreadPoolExecutor(max_workers=3) as ex:
    for g, rc, msg, dt in ex.map(work, list(t1.itertuples())): print(g, rc, f'{dt:.0f}s', msg)
