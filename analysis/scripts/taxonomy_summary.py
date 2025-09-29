#!/usr/bin/env python3
import argparse, csv, re, sys, pathlib

ap = argparse.ArgumentParser(description='Merge Kraken2 and Kaiju reports to summarize taxa and flag likely viruses')
ap.add_argument('--k2_reads', required=False, help='Kraken2 report on reads')
ap.add_argument('--k2_contigs', required=False, help='Kraken2 report on contigs')
ap.add_argument('--kaiju', required=False, help='Kaiju table report (species)')
ap.add_argument('--out_tsv', required=True)
ap.add_argument('--out_targets', required=True, help='YAML listing selected viruses if detected')
args = ap.parse_args()

def parse_kraken_report(path):
    rows = []
    if not path: return rows
    with open(path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 5: 
                parts = re.split(r'\s{2,}', line.strip())
                if len(parts) < 5: 
                    continue
            try:
                pct = float(parts[0])
            except:
                continue
            taxid = parts[-2] if parts[-2].isdigit() else ''
            name = parts[-1].strip() if parts[-1] else parts[3].strip()
            rank = parts[3].strip()[0] if len(parts[3].strip())>0 else ''
            rows.append({'source':'kraken2','name':name, 'pct':pct, 'rank':rank, 'taxid':taxid})
    return rows

def parse_kaiju_table(path):
    rows = []
    if not path: return rows
    with open(path) as f:
        reader = csv.reader(f, delimiter='\t')
        header = next(reader, None)
        for r in reader:
            if not r or len(r) < 7: 
                continue
            name = r[1]
            try:
                reads = int(r[2])
            except:
                continue
            pct = float(r[3]) if r[3] not in ('', 'NA') else 0.0
            rows.append({'source':'kaiju','name':name, 'pct':pct, 'rank':'S', 'taxid':r[6] if len(r)>6 else ''})
    return rows

rows = []
rows += parse_kraken_report(args.k2_reads)
rows += parse_kraken_report(args.k2_contigs)
rows += parse_kaiju_table(args.kaiju)

# Aggregate species level
from collections import defaultdict
agg = defaultdict(lambda: {'pct':0.0, 'sources':set()})
for r in rows:
    key = r['name']
    agg[key]['pct'] += r['pct']
    agg[key]['sources'].add(r['source'])

ordered = sorted(agg.items(), key=lambda x: x[1]['pct'], reverse=True)

with open(args.out_tsv, 'w') as out:
    out.write('name\tcombined_pct\tsources\n')
    for name, val in ordered:
        out.write(f"{name}\t{val['pct']:.3f}\t{','.join(sorted(val['sources']))}\n")

# Detect targets
targets = []
names_norm = [n.lower() for n,_ in ordered]
def present_any(substrs):
    s = '|'.join(substrs)
    return any(re.search(s, n) for n in names_norm)

if present_any(['monkeypox', 'mpox', 'orthopoxvirus']):
    targets.append('mpox')
if present_any(['varicella', 'vzv', 'varicella-zoster']):
    targets.append('vzv')

import yaml
with open(args.out_targets, 'w') as y:
    yaml.safe_dump({'selected_targets': targets}, y)
print('Wrote', args.out_tsv, 'and', args.out_targets)
