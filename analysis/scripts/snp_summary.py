#!/usr/bin/env python3
import argparse, gzip
ap = argparse.ArgumentParser(description='Summarize SNP count from a bgzipped VCF')
ap.add_argument('--vcf', required=True)
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()
def open_auto(p):
    return gzip.open(p, 'rt') if p.endswith('.gz') else open(p, 'r')
snps = 0; ins = 0; dels = 0
with open_auto(a.vcf) as f:
    for line in f:
        if line.startswith('#'): continue
        parts = line.strip().split('\t')
        ref = parts[3]; alt = parts[4].split(',')[0]
        if len(ref)==1 and len(alt)==1:
            snps += 1
        elif len(ref) < len(alt):
            ins += 1
        else:
            dels += 1
open(a.out_tsv,'w').write(f'total_snps\t{snps}\ninsertions\t{ins}\ndeletions\t{dels}\n')
print('Wrote', a.out_tsv)
