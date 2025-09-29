#!/usr/bin/env python3
import argparse
from pathlib import Path
def fasta_iter(path):
    h = None; seq = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith('>'):
                if h is not None:
                    yield h, ''.join(seq)
                h = line[1:].split()[0]; seq = []
            else:
                seq.append(line)
        if h is not None:
            yield h, ''.join(seq)
def n50(lengths):
    tot = sum(lengths); acc = 0
    for L in sorted(lengths, reverse=True):
        acc += L
        if acc >= tot/2: return L
    return 0
ap = argparse.ArgumentParser(description='QC SPAdes contigs and report stats')
ap.add_argument('--in', dest='contigs', required=True)
ap.add_argument('--min_len', type=int, default=300)
ap.add_argument('--out_tsv', required=True)
a = ap.parse_args()
lengths = []
with open(a.out_tsv, 'w') as out:
    out.write('contig\tlength\tN_pct\tGC_pct\n')
    for h, s in fasta_iter(a.contigs):
        L = len(s)
        if L < a.min_len: continue
        n_pct = 100.0 * (s.count('N') + s.count('n')) / L
        gc = s.upper().count('G') + s.upper().count('C')
        gc_pct = 100.0 * gc / L if L else 0.0
        out.write(f'{h}\t{L}\t{n_pct:.2f}\t{gc_pct:.2f}\n')
        lengths.append(L)
if lengths:
    Path(a.out_tsv).with_suffix('.summary.txt').write_text(
        f'contigs>=min_len\t{len(lengths)}\n'
        f'sum_len\t{sum(lengths)}\n'
        f'N50\t{n50(lengths)}\n', encoding='utf-8'
    )
