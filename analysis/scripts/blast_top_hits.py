#!/usr/bin/env python3
import argparse, os, subprocess, csv, tempfile
from pathlib import Path

ap = argparse.ArgumentParser(description='BLAST contigs against nt and extract top hits per contig')
ap.add_argument('--contigs', required=True, help='FASTA of contigs')
ap.add_argument('--out_tsv', required=True, help='Output TSV of top hits')
ap.add_argument('--remote', action='store_true', help='Use NCBI remote BLAST')
ap.add_argument('--db', default='nt', help='Local BLAST db name when not remote')
ap.add_argument('--max_hits', type=int, default=5)
a = ap.parse_args()

cmd = ["blastn", "-query", a.contigs, "-outfmt", "6 qseqid sacc pident length evalue bitscore staxids sscinames sskingdoms"]
if a.remote:
    cmd += ["-remote", "-db", "nt"]
else:
    cmd += ["-db", a.db]
res = subprocess.run(cmd, capture_output=True, text=True, check=True)
lines = [l for l in res.stdout.strip().splitlines() if l.strip()]
rows = {}
for l in lines:
    parts = l.split("\t")
    qid = parts[0]
    rows.setdefault(qid, []).append(parts)
with open(a.out_tsv, "w") as out:
    out.write("contig	acc	pident	length	evalue	bitscore	taxid	sciname	skingdom
")
    for qid, hits in rows.items():
        for h in hits[:a.max_hits]:
            out.write("\t".join([qid] + h[1:]) + "\n")
print(f"Wrote {a.out_tsv}")
