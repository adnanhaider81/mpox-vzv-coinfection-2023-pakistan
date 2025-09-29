#!/usr/bin/env python3
import argparse, time, os
from Bio import Entrez

def main():
    ap = argparse.ArgumentParser(description='Fetch GenBank FASTA by accession list')
    ap.add_argument('--email', required=False, help='Entrez email, or set env NCBI_EMAIL')
    ap.add_argument('--api_key', required=False, help='NCBI API key, or set env NCBI_API_KEY')
    ap.add_argument('--acc', required=True, help='Text file with one accession per line')
    ap.add_argument('--out_fasta', required=True)
    args = ap.parse_args()

    email = args.email or os.getenv('NCBI_EMAIL')
    if not email:
        raise SystemExit('Set --email or env NCBI_EMAIL to comply with NCBI usage policy')
    Entrez.email = email
    api_key = args.api_key or os.getenv('NCBI_API_KEY')
    if api_key:
        Entrez.api_key = api_key

    accs = [x.strip() for x in open(args.acc) if x.strip() and not x.startswith('#')]
    with open(args.out_fasta, 'w') as out:
        for acc in accs:
            h = Entrez.efetch(db='nucleotide', id=acc, rettype='fasta', retmode='text')
            out.write(h.read())
            h.close()
            time.sleep(0.34)
    print(f'Wrote {args.out_fasta}')

if __name__ == '__main__':
    main()
