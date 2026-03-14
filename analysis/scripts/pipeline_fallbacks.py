#!/usr/bin/env python3
import argparse
import gzip
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


def open_maybe_gzip(path, mode="rt"):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode, encoding=None if "b" in mode else "utf-8")


def copy_text_stream(src, dst):
    with open_maybe_gzip(src, "rt") as src_handle, open_maybe_gzip(dst, "wt") as dst_handle:
        shutil.copyfileobj(src_handle, dst_handle)


def fasta_iter(path):
    header = None
    seq = []
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq)
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header is not None:
            yield header, "".join(seq)


def write_fasta(records, out_path):
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def read_fastq_sequences(path, limit=4):
    sequences = []
    with open_maybe_gzip(path, "rt") as handle:
        while len(sequences) < limit:
            header = handle.readline()
            if not header:
                break
            seq = handle.readline().strip()
            plus = handle.readline()
            qual = handle.readline()
            if not plus or not qual:
                break
            sequences.append(seq)
    return sequences


def read_first_fasta(path):
    for header, seq in fasta_iter(path):
        return header, seq
    raise ValueError(f"No FASTA records found in {path}")


def fetch_accessions(fetcher, accessions, out_fasta):
    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        acc_path = tmpdir / "acc.txt"
        acc_path.write_text("\n".join(accessions) + "\n", encoding="utf-8")
        subprocess.run(
            [sys.executable, fetcher, "--acc", str(acc_path), "--out_fasta", str(out_fasta)],
            check=True,
        )


def cmd_trim(args):
    Path(args.out_r1).parent.mkdir(parents=True, exist_ok=True)
    copy_text_stream(args.r1, args.out_r1)
    copy_text_stream(args.r2, args.out_r2)


def cmd_assemble(args):
    r1_reads = read_fastq_sequences(args.r1)
    r2_reads = read_fastq_sequences(args.r2)
    mpox_seq = "".join(r1_reads[:2]) or "ATGCCGTAACCGTTAGGCTAACCGTTAGGCTA"
    vzv_seq = "".join(r2_reads[:2]) or "TTGACCGGTTACCAAGGTTACCGGAATTCGTA"
    write_fasta(
        [
            ("mpox_contig_1", mpox_seq),
            ("vzv_contig_1", vzv_seq),
        ],
        args.out_contigs,
    )
    log_path = Path(args.log)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text("Fallback assembly used because SPAdes was not available.\n", encoding="utf-8")


def copy_or_fetch(source, out_path, fetcher):
    source_path = Path(source)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if source_path.exists():
        shutil.copyfile(source_path, out_path)
        return
    fetch_accessions(fetcher, [source], out_path)


def cmd_prepare_refs(args):
    copy_or_fetch(args.mpox_source, args.out_mpox, args.fetcher)
    copy_or_fetch(args.vzv_source, args.out_vzv, args.fetcher)


def cmd_prepare_context(args):
    records = []
    for source in args.sources:
        source_path = Path(source)
        if source_path.exists():
            records.extend(list(fasta_iter(source_path)))
        else:
            with tempfile.TemporaryDirectory() as tmpdir:
                tmp_fasta = Path(tmpdir) / "context.fasta"
                fetch_accessions(args.fetcher, [source], tmp_fasta)
                records.extend(list(fasta_iter(tmp_fasta)))
    write_fasta(records, args.out)


def cmd_mock_map(args):
    for out_path, label in ((args.sorted_bam, "sorted"), (args.dedup_bam, "deduplicated")):
        out_file = Path(out_path)
        out_file.parent.mkdir(parents=True, exist_ok=True)
        out_file.write_text(
            f"Fallback {label} BAM placeholder for {Path(args.ref).name}\n",
            encoding="utf-8",
        )


def cmd_write_depth(args):
    header, seq = read_first_fasta(args.ref)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as handle:
        for pos in range(1, len(seq) + 1):
            handle.write(f"{header.split()[0]}\t{pos}\t{args.depth_value}\n")


def cmd_write_consensus(args):
    header, seq = read_first_fasta(args.ref)
    Path(args.mask).parent.mkdir(parents=True, exist_ok=True)
    Path(args.mask).write_text("", encoding="utf-8")
    Path(args.cons).parent.mkdir(parents=True, exist_ok=True)
    write_fasta([(header, seq)], args.cons)
    Path(args.vcf).parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.vcf, "wt", encoding="utf-8") as handle:
        handle.write("##fileformat=VCFv4.2\n")
        handle.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")


def cmd_align(args):
    records = []
    for source in args.inputs:
        records.extend(list(fasta_iter(source)))
    write_fasta(records, args.out)


def cmd_tree(args):
    ids = [header.split()[0] for header, _ in fasta_iter(args.aln)]
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    if not ids:
        newick = "();\n"
    elif len(ids) == 1:
        newick = f"{ids[0]};\n"
    else:
        newick = "(" + ",".join(f"{seq_id}:1.0" for seq_id in ids) + ");\n"
    out_path.write_text(newick, encoding="utf-8")


def build_parser():
    parser = argparse.ArgumentParser(description="Fallback helpers for smoke-testing the Snakemake workflow")
    sub = parser.add_subparsers(dest="cmd", required=True)

    trim = sub.add_parser("trim")
    trim.add_argument("--r1", required=True)
    trim.add_argument("--r2", required=True)
    trim.add_argument("--out-r1", dest="out_r1", required=True)
    trim.add_argument("--out-r2", dest="out_r2", required=True)
    trim.set_defaults(func=cmd_trim)

    assemble = sub.add_parser("assemble")
    assemble.add_argument("--r1", required=True)
    assemble.add_argument("--r2", required=True)
    assemble.add_argument("--out-contigs", dest="out_contigs", required=True)
    assemble.add_argument("--log", required=True)
    assemble.set_defaults(func=cmd_assemble)

    refs = sub.add_parser("prepare-refs")
    refs.add_argument("--mpox-source", required=True)
    refs.add_argument("--vzv-source", required=True)
    refs.add_argument("--out-mpox", required=True)
    refs.add_argument("--out-vzv", required=True)
    refs.add_argument("--fetcher", required=True)
    refs.set_defaults(func=cmd_prepare_refs)

    context = sub.add_parser("prepare-context")
    context.add_argument("--sources", nargs="*", default=[])
    context.add_argument("--out", required=True)
    context.add_argument("--fetcher", required=True)
    context.set_defaults(func=cmd_prepare_context)

    mock_map = sub.add_parser("mock-map")
    mock_map.add_argument("--ref", required=True)
    mock_map.add_argument("--sorted-bam", dest="sorted_bam", required=True)
    mock_map.add_argument("--dedup-bam", dest="dedup_bam", required=True)
    mock_map.set_defaults(func=cmd_mock_map)

    depth = sub.add_parser("write-depth")
    depth.add_argument("--ref", required=True)
    depth.add_argument("--out", required=True)
    depth.add_argument("--depth-value", dest="depth_value", type=int, default=30)
    depth.set_defaults(func=cmd_write_depth)

    consensus = sub.add_parser("write-consensus")
    consensus.add_argument("--ref", required=True)
    consensus.add_argument("--vcf", required=True)
    consensus.add_argument("--mask", required=True)
    consensus.add_argument("--cons", required=True)
    consensus.set_defaults(func=cmd_write_consensus)

    align = sub.add_parser("align")
    align.add_argument("--inputs", nargs="+", required=True)
    align.add_argument("--out", required=True)
    align.set_defaults(func=cmd_align)

    tree = sub.add_parser("tree")
    tree.add_argument("--aln", required=True)
    tree.add_argument("--out", required=True)
    tree.set_defaults(func=cmd_tree)

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
