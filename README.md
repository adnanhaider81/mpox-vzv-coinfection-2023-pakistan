# Genomic characterization of the first mpox and varicella-zoster co-infection in Pakistan (2023) through next-generation sequencing

Reproducible code and Snakemake workflow that mirror the study design and analyses reported in the Journal of Medical Virology paper on the first mpox and varicella-zoster virus (VZV) coinfection detected in Pakistan in 2023. DOI: 10.1002/jmv.29037

The default configuration is wired to bundled synthetic example inputs so the workflow can be smoke-tested from a fresh checkout. Replace those paths with your own FASTQs and preferred references for a real analysis run.

## Pipeline Overview

<p align="center">
  <img src="docs/figures/pipeline_overview.svg" alt="Publication-style overview of the mpox and VZV coinfection analysis workflow" width="100%">
</p>

[Open the full-size SVG](docs/figures/pipeline_overview.svg)

## Program summary
One end to end pipeline organized under Snakemake. It supports both metagenomic assembly driven discovery and reference based analysis, and it produces separate consensus genomes, phylogenies, and mutation reports for MPXV and VZV.

Components under one roof:
- Pathogen discovery
  - De novo SPAdes contigs from trimmed reads
  - Kraken2 classification of reads and contigs using a local database that you supply
  - Kaiju protein level classification of contigs using a local kaiju database
  - Cross check top results with NCBI BLAST either in remote mode against nt or against your local nt
  - Write a combined taxonomy summary per sample and flag likely targets such as mpox or VZV
- Targeted virus analysis
- Inputs
  - Paired end FASTQ from metagenomic DNA libraries sequenced on Illumina iSeq 2x150.
- Quality control and trimming
  - FastQC for initial QC.
  - Trimmomatic for adapter and quality trimming with sliding window 4:30 and minimum length 50, leading and trailing quality 3.
- Duplicate handling
  - Picard MarkDuplicates to flag PCR duplicates on filtered reads.
  - QC table with N50, GC percent, and percent N for contigs. Length cutoff default 300 nt.
- Closest match search and reference selection
  - BLASTN against NCBI nt using remote mode to avoid a local database, or point to a local nt database if available.
  - Parse top hits to select a best reference for MPXV and VZV. Fallback references are provided in config if BLAST is unavailable.
- Reference based mapping and consensus
  - BWA MEM mapping to the selected MPXV and VZV references, coordinate sort and index with SAMtools.
  - Depth mask at a minimum coverage threshold default 10, variant calling with bcftools mpileup and call, and masked bcftools consensus.
- Phylogeny and context
  - Fetch context accessions from GenBank using Entrez E-utilities and append study consensuses.
  - MAFFT alignment and IQ-TREE maximum likelihood tree with 1000 ultrafast bootstraps. Optional ModelFinder can be enabled.
- Mutation and clade analysis
  - Nextclade for hMPXV dataset to assign clade and list mutations for MPXV. For VZV, a simple SNP summary from VCF is reported.
- Outputs and reporting
  - Separate result trees and tables are written for MPXV and VZV.
  - Example plotting script and CI smoke test are included.

The defaults align with the manuscript methods, including iSeq 2x150, Trimmomatic parameters, duplicate handling, SPAdes, and phylogenetics with MAFFT and IQ-TREE. Where the paper used Geneious for consensus, this workflow uses bcftools with a depth mask that emulates minimum coverage greater than 10.



### Kraken2 and Kaiju databases
This repo does not ship databases. You must download and point the config at your local copies.

Authoritative sources
- Kraken2 homepage and docs: https://ccb.jhu.edu/software/kraken2/  
- Kraken2 GitHub: https://github.com/DerrickWood/kraken2  
- Prebuilt Kraken2 plus Bracken indices by the Langmead Lab (RefSeq-based): https://benlangmead.github.io/aws-indexes/k2  
- Kaiju homepage: https://kaiju.binf.ku.dk/  
- Kaiju prebuilt indexes: https://bioinformatics-centre.github.io/kaiju/downloads.html  
- NCBI Taxonomy dump (names.dmp and nodes.dmp): https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/

Quick start examples
```bash
# Kraken2: download a RefSeq-based bundle from the Langmead AWS indexes page
# Pick a tarball that matches your RAM, then untar and set kraken2.db to the extracted folder
mkdir -p /data/db/kraken2 && cd /data/db/kraken2
# example only - visit the page above to choose the exact URL
curl -O https://benlangmead.github.io/aws-indexes/k2/2025-07-xx/k2_standard_202507.tar.gz
tar -xzf k2_standard_202507.tar.gz

# Kaiju: download a prebuilt index that already includes names.dmp and nodes.dmp
mkdir -p /data/db/kaiju && cd /data/db/kaiju
# example only - visit the downloads page to choose the exact file
curl -O https://kaiju.binf.ku.dk/database/kaiju_db_nr_euk_2023-05-16.tgz
tar -xzf kaiju_db_nr_euk_2023-05-16.tgz

# Or build Kaiju DB locally using kaiju-makedb
kaiju-makedb -s refseq -t /data/db/kaiju/taxdump  # downloads taxonomy and builds the .fmi index
```

Configure paths in `config/config.yaml`:
```yaml
kraken2:
  enable: true
  db: /data/db/kraken2/k2_standard_202507   # directory that contains hash.k2d and opts.k2d

kaiju:
  enable: true
  db_fmi: /data/db/kaiju/kaiju_db_nr_euk.fmi
  nodes: /data/db/kaiju/nodes.dmp
  names: /data/db/kaiju/names.dmp
```

## Requirements
- Python 3.11 or newer
- Option A: pip and virtualenv
- Option B: conda or mamba
- Snakemake for the full pipeline
- System tools installed through conda: fastqc, trimmomatic, picard, spades, blast, bwa, samtools, bcftools, mafft, iqtree, nextclade, entrez-direct

### NCBI usage note
Set a contact email once per shell for E-utilities. Optional API key improves rate limits.
```bash
export NCBI_EMAIL="you@example.com"
export NCBI_API_KEY="xxxxxxxxxxxxxxxxxxxxxxxxxxxx"   # optional
```

## Quick verification
```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r env/requirements.txt
python analysis/scripts/example_qc_plot.py --in data-example/example_counts.tsv --out results-example/example_plot.png
```

## One command run
```bash
export NCBI_EMAIL="you@example.com"
conda env create -f env/environment.yml
conda activate mpox-vzv-env
snakemake -s workflow/Snakefile -c 4 --printshellcmds
```

## Configuration
Edit `config/config.yaml`. Minimal example:
```yaml
pairs:
  - sample: MPXV_2023_02
    r1: data-example/MPXV_2023_02_R1.fastq
    r2: data-example/MPXV_2023_02_R2.fastq

references:
  mpox: data-example/refs/MPXV.fasta
  vzv: data-example/refs/VZV.fasta

context_acc:
  mpox:
    - data-example/context/mpox_context.fasta
  vzv:
    - data-example/context/vzv_context.fasta

params:
  threads: 4
  min_len_contig: 300
  min_depth_consensus: 10
  min_qual: 20
  use_model_finder: false
  iqtree_model: GTR+G+I
  bootstrap: 1000
  blast_remote: true
```

## Outputs
- `results/consensus/mpox/<sample>.fa` and `results/consensus/vzv/<sample>.fa` - masked consensuses
- `results/aln/mpox/aln.fasta` and `results/aln/vzv/aln.fasta` - MAFFT alignments with context
- `results/iqtree/mpox.treefile` and `results/iqtree/vzv.treefile` - ML trees
- `results/blast/<sample>.contig_top_hits.tsv` - contig to best hit summary
- `results/coverage/<sample>.<virus>.depth.txt` - per base depth
- `results/mutations/mpox_nextclade.csv` - MPXV clade and mutation list if Nextclade is available
- `results/mutations/vzv_snp_summary.tsv` - per sample SNP summary relative to VZV reference
- Optional: SPAdes contig QC tables and summary under `results/spades/`


Additional discovery outputs
- `results/kraken2/<sample>.reads.report.txt` and `.contigs.report.txt` reports
- `results/kaiju/<sample>.kaiju.report.tsv` species table
- `results/taxonomy/<sample>_summary.tsv` merged summary from Kraken2 and Kaiju
- `results/taxonomy/<sample>_selected_targets.yaml` suggested targets for downstream analysis

## SPAdes contigs: QC and identification
1. Review `results/spades/<sample>/contigs.qc.tsv` and `contigs.summary.txt`. Default minimum contig length is 300 nt. N50, GC percent, and N percent are reported.
2. Identify likely MPXV and VZV contigs with BLASTN. By default the workflow runs `blastn -remote -db nt` for the top 5 hits per contig and writes `results/blast/<sample>.contig_top_hits.tsv`. Set `params.blast_remote: false` and provide a local `nt` path if you maintain a local database.
3. The workflow selects a best reference per virus using either BLAST hits or the provided fallback accessions in `config.yaml`. You can override by editing `config/selected_refs.yaml` after a first run.

## Closest references and phylogeny
1. The pipeline fetches the selected references and any additional context sequences you list in `config/config.yaml` using Entrez and combines them with your sample consensus.
2. MAFFT builds an alignment per virus. IQ-TREE then infers ML trees with 1000 ultrafast bootstraps. You can enable automatic model selection by setting `use_model_finder: true` in config.
3. The paper used models T92+G for MPXV and T92+G+I for VZV in MEGA. If you want to force these in IQ-TREE, set `iqtree_model: T92+G` or `T92+G+I` in config. ModelFinder is a good default if you are unsure.

## Mutation analysis
- MPXV: Nextclade hMPXV dataset can be used to assign clade and list substitutions and deletions. Enable by installing `nextclade` from conda and it will run automatically if present. You can also use Nextclade Web when offline and place the CSV output under `results/mutations/`.
- VZV: The workflow outputs a simple SNP table from the VCF. For amino acid effects you can extend with SnpEff if you add a VZV database.

## How to cite
- Paper: Umair M, Jamal Z, Haider SA, Hakim R, Ammar M, Ali Q, Akhtar N, Ikram A, Salman M. Genomic characterization of the first mpox and varicella-zoster co-infection in Pakistan (2023) through next-generation sequencing. Journal of Medical Virology. 2023. https://doi.org/10.1002/jmv.29037
- Software: Haider SA. Mpox and VZV coinfection Pakistan 2023 analysis. Version 3.0.0. GitHub repository.
## References
- Andrews S. 2010. FastQC. Babraham Bioinformatics.
- Bolger AM, Lohse M, Usadel B. 2014. Trimmomatic. Bioinformatics 30:2114-2120.
- Institute B. 2019. Picard Toolkit. Broad Institute.
- Bankevich A, Nurk S, Antipov D, et al. 2012. SPAdes. J Comput Biol 19:455-477.
- Camacho C, et al. 2009. BLAST+. BMC Bioinformatics 10:421.
- Li H, Durbin R. 2010. BWA. Bioinformatics 26:589-595.
- Li H, et al. 2009. SAMtools. Bioinformatics 25:2078-2079.
- Danecek P, et al. 2021. BCFtools. GigaScience 10:giab008.
- Katoh K, Standley DM. 2013. MAFFT. Mol Biol Evol 30:772-780.
- Minh BQ, et al. 2020. IQ-TREE 2. Mol Biol Evol 37:1530-1534.
- Aksamentov I, Roemer C, Hodcroft E, Neher R. 2021. Nextclade. JOSS 6:3773.
- Kans J. NCBI E-utilities Help.
- Köster J, Rahmann S. 2012. Snakemake. Bioinformatics 28:2520-2522.


- Wood DE, Lu J, Langmead B. 2019. Improved metagenomic analysis with Kraken 2. Genome Biology 20:257.
- Menzel P, Ng KL, Krogh A. 2016. Fast and sensitive taxonomic classification for metagenomics with Kaiju. Nature Communications 7:11257.
## Contributing
See `CONTRIBUTING.md`. Open an issue for questions. Do not commit restricted data.

## License
MIT. See `LICENSE` for details.
