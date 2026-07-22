# Nanopore Bacterial Genome Assembly & Comparative Genomics Pipeline

Working notes and scripts for a Nanopore long-read bacterial genome assembly
workflow, covering everything from raw read transfer to genome annotation,
phylogenetics, and pangenome/comparative genomics analysis. Compiled during
my master's thesis research.

**Status: work in progress.** Some sections (marked 🚧 below) are still
being completed as the analysis develops.

All standalone scripts referenced below live in [`scripts/`](scripts/),
organized by language (`bash/`, `python/`, `R/`).

## Table of contents

1. [Transferring read data with gdown](#1-transferring-read-data-with-gdown)
2. [Genome assembly with Flye](#2-genome-assembly-with-flye)
3. [Assembly polishing with Medaka](#3-assembly-polishing-with-medaka)
4. [Quality assessment with QUAST](#4-quality-assessment-with-quast)
5. [Extracting individual contigs](#5-extracting-individual-contigs)
6. [Evaluating contig sizes](#6-evaluating-contig-sizes)
7. [Reorienting circular genomes with dnaapler](#7-reorienting-circular-genomes-with-dnaapler)
8. [Genome annotation with Bakta](#8-genome-annotation-with-bakta)
9. [16S rRNA extraction with barrnap](#9-16s-rrna-extraction-with-barrnap)
10. [Multiple sequence alignment with MAFFT](#10-multiple-sequence-alignment-with-mafft)
11. [Phylogenetic tree inference with IQ-TREE](#11-phylogenetic-tree-inference-with-iq-tree)
12. [Pangenome analysis with Roary](#12-pangenome-analysis-with-roary)
13. [Functional annotation with eggNOG-mapper](#13-functional-annotation-with-eggnog-mapper)
14. [KEGG pathway visualization with KEGGaNOG](#14-kegg-pathway-visualization-with-kegganog)
15. [🚧 OrthoVenn](#15-orthovenn)
16. [PPanGGOLiN](#16-ppanggolin)
17. [Genomic island synteny with clinker](#17-genomic-island-synteny-with-clinker)
18. [Plant growth-promoting gene detection with PGPT-finder](#18-plant-growth-promoting-gene-detection-with-pgpt-finder)

---

## 1. Transferring read data with gdown

[`gdown`](https://github.com/wkentaro/gdown) downloads public files/folders
from Google Drive. Compared to `curl`/`wget`, it:

- Bypasses the "can't scan large file for viruses" notice that blocks
  `curl`/`wget` downloads.
- Recursively downloads every file in a folder (up to 50 files per folder).
- Can export Google Docs/Sheets/Slides directly to PDF/XML/CSV.

### Installation

```bash
pip install gdown
```

If you see:

```
WARNING: The script gdown is installed in '/home/aguerra/.local/bin' which is not on PATH.
```

it means `gdown` was installed correctly, but its location isn't on your
shell's `PATH`, so the `gdown` command won't be found. Add it permanently by
appending this line to `~/.bashrc`:

```bash
export PATH="$HOME/.local/bin:$PATH"
```

Then reload the shell config:

```bash
source ~/.bashrc
```

If you're working over SSH, closing and reopening the session also works.
Verify with:

```bash
which gdown
# should print: /home/aguerra/.local/bin/gdown
```

### Downloading a folder by its Google Drive ID

```bash
gdown LINK_DRIVE_FOLDER -O OUTPUT_DIR_NAME --folder
```

This creates `OUTPUT_DIR_NAME` containing every file from that Drive folder.

## 2. Genome assembly with Flye

Once the reads are verified and merged into a single FASTQ file,
[Flye](https://github.com/mikolmogorov/Flye) assembles the genome.

### Installation

```bash
conda create -n flye
conda activate flye
conda install flye
```

### Usage

```bash
flye --nano-hq all_reads.fastq --out-dir flye_assembly --genome-size 5.0m
```

## 3. Assembly polishing with Medaka

After assembly, [Medaka](https://github.com/nanoporetech/medaka) polishes the
draft assembly, producing a `consensus.fasta` with sequencing errors
corrected.

### Installation

```bash
conda create -n medaka -c conda-forge -c nanoporetech -c bioconda medaka
conda activate medaka
```

### Usage

For native bacterial data (isolates, metagenomic samples, or bacterially
expressed plasmids), a bacteria-specific model improves consensus accuracy
and is compatible with several R10 chemistry basecaller versions. The
`--bacteria` flag selects it automatically when compatible with the input
basecaller:

```bash
medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} --bacteria
```

`${BASECALLS}` is the merged `all_reads.fastq` file from step 1. Medaka
writes several files to `${OUTDIR}`; the one used downstream is
`consensus.fasta`.

## 4. Quality assessment with QUAST

[QUAST](https://github.com/ablab/quast) reports assembly quality metrics,
including:

- Number and total length of large contigs (> 500 bp).
- Length of the largest contig.
- N50 (the contig length such that all contigs at least that long together
  cover ≥ 50% of the assembly).
- Number of predicted genes (via GeneMark.hmm for prokaryotes,
  GeneMark-ES/GlimmerHMM for eukaryotes, or MetaGeneMark for metagenomes).

### Installation

```bash
conda create -n quast
conda activate quast
conda install -c bioconda quast
```

### Usage

```bash
quast.py consensus.fasta -o quast_report
```

This generates a `quast_report/` folder with the quality reports.

## 5. Extracting individual contigs

When assembly produces multiple contigs (e.g., chromosome plus plasmids),
individual contigs can be extracted from the multi-FASTA with `awk`. It
starts printing when it reaches the target header and stops at the next
`>` header:

```bash
awk '/^>contig_1$/{p=1} /^>/ && !/^>contig_1$/{p=0} p' consensus.fasta > contig_1.fasta
```

To extract every contig in the multi-FASTA at once:

```bash
awk '/^>/ {out = substr($1, 2) ".fasta"} {print >> out}' consensus.fasta
```

This writes one FASTA file per contig, named after its header, to the same
directory as `consensus.fasta`.

## 6. Evaluating contig sizes

```bash
for f in *.fasta; do awk '/^>/ {next} {l += length($0)} END {print FILENAME, l}' $f; done
```

Example output:

```
consensus_chromo.fasta 4459216
contig_1.fasta 4194028
contig_2.fasta 146552
contig_3.fasta 107201
contig_5.fasta 4878
contig_7.fasta 4673
contig_9.fasta 1884
```

## 7. Reorienting circular genomes with dnaapler

[dnaapler](https://github.com/gbouras13/dnaapler) reorients complete
circular genomes to a consistent start position.

### Installation

```bash
conda create -n dnaapler_env
conda activate dnaapler_env
conda install -c bioconda dnaapler
```

### Usage

```bash
# from a mixed-contig FASTA
dnaapler all -i input_mixed_contigs.fasta -o output_directory_path -p my_bacteria_name -t 8

# from a GFA file (e.g., from Flye, Unicycler, or Autocycler)
dnaapler all -i assembly.gfa -o output_directory_path -p my_bacteria_name -t 8
```

Key outputs:

- `dnaapler_<timestamp>.log`
- `logs/`
- `<prefix>_all_reorientation_summary.tsv`
- `<prefix>_MMseqs2_output.txt`
- `<prefix>_reoriented.fasta` — the reoriented assembly used downstream.

## 8. Genome annotation with Bakta

### Installation

```bash
conda create -n bakta
conda activate bakta
conda install -c conda-forge -c bioconda bakta
```

### Database setup

[Bakta](https://github.com/oschwengers/bakta) requires a database hosted on
Zenodo, available in `full` and `light` variants. `full` gives the best
annotation results; `light` trades some completeness for a smaller
footprint and faster download.

List available versions:

```bash
bakta_db list
```

Download and set up the latest compatible version:

```bash
bakta_db download --output <output-path> --type [light|full]
```

Or install manually:

```bash
wget https://zenodo.org/record/14916843/files/db-light.tar.xz
bakta_db install -i db-light.tar.xz
```

Update the AMRFinderPlus DB independently if needed:

```bash
amrfinder_update --force_update --database db-light/amrfinderplus-db/
```

### Usage — single genome

```bash
bakta --db path/to/dbdir --prefix PREFIX_NAME --output outputdir/name --verbose "$file"
```

### Usage — batch over multiple genomes

See [`scripts/bash/bakta_batch_annotation.sh`](scripts/bash/bakta_batch_annotation.sh).

## 9. 16S rRNA extraction with barrnap

[barrnap](https://github.com/tseemann/barrnap) predicts ribosomal RNA genes.

### Installation

```bash
conda create -n barrnap
conda activate barrnap
conda install -c bioconda -c conda-forge barrnap
```

### Usage

```bash
barrnap --kingdom bac --threads 8 genome.fasta --outseq rrna_sequences.fasta > rrna_locations.gff
```

The `--kingdom` model determines which genes get extracted (e.g., 16S/23S
for bacteria via `bac`, 18S for eukaryotes via `euk`).

`rrna_sequences.fasta` typically contains several 16S/23S hits, since these
genes are usually present in multiple copies across the genome:

![16S/23S rRNA sequences extracted by barrnap](images/16s_23s_rrna_sequences.png)

`rrna_locations.gff` reports the genomic location of each hit (in this
example, all hits fall on contig 1):

![rRNA gene locations in the GFF output](images/rrna_locations_gff.png)

### Renaming sequences to genus + species

NCBI FASTA headers follow the order:
`ACCESSION GENUS EPITHET <additional information>`. This command keeps only
genus and epithet for renaming — always spot-check the result manually:

```bash
sed -E '/^>/ s/^>[^ ]+ ([^ ]+) ([^ ]+).*/>\1_\2/' 16s_sequences_raw.fasta > 16s_sequences_final.fasta
```

## 10. Multiple sequence alignment with MAFFT

### Installation

```bash
sudo apt install mafft
```

Key flags:

- `--op #` — gap opening penalty (default: 1.53)
- `--ep #` — offset, acts like a gap extension penalty (default: 0.0)
- `--maxiterate #` — max iterative refinement passes (default: 0)
- `--clustalout` — output in Clustal format (default: FASTA)
- `--reorder` — output sequences in aligned order (default: input order)
- `--quiet` — suppress progress output
- `--thread #` — number of threads (`-1` lets MAFFT decide)
- `--dash` — add structural information (Rozewicki et al., submitted)

If unsure which settings to use, let MAFFT choose automatically:

```bash
mafft --auto in.fasta > out.fasta
```

Standard alignment used in this pipeline:

```bash
mafft --reorder --adjustdirection --auto input.fasta > output.fasta
```

- `--reorder` — output in aligned order
- `--adjustdirection` — reorients sequences relative to the first sequence
- `--auto` — automatic strategy selection

## 11. Phylogenetic tree inference with IQ-TREE

### Installation

```bash
conda create -n iqtree
conda activate iqtree
conda install -c bioconda iqtree
```

### Usage — general maximum-likelihood tree

Combines ModelFinder, tree search, SH-aLRT, and 1000 ultrafast bootstrap
replicates:

```bash
iqtree -s example.phy -B 1000 -alrt 1000
```

## 12. Pangenome analysis with Roary

[Roary](https://github.com/sanger-pathogens/Roary) builds the pangenome
from a set of annotated genomes.

### Usage

Run on GFF files produced by Bakta/PROKKA annotation:

```bash
roary -f ./roary_output -e -r -n -v -p 7 *.gff
```

### Plots

- **Rarefaction curve** (conserved vs. total genes across genomes):
  [`scripts/R/pangenome_rarefaction_curve.R`](scripts/R/pangenome_rarefaction_curve.R)
- **UpSet plot** of genes shared with a genome of interest:
  [`scripts/R/pangenome_upset_plot.R`](scripts/R/pangenome_upset_plot.R)
- **Frequency histogram, tree + presence/absence matrix, and composition
  donut chart** — a modified version of `roary_plots.py`:
  [`scripts/python/roary_plots_custom.py`](scripts/python/roary_plots_custom.py)

  ```bash
  python scripts/python/roary_plots_custom.py <tree.newick> <gene_presence_absence.csv>
  ```

  The Newick tree can be generated with IQ-TREE from Roary's
  `core_gene_alignment.aln`.

### 🚧 Extracting singletons of a genome of interest

Not yet implemented.

## 13. Functional annotation with eggNOG-mapper

### COG extraction with eggnogCOGextractor

[eggnogCOGextractor](https://github.com/vsmicrogenomics/eggnogCOGextractor)
summarizes COG categories from eggNOG-mapper output.

**Installation:** download `eggnogCOGextractor.py` — no further install
needed.

**Usage:** create an `input/` and an `output/` folder. Place the
`*.emapper.annotations.tsv` file in `input/`, but rename it to drop the
`.tsv` extension so it ends in `.emapper.annotations`. Then run:

```bash
python eggnogCOGextractor.py
```

This writes a summary TSV to `output/`:

![eggnogCOGextractor summary output](images/eggnogcogextractor_output.png)

This file feeds the COG bar plots below.

### COG frequency bar plot

[`scripts/python/cog_barplot_overview.py`](scripts/python/cog_barplot_overview.py)
plots gene counts per COG category from the `.emapper.cog.tsv` output:

![COG category bar chart](images/cog_bar_chart_final.png)

Adjust the color palette and other visual details as needed for your data.

### Detailed bar plot for specific COG categories

[`scripts/python/cog_barplot_category_detail.py`](scripts/python/cog_barplot_category_detail.py)
plots the top annotated gene functions within selected COG categories (e.g.,
`Q` — secondary metabolites, `V` — defense mechanisms) from the raw
`.emapper.annotations` file:

![Detailed function comparison for the Q and V COG categories](images/cog_detailed_plot_example.png)

## 14. KEGG pathway visualization with KEGGaNOG

[KEGGaNOG](https://github.com/iliapopov17/KEGGaNOG) visualizes KEGG pathway
completeness from eggNOG-mapper output.

### Installation

```bash
conda create -n kegganog pip -y
conda activate kegganog
pip install kegganog
```

### Usage — single-mode completeness plot

```bash
KEGGaNOG -i MM_rqbj1453.emapper.annotations.tsv -o output_dir
```

`-i` takes the eggNOG-mapper output file, `-o` the output directory
(created automatically if missing).

### 🚧 Completeness bar plot

[`scripts/python/kegg_completeness_barplot.py`](scripts/python/kegg_completeness_barplot.py)
plots completeness from the `*_pathways.tsv` file generated by KEGGaNOG.
**Still needs SVG export and a configurable color palette** — see the TODO
in the script.

## 15. 🚧 OrthoVenn

Not yet documented — analysis pending.

## 16. PPanGGOLiN

[PPanGGOLiN](https://github.com/labgem/PPanGGOLiN) identifies genomic
islands (regions of genomic plasticity, or "spots") that may be shared
across multiple genomes.

🚧 Usage details pending — analysis in progress.

## 17. Genomic island synteny with clinker

The GenBank files produced by Bakta are used as input to
[clinker](https://github.com/gamcil/clinker), restricted to the coordinate
ranges of each genomic island. Those ranges come from PPanGGOLiN's
`spot_X_identical_rgps.tsv` files, which list the coordinates of each
region of genomic plasticity (RGP).

Once clinker generates its HTML output, open it in a browser to see the
full genomes with coordinates, then zoom to the region of interest (e.g.,
around position 1,179,962 in a given genome) to inspect a specific island.

## 18. Plant growth-promoting gene detection with PGPT-finder

[PGPg_finder](https://github.com/tpellegrinetti/PGPg_finder) screens
genomes for plant growth-promoting genes.

### Installation

```bash
git clone https://github.com/tpellegrinetti/PGPg_finder/
bash PGPg_finder/install.sh
```

This creates a dedicated `PGPg_finder` conda environment with all
dependencies.

### Usage — single genome

```bash
python PGPg_finder/PGPg_finder.py -w genome_wf -i GENOME_FILE.fna -o pgp_chromo -t 20
```

### Heatmap and category summary

Get `heatmap.py` from the PGPg_finder repository, then run:

```bash
python heatmap.py gene_counts.txt heatmap_output/ pathways_plabase.txt summary.txt
```

`pathways_plabase.txt` and `summary.txt` are fetched via `wget` from the
link in the PGPg_finder repository, or found in its install directory.