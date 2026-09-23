# Kmask

**Entropy-based masking of low-complexity sequences for Kraken2 / KrakenUniq databases**

Kmask is a fast, multithreaded C++ tool that finds low-complexity regions in genomic FASTA files using Shannon entropy and replaces them with `N`. It is designed as a preprocessing step before building [Kraken2](https://github.com/DerrickWood/kraken2) or [KrakenUniq](https://github.com/fbreitwieser/krakenuniq) databases: removing low-complexity k-mers reduces spurious matches and false-positive classifications in metagenomic classification.

## Installation

### Build from source

Requirements: `make`, a C++17 compiler (GCC 8 or newer, or a recent clang), and POSIX threads.

```bash
git clone https://github.com/yge15/kmask.git
cd kmask
make
```

This produces a `kmask` executable in the current directory. To install it elsewhere:

```bash
make install PREFIX=$HOME/.local   # installs to $HOME/.local/bin/kmask
```

<!-- Uncomment once the Bioconda recipe has been merged:

### Bioconda

```bash
conda install -c conda-forge -c bioconda kmask
```
-->

## Usage

```
kmask -k <kmer_len> -l <lmer_len> -s <threshold> -t <threads>
      [-o <output_dir>] [-b] [-v] <fasta1> [fasta2 ...]
```

### Options

| Option | Long form | Description | Default |
|--------|-----------|-------------|---------|
| `-k` | `--kmer_len <int>` | Window (k-mer) length | `31` |
| `-l` | `--lmer_len <int>` | l-mer length used to compute entropy | `3` |
| `-s` | `--threshold <float>` | Entropy threshold in bits; windows below it are masked | `3.3` |
| `-t` | `--threads <int>` | Number of threads | `1` |
| `-o` | `--output-dir <path>` | Output directory (created if it does not exist) | `.` |
| `-b` | `--bed` | Write a BED file of newly masked regions to standard output | off |
| `-v` | `--verbose` | Per-file summaries and global statistics | off |
| `-h` | `--help` | Print usage information | |

### Examples

Mask a genome with the default parameters:

```bash
kmask -t 8 -o masked/ genome.fa
```

Mask several files with explicit parameters and save the masked regions as BED:

```bash
kmask -k 31 -l 3 -s 3.3 -t 8 -o masked/ -b -v genomes/*.fna > masked_regions.bed
```

## Output

**Masked FASTA.** For each input file, Kmask writes a masked copy to the output directory, named `<input stem>-kmasked<original extension>` (for example, `genome.fa` becomes `genome-kmasked.fa`). With `-v`, the parameters are included in the name, as in `genome-k31-l3-s3.30-kmasked.fa`. The input files are never modified.

**BED (`-b`).** Masked regions are written to standard output, so redirect them to a file. The first line is a comment containing the full command line. Each following line has three tab-separated columns: the FASTA header, the start (0-based), and the end (exclusive). Adjacent masked bases are merged into a single interval, and bases that were already `N` in the input are not reported. Because the first column is the entire FASTA header line, headers containing spaces or tabs will not be valid BED as-is. When running with several threads, the order of records is not guaranteed.

**Logs (`-v`).** Summaries (masked bases out of total bases, per file) and global statistics are written to standard error.

### Behavior to be aware of

- Output sequences are written in uppercase with 80 bases per line, so any lowercase (soft-masked) formatting in the input is not preserved.
- Only `A`, `C`, `G`, and `T` are scored. A window containing any other character (including `N`) is not evaluated, and existing `N`s are kept as they are.
- Sequences shorter than `k` are written unchanged.
- Masking replaces bases with `N`. It does not lowercase them.

## Using Kmask with KrakenUniq

Mask your genomes directly into the `library/` directory of a custom KrakenUniq database, then build the database as usual:

```bash
DB=mydb
mkdir -p $DB/library

# 1. Mask the genomes into the database library directory
kmask -k 31 -l 3 -s 3.3 -t 8 -o $DB/library genomes/*.fna

# 2. Add the sequence ID -> taxonomy ID map
cp seqid2taxid.map $DB/

# 3. Download the taxonomy and build the database
krakenuniq-download --db $DB taxonomy
krakenuniq-build --db $DB --kmer-len 31 --threads 8
```

Notes:

- KrakenUniq reads the FASTA files in `$DB/library/`. Kmask keeps each input file's extension (`.fna`, `.fa`, or `.fasta`), so the masked files are picked up as they are.
- `seqid2taxid.map` has two tab-separated columns: the sequence ID (the FASTA header up to the first space) and its taxonomy ID. Kmask does not change FASTA headers, so the same map works for the masked files.
- Kmask's default `-k 31` matches KrakenUniq's default `--kmer-len 31`.
- `krakenuniq-download` has a `--dust` option that masks downloaded genomes with dustmasker. Leave it off for sequences you are masking with Kmask instead.

## Citation

If you use Kmask in published work, please cite:

> Ge et al. *Improving Metagenomic Classification with Kmask: Entropy-Based Masking of Low-Complexity Regions.* 2025 (preprint forthcoming).

## License

Kmask is distributed under the terms of the GNU General Public License v3.0. See [LICENSE](LICENSE) for the full text.

## Issues and contributions

Bug reports and feature requests are welcome on the [issue tracker](https://github.com/yge15/kmask/issues).