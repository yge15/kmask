# Kmask

**Entropy-based masking of low-complexity sequences**

Kmask is a fast, multithreaded C++ tool that masks low-complexity regions in genomic sequences prior to building Kraken2 or KrakenUniq databases. Low-complexity k-mers (e.g. simple repeats, homopolymer runs) are a common source of false-positive taxonomic classifications in metagenomics. Kmask identifies these regions using a Shannon entropy measure computed over l-mer frequencies and replaces low-entropy k-mers with `N`, improving classification specificity and reducing spurious hits.

## How it works

For each k-mer in the input sequence, Kmask computes an entropy score based on the frequency distribution of its constituent l-mers. K-mers whose entropy falls below a user-defined threshold are treated as low-complexity and masked (replaced with `N`). The masked FASTA files are suitable for downstream Kraken2/KrakenUniq database construction, with genuinely informative k-mers preserved and repetitive/low-information regions removed to avoid misclassifications.

## Requirements

- A C++17 compiler (`g++` recommended)
- A CPU with AVX2 support (the default build flags target `core-avx2`)
- POSIX threads (`pthread`)

## Building

```bash
git clone https://github.com/yge15/kmask.git
cd kmask
make
```

This produces a `kmask` binary in the repository root. To remove build artifacts:

```bash
make clean
```

> **Note:** The Makefile compiles with `-march=core-avx2`. If you're building on a CPU without AVX2 support, edit the `CXXFLAGS` in the `Makefile` to target your architecture (e.g. `-march=native` or a specific `-march` value).

## Usage

```bash
./kmask -k <kmer_len> -l <lmer_len> -s <threshold> -t <threads> [-o <output_dir>] [-b] [-v] <fasta1> [fasta2...]
```

### Options

| Flag | Long form | Description | Default |
|------|-----------|--------------|---------|
| `-k` | `--kmer_len <int>` | Length of k-mers | `31` |
| `-l` | `--lmer_len <int>` | Length of l-mers used for entropy calculation | `3` |
| `-s` | `--threshold <float>` | Entropy threshold below which a k-mer is masked | `3.4` |
| `-t` | `--threads <int>` | Number of worker threads | `1` |
| `-o` | `--output-dir <path>` | Directory to write masked output | `.` |
| `-b` | `--bed` | Also output a BED file of masked regions | off |
| `-v` | `--verbose` | Print per-file summaries and global statistics | off |
| `-h` | `--help` | Display usage information | — |

### Example

Mask a genome using 31-mers, 3-mers for entropy estimation, an entropy threshold of 3.4, and 8 threads:

```bash
./kmask -k 31 -l 3 -s 3.4 -t 8 -o masked/ genome.fa
```

Multiple FASTA files can be passed at once and are processed in parallel across the given thread count:

```bash
./kmask -k 31 -l 3 -s 3.4 -t 16 -o masked/ -v genome1.fa genome2.fa genome3.fa
```

## Output

Kmask writes a masked FASTA file for each input, with low-complexity k-mers replaced by `N`. These outputs are intended to be used directly as input to `kraken2-build`. When run with `-b`/`--bed`, Kmask additionally writes a BED file per input describing the coordinates of masked regions.

## Citation

If you use Kmask in published work, please cite:

> Ge et al., *Improving Metagenomic Classification with Kmask: Entropy-Based Masking of Low-Complexity Regions*, 2025 (preprint forthcoming)

## License

Kmask is distributed under the [GNU General Public License v3.0](LICENSE).