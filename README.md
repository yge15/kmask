# Kmask: Entropy-Based Masking of Low-Complexity Sequences

**Kmask** is a fast, multithreaded C++ tool for masking low-complexity regions in genomic sequences using Shannon entropy. It is designed as a preprocessing step for building **Kraken2 / KrakenUniq** databases, reducing false-positive classifications by replacing low-entropy k-mers with `N`s.

**Kmask** evaluates each k-mer using an entropy measure based on l-mer frequencies, selectively masking regions that fall below a user-defined entropy threshold. This improves taxonomic specificity and minimizes spurious mappings in metagenomic classification.

---

## Usage
./kmask -k \<kmer_len\> -l \<lmer_len\> -s \<threshold\> -t \<threads\> \
        [-o \<output_dir\>] [-b] [-v] \<fasta1\> [fasta2...]

## Options
- `-k`, `--kmer_len <int>`  
  Length of k-mers (default: `31`)

- `-l`, `--lmer_len <int>`  
  Length of l-mers (default: `3`)

- `-s`, `--threshold <float>`  
  Entropy threshold (default: `3.4`)

- `-t`, `--threads <int>`  
  Number of threads (default: `1`)

- `-o`, `--output-dir <path>`  
  Output directory (default: `.`)

- `-b`, `--bed`  
  Output BED of masked regions

- `-v`, `--verbose`  
  Enable detailed logging: per-file summaries and global statistics

- `-h`, `--help`  
  Display usage information

## Example
./kmask -k 31 -l 3 -s 3.4 -t 8 -o masked/ genome.fa

---

## Citation

If you use **Kmask** in published work, please cite:

> *Improving Metagenomic Classification with Kmask: Entropy-Based Masking of Low-Complexity Regions*  
> Ge et al., 2025 (preprint forthcoming)

---

## License

**Kmask** is distributed under the terms of the [GNU General Public License (GPL)](LICENSE)
