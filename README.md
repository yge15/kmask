# kmask: Entropy-Based Low-Complexity Masking for Kraken Databases

`kmask` is a fast, multithreaded C++ tool for masking low-complexity regions in genomic sequences using Shannon entropy. It is designed as a preprocessing step for building **Kraken / Kraken2** databases, reducing false-positive classifications by replacing low-entropy k-mers with `N`s.

`kmask` evaluates each k-mer using an entropy measure based on l-mer frequencies, selectively masking regions that fall below a user-defined entropy threshold. This improves taxonomic specificity and minimizes spurious mappings in metagenomic classification.

---

## **Features**

- **Entropy-based masking** of low-complexity DNA regions  
- **Fully multithreaded** FASTA processing  
- **BED output** of masked intervals (optional)  
- **Multi-FASTA support** (any number of records per file)  
- **Progress bar** for large database builds  
- **Verbose mode** for summary statistics per file  
- **Produces Kraken-compatible masked FASTA files**  
- Robust to filtering already-masked (`N`-containing) k-mers  

---

## **Installation**

### **Requirements**
- C++17-compatible compiler  
- GNU Make  
- g++  
- pthreads  

### **Build**

```
make
```

This produces a binary named **`kmask`**.

To clean build artifacts:

```
make clean
```

---

## **Usage**

```
kmask -k <kmer_len> -l <lmer_len> -s <entropy_cutoff> \
      -t <threads> [-b] [-v] [-o <output_dir>] <fasta1> [fasta2 ...]
```

### **Required Parameters**
| Flag | Meaning |
|------|---------|
| `-k <int>` | K-mer window size (e.g., 31 for Kraken databases) |
| `-l <int>` | Entropy unit size (e.g., l=3 for 3-mer entropy) |
| `-s <float>` | Entropy threshold for masking |
| `-t <int>` | Number of threads |

### **Optional Flags**
| Flag | Meaning |
|------|---------|
| `-b` | Output BED file of masked regions |
| `-o <dir>` | Output directory (default: current directory) |
| `-v` | Verbose mode: per-file summary + global stats |

---

## **Example**

```
kmask \
    -k 31 -l 3 -s 3.40 \
    -t 16 \
    -b -v \
    -o masked_output \
    genomes/*.fna
```

Output files will look like:

```
masked_output/
    genome1-k31-l3-s3.40-kmasked.fna
    genome2-k31-l3-s3.40-kmasked.fna
    genome1.bed   (if -b)
    genome2.bed   (if -b)
```

---

## **What kmask Produces**

### **Masked FASTA Files**
Low-entropy windows are replaced with `N`s:

```
ATGCCGATGAAAATTT...
         ^^^^^ k=31 masked
```

### **BED File (optional)**
Outputs all masked regions with the command-line call recorded:

```
# kmask -k 31 -l 3 -s 3.4 -t 8 -b input.fa
chr1    120     151
chr1    550     581
```

### **Verbose Summary**
If `-v` is used:

```
[SUMMARY] genome1.fa | Masked 18294 / 430122 (4.25%)
[MANIFEST] genome1.fa → masked_output/genome1-k31-l3-s3.40-kmasked.fna
...
========== GLOBAL SUMMARY ==========
Total bases processed : 4,300,122,000
Total bases masked    : 182,940,220
Masked percent        : 4.25%
```

---

## **Performance Notes**

- Designed for **large database builds**, including Microbial2025 (~66k genomes).  
- Thread-safe processing of many FASTA files simultaneously.  
- Skips any k-mer containing `N`s to avoid cascading masking.  
- Outperforms DUST masking for false positive suppression in Kraken classification.

---

## **Citation**

If you use **kmask** in published work, please cite:

> *Improving Metagenomic Classification with Kmask: Entropy-Based Masking of Low-Complexity Regions*  
> Ge et al., 2025 (preprint forthcoming)

---

## **License**

MIT License (or specify your preferred license)
