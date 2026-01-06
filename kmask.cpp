// kmask.cpp

#include "kmask.hpp"
#include "circular_queue.hpp"
#include <getopt.h>
#include <assert.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <mutex>
#include <sstream>
#include <vector>
#include <string>

std::mutex bed_mutex;       // Mutex for BED stdout writing
std::mutex stats_mutex;     // Mutex for stderr summary/logging

// -------------------- FASTA I/O --------------------

// Read one FASTA entry at a time from an existing stream
bool read_fasta(std::istream& in, std::string& header, std::string& sequence) {
    header.clear();
    sequence.clear();
    std::string line;

    // 1) Find next header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            header = line.substr(1);  // remove '>'
            break;
        }
    }
    if (header.empty()) {
        return false;
    }

    // 2) Read sequence lines until next header or EOF
    while (true) {
        std::streampos pos = in.tellg();
        if (!std::getline(in, line)) {
            // EOF
            break;
        }
        if (!line.empty() && line[0] == '>') {
            // Next header
            in.seekg(pos);  // rewind
            break;
        }

        // Allocate more space for the new line
        sequence.reserve(sequence.size() + line.size());
        for (size_t i = 0; i < line.size(); ++i) {
            unsigned char c = static_cast<unsigned char>(line[i]);
            if (c <= 32) continue;          // Skip whitespace
            if (c >= 'a' && c <= 'z') c -= 32; // ASCII uppercase
            sequence.push_back(static_cast<char>(c));
        }
    }
    return true;
}

// -------------------- ENTROPY HELPERS --------------------

// Convert base char to 0–3 code for A,C,G,T or -1 for others
static inline int base_to_code(char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default:  return -1;
    }
}

// Counts l-mers (substrings of length l) in a kmer
static inline bool count_next_lmer(const std::string& sequence,
                                   size_t start,
                                   size_t k,
                                   size_t l,
                                   std::vector<int>& counts,
                                   CircularQueue<size_t>& lmer_queue,
                                   size_t& len,
                                   size_t& bit_code) {
    // Rolling encoding of l-mers
    size_t num_states = counts.size() - 1;         // 4^l = 2^(2*l) possible l-mers
    assert(num_states == (static_cast<size_t>(1ULL) << (2 * l)));

    size_t bit_mask = num_states - 1;          // bitmask with last 2*l bits = 1

    // Encode next base into 2 bits
    int b = base_to_code(sequence[start]);
    if (b < 0) {
        len = 0;
        bit_code = num_states;
    } else {
        // Shift previous bits left by 2 and OR with the new base
        bit_code = ((bit_code << 2) | static_cast<size_t>(b)) & bit_mask;
        len = std::min(len + 1, k);
    }

    // first l-mer hasn't been encountered yet
    if (start < l - 1)
        return false;

    // we have at least one l-mer
    --counts[lmer_queue.front()];

    size_t code_to_push = len >= l ? bit_code : num_states;

    lmer_queue.push_back(code_to_push);
    ++counts[code_to_push];

    return len == k; // valid l-mers counted;
}

// Compute Shannon entropy from l-mer counts
static inline double compute_shannon_entropy(const std::vector<int>& counts, int total_lmers) {
    if (total_lmers <= 0) {
        return 0.0;
    }

    double entropy = 0.0;
    assert(counts.size());
    assert(counts.back() == 0);
    for (int c : counts) {
        if (c <= 0) continue;  // skip unused l-mer states
        double p = static_cast<double>(c) / static_cast<double>(total_lmers);
        entropy -= p * std::log2(p);
    }
    return entropy;
}

// -------------------- MASKING & COUNTING --------------------

// Mask low entropy regions: replace kmers with 'N's if entropy < threshold
std::string mask_low_entropy_regions(const std::string& sequence,
                                     size_t k, size_t l,
                                     double threshold) {
    const size_t n = sequence.size();
    if (k == 0 || l == 0 || l > k || n < k) {
        return sequence;
    }

    std::string masked = sequence;

    // Pre-allocate l-mer count array once
    // Number of possible l-mers = 4^l = 1 << (2*l)
    const size_t num_states = static_cast<size_t>(1ULL) << (2 * l);

    // add extra position for invalid lmers
    std::vector<int> counts(num_states + 1);

    const int total_lmers = static_cast<int>(k - l + 1);

    size_t bit_code = num_states;                       // rolling 2-bit-encoded l-mer
    CircularQueue<size_t> lmer_queue(total_lmers, bit_code);

    // initialize with dummy l-mers
    counts[bit_code] = total_lmers;

    // num bases encoded
    size_t len = 0;

    // Slide a k-mer window across the sequence
    for (size_t i = 0; i < n; ++i) {
        // Fill counts for all l-mers within this k-mer
        if (!count_next_lmer(sequence, i, k, l, counts, lmer_queue, len, bit_code)) {
            continue;
        }

        assert(i + 1 >= k);

        // Compute Shannon entropy (bits) for this k-mer
        double entropy = compute_shannon_entropy(counts, total_lmers);

        // Mask low-entropy k-mers with 'N'
        if (entropy < threshold) {
            masked.replace(i - k + 1, k, std::string(k, 'N'));
        }
    }

    assert(masked.size() == sequence.size());

    return masked;
}

// Count how many bases are masked (Ns)
static inline size_t count_masked_bases(const std::string& original_sequence,
                                        const std::string& masked_sequence) {
    const size_t len = std::min(original_sequence.size(), masked_sequence.size());
    size_t num_masked = 0;

    // Only account for the newly introduced Ns
    for (size_t i = 0; i < len; ++i) {
        if (original_sequence[i] != 'N' && masked_sequence[i] == 'N') {
            ++num_masked;
        }
    }
    return num_masked;
}

// -------------------- TOP-LEVEL PROCESSING --------------------

// Write masked FASTA entries to output file
void write_fasta(std::ofstream& out,
                 const std::string& header,
                 const std::string& sequence) {
    constexpr size_t line_width = 80;
    out << ">" << header << "\n";
    for (size_t i = 0; i < sequence.size(); i += line_width) {
        out << sequence.substr(i, line_width) << "\n";
    }
}

// Write masked regions in BED format to stdout, with thread-safe printing
void write_bed(const std::string& header,
               const std::string& original_sequence,
               const std::string& masked_sequence,
               const std::string& full_command) {
    static std::once_flag header_printed;

    // Print BED header once with the full command line
    std::call_once(header_printed, [&]() {
        std::lock_guard<std::mutex> lock(bed_mutex);
        std::cout << "# " << full_command << "\n";
    });

    std::vector<std::pair<size_t, size_t>> masked_regions;
    bool in_mask = false;   // true when currently inside a newly-masked stretch
    size_t start = 0;       // start coordinate of the current masked stretch

    const size_t len = std::min(original_sequence.size(), masked_sequence.size());

    for (size_t i = 0; i < len; ++i) {
        bool newly_masked = (original_sequence[i] != 'N' && masked_sequence[i] == 'N');

        if (newly_masked) {
            if (!in_mask) {
                in_mask = true;
                start = i;
            }
        } else {
            if (in_mask) {
                masked_regions.emplace_back(start, i);
                in_mask = false;
            }
        }
    }
    // When the masked stretch reaches the end of the sequence
    if (in_mask) {
        masked_regions.emplace_back(start, len);
    }

    std::lock_guard<std::mutex> lock(bed_mutex);
    for (const auto& region : masked_regions) {
        std::cout << header << "\t" << region.first << "\t" << region.second << "\n";
    }
}

// Process FASTA files and optionally outputs BED of masked regions to stdout
std::pair<std::size_t, std::size_t> process_fasta(const std::string& input_file,
                                                  size_t k, size_t l, double threshold,
                                                  const std::string& output_dir,
                                                  bool verbose,
                                                  bool output_bed,
                                                  const std::string& full_command) {
    // Derive output filename
    std::filesystem::path in_path(input_file);
    std::string stem = in_path.stem().string();   // e.g. "abc"
    std::string ext  = in_path.extension().string(); // e.g. ".fna", ".fa", ".fasta"
    std::string out_path;

    if (verbose) {
        std::ostringstream name;
        name << stem
             << "-k" << k
             << "-l" << l
             << "-s" << std::fixed << std::setprecision(2) << threshold
             << "-kmasked" << ext;   // preserve original extension
        out_path = (std::filesystem::path(output_dir) / name.str()).string();
    } else {
        out_path = (std::filesystem::path(output_dir) /
                    (stem + "-kmasked" + ext)).string();
    }

    // Open input FASTA
    std::ifstream in(input_file);
    if (!in.is_open()) {
        throw std::runtime_error("[ERROR] Cannot open file: " + input_file);
    }

    // Open output FASTA
    std::ofstream out(out_path);
    if (!out.is_open()) {
        throw std::runtime_error("[ERROR] Cannot write to file: " + out_path);
    }

    // Calculate summary stats
    size_t masked_count = 0;
    size_t total_count  = 0;

    std::string header;
    std::string sequence;

    // Stream one FASTA entry at a time
    while (read_fasta(in, header, sequence)) {
        // Mask low-entropy regions in this entry
        std::string masked_sequence = mask_low_entropy_regions(sequence, k, l, threshold);

        // Optionally write BED intervals for this entry (thread-safe inside write_bed)
        if (output_bed) {
            write_bed(header, sequence, masked_sequence, full_command);
        }

        // Write masked FASTA entry
        write_fasta(out, header, masked_sequence);

        // Accumulate stats only when verbose
        if (verbose) {
            masked_count += count_masked_bases(sequence, masked_sequence);
            total_count  += sequence.size();
        }
    }

    // Thread-safe summary + manifest
    if (verbose) {
        double percent = (total_count > 0)
            ? (100.0 * static_cast<double>(masked_count) / static_cast<double>(total_count))
            : 0.0;

        std::lock_guard<std::mutex> lock(stats_mutex);
        std::cerr << "[SUMMARY] " << input_file
                  << " | Masked " << masked_count
                  << " / " << total_count
                  << " (" << percent << "%)\n";
        std::cerr << "[MANIFEST] " << input_file
                  << " → " << out_path << "\n";
    }

    return {masked_count, total_count};
}
