// kmask.cpp

#include "kmask.hpp"
#include <getopt.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <filesystem>
#include <cmath>
#include <algorithm>
#include <mutex>

std::mutex bed_mutex;       // Mutex for BED stdout writing
std::mutex stats_mutex;     // Mutex for stderr summary/logging

// Reads multi-entry FASTA and returns vector of <header, sequence>
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("[ERROR] Cannot open file: " + filename);
    }

    std::vector<std::pair<std::string, std::string>> entries;
    std::string line, header, sequence;

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        if (line[0] == '>') {
            if (!header.empty()) {
                entries.emplace_back(header, sequence);
                sequence.clear();
            }
            header = line.substr(1);  // Remove '>'
        } else {
            for (char c : line) sequence += std::toupper(c);
        }
    }

    if (!header.empty()) {
        entries.emplace_back(header, sequence);
    }
    return entries;
}

// Counts l-mers (substrings of length l) in a kmer
std::unordered_map<std::string, int> count_lmers(const std::string& kmer, size_t l) {
    std::unordered_map<std::string, int> counts;
    if (kmer.size() < l) return counts;
    for (size_t i = 0; i <= kmer.size() - l; ++i) {
        counts[kmer.substr(i, l)]++;
    }
    return counts;
}

// Compute Shannon entropy from l-mer counts
double compute_shannon_entropy(const std::unordered_map<std::string, int>& counts) {
    double entropy = 0.0;
    int total = 0;
    for (const auto& p : counts) total += p.second;
    for (const auto& p : counts) {
        double p_i = static_cast<double>(p.second) / total;
        if (p_i > 0) entropy -= p_i * std::log2(p_i);
    }
    return entropy;
}

// Mask low entropy regions: replace kmers with 'N's if entropy < threshold
std::string mask_low_entropy_regions(const std::string& sequence, size_t k, size_t l, double threshold) {
    std::string masked = sequence;
    for (size_t i = 0; i <= sequence.size() - k; ++i) {
        std::string kmer = sequence.substr(i, k);
        
        // Skip this kmer if it contains at least one 'N'
        if (kmer.find('N') != std::string::npos) {
            continue;
        }

        auto counts = count_lmers(kmer, l);
        double entropy = compute_shannon_entropy(counts);
        if (entropy < threshold) {
            masked.replace(i, k, std::string(k, 'N'));
        }
    }
    return masked;
}

// Count how many bases are masked (N)
size_t count_masked_bases(const std::string& sequence) {
    return std::count(sequence.begin(), sequence.end(), 'N');
}

// Write masked FASTA entries to output file
void write_fasta(const std::string& filename, 
                 const std::vector<std::pair<std::string, std::string>>& entries) {
    std::ofstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("[ERROR] Cannot write to file: " + filename);
    }

    constexpr size_t line_width = 80;
    for (const auto& [header, sequence] : entries) {
        file << ">" << header << "\n";
        for (size_t i = 0; i < sequence.size(); i += line_width) {
            file << sequence.substr(i, line_width) << "\n";
        }
    }
}

// Write masked regions in BED format to stdout, with thread-safe printing
void write_bed(const std::string& header,
               const std::string& masked_sequence,
               const std::string& full_command) {
    static std::once_flag header_printed;
    std::call_once(header_printed, [&]() {
        std::lock_guard<std::mutex> lock(bed_mutex);
        std::cout << "# " << full_command << "\n";
    });

    std::vector<std::pair<size_t, size_t>> masked_regions;
    bool in_mask = false;
    size_t start = 0;

    for (size_t i = 0; i < masked_sequence.size(); ++i) {
        if (masked_sequence[i] == 'N') {
            if (!in_mask) {
                start = i;
                in_mask = true;
            }
        } else if (in_mask) {
            masked_regions.emplace_back(start, i);
            in_mask = false;
        }
    }
    if (in_mask) masked_regions.emplace_back(start, masked_sequence.size());

    std::lock_guard<std::mutex> lock(bed_mutex);
    for (const auto& [start, end] : masked_regions) {
        std::cout << header << "\t" << start << "\t" << end << "\n";
    }
}

// Process FASTA files and optionally outputs BED of masked regions to stdout
std::pair<std::size_t, std::size_t> process_fasta(
    const std::string& input_file,
    size_t k, size_t l, double threshold,
    const std::string& output_dir,
    bool verbose,
    bool output_bed,
    const std::string& full_command) 
{
    auto entries = read_fasta(input_file);
    std::vector<std::pair<std::string, std::string>> masked_entries;

    for (const auto& [header, sequence] : entries) {
        std::string masked_sequence = mask_low_entropy_regions(sequence, k, l, threshold);

        if (output_bed) {
            write_bed(header, masked_sequence, full_command);
        }

        masked_entries.emplace_back(header, masked_sequence);
    }

    std::filesystem::path in_path(input_file);
    std::string out_path;
    if (verbose) {
        std::ostringstream name;
        name << in_path.stem().string()
             << "-k" << k
             << "-l" << l
             << "-s" << std::fixed << std::setprecision(2) << threshold
             << "-kmasked.fna";
        out_path = (std::filesystem::path(output_dir) / name.str()).string();
    } else {
        out_path = (std::filesystem::path(output_dir) /
                    (in_path.stem().string() + "-kmasked.fna")).string();
    }

    write_fasta(out_path, masked_entries);

    // Calculate summary stats
    size_t masked_count = 0, total_count = 0;
    for (const auto& [_, masked] : masked_entries) {
        masked_count += count_masked_bases(masked);
        total_count += masked.size();
    }
    double percent = (total_count > 0) ? (100.0 * masked_count / total_count) : 0.0;

    // Print summary and manifest info to stderr (thread-safe)
    {
        std::lock_guard<std::mutex> lock(stats_mutex);
        std::cerr << "[SUMMARY] " << input_file << " | Masked " << masked_count
                  << " / " << total_count << " (" << percent << "%)\n";
        std::cerr << "[MANIFEST] " << input_file << " → " << out_path << "\n";
    }

    return {masked_count, total_count};
}