// kmask.hpp

#ifndef KMASK_HPP
#define KMASK_HPP

#include <string>
#include <unordered_map>
#include <vector>
#include <mutex>
#include <atomic>

// FASTA I/O
std::vector<std::pair<std::string, std::string>> read_fasta(const std::string& filename);
void write_fasta(const std::string& filename,
                 const std::vector<std::pair<std::string, std::string>>& entries);

// Entropy & Masking
std::vector<std::string> get_kmers(const std::string& sequence, size_t k);
std::unordered_map<std::string, int> count_lmers(const std::string& kmer, size_t l);
double compute_shannon_entropy(const std::unordered_map<std::string, int>& counts);
std::string mask_low_entropy_regions(const std::string& sequence, size_t k, size_t l, double threshold);

// Processing
std::pair<std::size_t, std::size_t> process_fasta(
    const std::string& input_file,
    size_t k, size_t l, double threshold,
    const std::string& output_dir,
    bool verbose, bool output_bed,
    const std::string& full_command);

// BED Output
void write_bed(const std::string& header,
               const std::string& masked_sequence,
               const std::string& full_command);

#endif // KMASK_HPP