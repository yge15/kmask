// kmask.hpp

#ifndef KMASK_HPP
#define KMASK_HPP

#include <string>
#include <utility>
#include <cstddef>
#include <istream>
#include <fstream>

// FASTA I/O
bool read_fasta(std::istream& in, std::string& header, std::string& sequence);
void write_fasta(std::ofstream& out, const std::string& header, const std::string& sequence);

// Entropy
// double compute_shannon_entropy(const std::vector<int>& counts, int total_lmers)

// Masking
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
               const std::string& original_sequence,
               const std::string& masked_sequence,
               const std::string& full_command);

#endif // KMASK_HPP