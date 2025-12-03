// main.cpp

#include "kmask.hpp"
#include "progress.hpp"

#include <getopt.h>
#include <filesystem>
#include <iostream>
#include <thread>
#include <vector>
#include <mutex>
#include <atomic>
#include <iomanip>
#include <algorithm>

int main(int argc, char* argv[]) {
    // Default parameters
    std::size_t k = 31, l = 3;
    double entropy_threshold = 3.4;
    int num_threads = 1;
    bool output_bed = false;
    bool verbose = false;
    std::string output_dir = ".";

    // Store the full command line for BED header
    std::string full_command;
    for (int i = 0; i < argc; ++i) {
        if (i > 0) full_command += " ";
        full_command += argv[i];
    }

    // Parse command line arguments
    const char* short_opts = "k:l:s:t:bo:v";
    const option long_opts[] = {
        {"kmer_len",   required_argument, nullptr, 'k'},
        {"lmer_len",   required_argument, nullptr, 'l'},
        {"threshold",  required_argument, nullptr, 's'},
        {"threads",    required_argument, nullptr, 't'},
        {"bed",        no_argument,       nullptr, 'b'},
        {"output-dir", required_argument, nullptr, 'o'},
        {"verbose",    no_argument,       nullptr, 'v'},
        {nullptr, 0, nullptr, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, short_opts, long_opts, nullptr)) != -1) {
        switch (opt) {
            case 'k': k = static_cast<std::size_t>(std::stoul(optarg)); break;
            case 'l': l = static_cast<std::size_t>(std::stoul(optarg)); break;
            case 's': entropy_threshold = std::stod(optarg); break;
            case 't': num_threads = std::stoi(optarg); break;
            case 'b': output_bed = true; break;
            case 'o': output_dir = optarg; break;
            case 'v': verbose = true; break;
            default:
                std::cerr << "Usage: " << argv[0]
                          << " -k <kmer_len> -l <lmer_len> -s <threshold> -t <threads> [-o <output_dir>] [-b] [-v] <fasta1> [fasta2...]\n";
                return 1;
        }
    }

    // Gather input FASTA files
    std::vector<std::string> fasta_files;
    for (int i = optind; i < argc; ++i) {
        fasta_files.emplace_back(argv[i]);
    }

    // Validate input
    if (fasta_files.empty() || k == 0 || l == 0 || l > k) {
        std::cerr << "[ERROR] Must provide valid -k, -l, -s and at least one FASTA file (l <= k).\n";
        return 1;
    }

    // Create output directory if it doesn't exist
    if (!std::filesystem::exists(output_dir)) {
        if (verbose) {
            std::cerr << "[INFO] Creating output directory: " << output_dir << "\n";
        }
        std::error_code ec;
        if (!std::filesystem::create_directories(output_dir, ec)) {
            std::cerr << "[ERROR] Failed to create output directory: " << ec.message() << "\n";
            return 1;
        }
    }

    // Shared state for multithreading
    std::size_t current_file = 0;
    std::atomic<std::size_t> total_bases_masked{0};
    std::atomic<std::size_t> total_bases_processed{0};
    std::atomic<std::size_t> files_done{0};
    const std::size_t files_total = fasta_files.size();

    // Mutexes for synchronization
    std::mutex io_mutex;    // for std::cerr output (logs + progress bar)
    std::mutex file_mutex;  // for safe file index increment

    // Initial progress render (0/N)
    render_progress_bar(0, files_total, io_mutex);

    // Thread worker
    auto worker = [&]() {
        while (true) {
            std::size_t idx;
            {
                std::lock_guard<std::mutex> lock(file_mutex);
                if (current_file >= fasta_files.size()) return;
                idx = current_file++;
            }

            const std::string& file = fasta_files[idx];
            if (verbose) {
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cerr << "\n[INFO] Processing: " << file << "\n";
            }

            // Process FASTA file and update global counters
            auto counts = process_fasta(file, k, l, entropy_threshold, output_dir, verbose, output_bed, full_command);
            total_bases_masked += counts.first;
            total_bases_processed += counts.second;

            if (verbose) {
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cerr << "[INFO] Finished: " << file << "\n";
            }

            // Update progress bar after finishing this file
            std::size_t done_now = ++files_done;
            render_progress_bar(done_now, files_total, io_mutex);
        }
    };

    // Launch threads
    num_threads = std::max(1, num_threads);     // in case user passes -t 0 or a negative value
    std::vector<std::thread> threads;
    threads.reserve(static_cast<std::size_t>(num_threads));
    for (int i = 0; i < num_threads; ++i) {
        threads.emplace_back(worker);
    }
    for (std::thread& t : threads) {
        t.join();
    }

    // Print global summary
    if (verbose) {
        std::lock_guard<std::mutex> lock(io_mutex);
        std::cerr << "\n========== GLOBAL SUMMARY ==========\n";
        std::cerr << "Total bases processed : " << total_bases_processed.load() << "\n";
        std::cerr << "Total bases masked    : " << total_bases_masked.load() << "\n";
        if (total_bases_processed > 0) {
            double percent = 100.0 * static_cast<double>(total_bases_masked.load()) /
                             static_cast<double>(total_bases_processed.load());
            std::cerr << "Masked percent        : " << percent << "%\n";
        }
    }

    return 0;
}