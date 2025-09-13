// progress.cpp

#include "progress.hpp"
#include <iostream>
#include <iomanip>

// Draws and updates an in-place progress bar showing files processed
void render_progress_bar(std::size_t done, std::size_t total, std::mutex& io_mutex) {
    const std::size_t bar_width = 60;  // longer bar
    double frac = total ? static_cast<double>(done) / static_cast<double>(total) : 1.0;
    std::size_t filled = static_cast<std::size_t>(frac * static_cast<double>(bar_width));

    std::lock_guard<std::mutex> lock(io_mutex);
    std::cerr << "\r[";

    for (std::size_t i = 0; i < bar_width; ++i) {
        std::cerr << (i < filled ? "=" : " ");
    }
    
    std::cerr << "] " << std::setw(3) << static_cast<int>(frac * 100.0) << "%  "
              << done << "/" << total << " files" << std::flush;

    if (done >= total) {
        std::cerr << std::endl;
    }
}