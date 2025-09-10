// progress.hpp

#ifndef PROGRESS_HPP
#define PROGRESS_HPP

#include <cstddef>
#include <mutex>

// Renders a simple in-place progress bar showing files processed.
void render_progress_bar(std::size_t done, std::size_t total, std::mutex& io_mutex);

#endif // PROGRESS_HPP