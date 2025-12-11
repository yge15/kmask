#ifndef __CIRCULAR_QUEUE_HPP__
#define __CIRCULAR_QUEUE_HPP__

#include <vector>
#include <cstdint>

template <typename T>
class CircularQueue {
  public:
    CircularQueue() : front_(0), back_(-1) {}
    CircularQueue(size_t size, T a = T{}) : data_(size, std::move(a)), front_(0), back_(size - 1) {}

    size_t size() const { return data_.size(); }

    void push_back(T a) {
        assert(data_.size());
        front_ = (front_ + 1) % data_.size();
        back_ = (back_ + 1) % data_.size();
        data_[back_] = std::move(a);
    }

    void push_front(T a) {
        assert(data_.size());
        front_ = (front_ + data_.size() - 1) % data_.size();
        back_ = (back_ + data_.size() - 1) % data_.size();
        data_[front_] = std::move(a);
    }

    const T& back() const { return data_[back_]; }
    T& back() { return data_[back_]; }

    const T& front() const { return data_[front_]; }
    T& front() { return data_[front_]; }

    const std::vector<T>& data() const { return data_; }
    std::vector<T>& data() { return data_; }

  private:
    std::vector<T> data_;
    size_t front_;
    size_t back_;
};

#endif // __CIRCULAR_QUEUE_HPP__
