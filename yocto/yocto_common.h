//
// # Yocto/Common: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support collection helpers,
// Python-like iteration, timing functions, and concurrency utilities.
//
//
// ## Collection helpers
//
// The library contains a set of helpers that make it easier to use the STL.
// Here we also define the vocabulary types used in the rest of Yocto/GL.
//
// 1. check whether a value is in a container with `contain()`
//
//
// ## Python-like iterators and collection helpers
//
// This library includes a set of functions to help use C++ collections with
// more ease, inspired by Python. All functions and operators are defined in
// the yocto namespace so they will not affect the code outside. But within
// the Yocto/GL collection they are the best way to do this.
//
// 1. use `range()` to iterato over an integer sequence
// 2. use `enumerate()` to iteratare over a vector and number its elements
// 3. use opeartors + to either concatenate two vectors or a vector and an
//    element
// 4. use operators += to append an element or a vector to a given vector
//
//
// ## Concurrency utilities
//
// C++ has very basic supprt for concurrency and most of it is still platform
// dependent. We provide here very basic support for concurrency utlities
// built on top of C++ low-level threading and synchronization.
//
// 1. use `concurrent_queue()` for communicationing values between threads
// 2. use `parallel_for()` for basic parallel for loops
//
//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//

#ifndef _YOCTO_COMMON_H_
#define _YOCTO_COMMON_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <chrono>
#include <deque>
#include <future>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <utility>
#include <vector>

// -----------------------------------------------------------------------------
// DICTIONARY TYPES
// -----------------------------------------------------------------------------
namespace yocto {

// Aliased typenames for readability
using std::array;
using std::pair;
using std::string;
using std::unordered_map;
using std::vector;
using namespace std::literals::string_literals;

// For dictionaries, we use hash tables with a using directive to switch between
// implementations.
template <typename K, typename V>
using hash_map = unordered_map<K, V>;

}  // namespace yocto

// -----------------------------------------------------------------------------
// TIMING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time();

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int min, int max);
inline auto range(int max);

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template <typename T>
inline auto enumerate(const vector<T>& vals);
template <typename T>
inline auto enumerate(vector<T>& vals);

// Vector append and concatenation
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const vector<T>& b);
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const T& b);
template <typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b);
template <typename T>
inline vector<T> operator+(const vector<T>& a, const T& b);

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
  concurrent_queue() {}
  concurrent_queue(const concurrent_queue& other) = delete;
  concurrent_queue& operator=(const concurrent_queue& other) = delete;

  bool empty();
  void clear();
  void push(const T& value);
  bool try_pop(T& value);

 private:
  std::mutex    mutex;
  std::deque<T> queue;
};

// Run a task asynchronously
template <typename Func, typename... Args>
inline auto run_async(Func&& func, Args&&... args);

// Check if an async task is ready
inline bool is_valid(const std::future<void>& result);
inline bool is_running(const std::future<void>& result);
inline bool is_ready(const std::future<void>& result);

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, Func&& func);
template <typename Func>
inline void parallel_for(int num, Func&& func);

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, Func&& func);
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, Func&& func);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// TIMING UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
  return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Range object to support Python-like iteration. Use with `range()`.
struct range_helper {
  struct iterator {
    int       pos = 0;
    iterator& operator++() {
      pos++;
      return *this;
    }
    bool operator!=(const iterator& other) const { return pos != other.pos; }
    int  operator*() const { return pos; }
  };
  int      begin_ = 0, end_ = 0;
  iterator begin() const { return {begin_}; }
  iterator end() const { return {end_}; }
};

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int min, int max) { return range_helper{min, max}; }
inline auto range(int max) { return range(0, max); }

// Enumerate object to support Python-like enumeration. Use with `enumerate()`.
template <typename T>
struct enumerate_helper {
  struct iterator {
    T*        data = nullptr;
    int       pos  = 0;
    iterator& operator++() {
      pos++;
      return *this;
    }
    bool operator!=(const iterator& other) const { return pos != other.pos; }
    pair<int&, T&> operator*() const { return {pos, *(data + pos)}; }
  };
  T*       data = nullptr;
  int      size = 0;
  iterator begin() const { return {data, 0}; }
  iterator end() const { return {data, size}; }
};

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template <typename T>
inline auto enumerate(const vector<T>& vals) {
  return enumerate_helper<const T>{vals.data(), vals.size()};
}
template <typename T>
inline auto enumerate(vector<T>& vals) {
  return enumerate_helper<T>{vals.data(), vals.size()};
}

// Vector append and concatenation
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const vector<T>& b) {
  a.insert(a.end(), b.begin(), b.end());
  return a;
}
template <typename T>
inline vector<T>& operator+=(vector<T>& a, const T& b) {
  a.push_back(b);
  return a;
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const vector<T>& b) {
  auto c = a;
  return c += b;
}
template <typename T>
inline vector<T> operator+(const vector<T>& a, const T& b) {
  auto c = a;
  return c += b;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
bool concurrent_queue<T>::empty() {
  std::lock_guard<std::mutex> lock(mutex);
  return queue.empty();
}
template <typename T>
void concurrent_queue<T>::clear() {
  std::lock_guard<std::mutex> lock(mutex);
  queue.clear();
}
template <typename T>
void concurrent_queue<T>::push(const T& value) {
  std::lock_guard<std::mutex> lock(mutex);
  queue.push_back(value);
}
template <typename T>
bool concurrent_queue<T>::try_pop(T& value) {
  std::lock_guard<std::mutex> lock(mutex);
  if (queue.empty()) return false;
  value = queue.front();
  queue.pop_front();
  return true;
}

// Run a task asynchronously
template <typename Func, typename... Args>
inline auto run_async(Func&& func, Args&&... args) {
  return std::async(std::launch::async, std::forward<Func>(func),
      std::forward<Args>(args)...);
}
// Check if an async task is ready
inline bool is_valid(const std::future<void>& result) { return result.valid(); }
inline bool is_running(const std::future<void>& result) {
  return result.valid() && result.wait_for(std::chrono::microseconds(0)) !=
                               std::future_status::ready;
}
inline bool is_ready(const std::future<void>& result) {
  return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                               std::future_status::ready;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, Func&& func) {
  auto             futures  = vector<std::future<void>>{};
  auto             nthreads = std::thread::hardware_concurrency();
  std::atomic<int> next_idx(begin);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, end]() {
          while (true) {
            auto idx = next_idx.fetch_add(1);
            if (idx >= end) break;
            func(idx);
          }
        }));
  }
  for (auto& f : futures) f.get();
}

template <typename Func>
inline void parallel_for(int num, Func&& func) {
  parallel_for(0, num, std::forward<Func>(func));
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, Func&& func) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); });
}

}  // namespace yocto

#endif
