//
// # Yocto/Parallel: Parallel utilities
//
// Yocto/Parallel is a collection of utilities helpful in implementing other
// Yocto/GL libraries. Yocto/Parallel is implemented in `yocto_parallel.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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

#ifndef _YOCTO_PARALLEL_H_
#define _YOCTO_PARALLEL_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <atomic>
#include <deque>
#include <future>
#include <mutex>
#include <thread>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::atomic;
using std::deque;
using std::future;
using std::vector;

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
  std::mutex mutex;
  deque<T>   queue;
};

// Run a task asynchronously
template <typename Func, typename... Args>
inline auto run_async(Func&& func, Args&&... args);

// Check if an async task is ready
inline bool is_valid(const future<void>& result);
inline bool is_running(const future<void>& result);
inline bool is_ready(const future<void>& result);

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for(T num, Func&& func);
// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the two integer indices.
template <typename T, typename Func>
inline void parallel_for(T num1, T num2, Func&& func);

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
inline bool is_valid(const future<void>& result) { return result.valid(); }
inline bool is_running(const future<void>& result) {
  return result.valid() && result.wait_for(std::chrono::microseconds(0)) !=
                               std::future_status::ready;
}
inline bool is_ready(const future<void>& result) {
  return result.valid() && result.wait_for(std::chrono::microseconds(0)) ==
                               std::future_status::ready;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename T, typename Func>
inline void parallel_for(T num, Func&& func) {
  auto      futures  = vector<future<void>>{};
  auto      nthreads = std::thread::hardware_concurrency();
  atomic<T> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, num]() {
          while (true) {
            auto idx = next_idx.fetch_add(1);
            if (idx >= num) break;
            func(idx);
          }
        }));
  }
  for (auto& f : futures) f.get();
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the two integer indices.
template <typename T, typename Func>
inline void parallel_for(T num1, T num2, Func&& func) {
  auto      futures  = vector<future<void>>{};
  auto      nthreads = std::thread::hardware_concurrency();
  atomic<T> next_idx(0);
  for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
    futures.emplace_back(
        std::async(std::launch::async, [&func, &next_idx, num1, num2]() {
          while (true) {
            auto j = next_idx.fetch_add(1);
            if (j >= num2) break;
            for (auto i = (T)0; i < num1; i++) func(i, j);
          }
        }));
  }
  for (auto& f : futures) f.get();
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
