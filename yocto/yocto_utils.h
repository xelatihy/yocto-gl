//
// # Yocto/Utils: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils is a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support simple path manipulation, 
// and concurrency utilities. These utilities are likely to change often and 
// are to be considered internal to Yocto.
//
//
// ## Path manipulation
//
// We define a few path manipulation utilities to split and join path
// components.
//
// 1. Get paths components with `get_dirname()`, `get_filename()` and
//   `get_extension()`
// 2. Remove parts of a path with get_noextension() and `get_basename()`
// 3. check if a file exists with `exists_file()`
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
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

#ifndef _YOCTO_UTILS_H_
#define _YOCTO_UTILS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

#include <ctype.h>
#include <stdio.h>

#include <chrono>
#include <deque>
#include <future>
#include <mutex>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::atomic;
using std::deque;
using std::future;
using std::runtime_error;
using std::string;
using std::string_view;
using std::vector;
using namespace std::string_literals;
using namespace std::string_view_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// PRINTING AND FORMATTING HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

// print information and returns a timer that will print the time when
// destroyed. Use with RIIA for scoped timing.
struct timer {
  timer() : start{get_time()} {}
  int64_t elapsed() { return get_time() - start; }
  string  elapsedf() {
    auto duration = get_time() - start;
    auto elapsed  = duration / 1000000;  // milliseconds
    auto hours    = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs  = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buffer[256];
    sprintf(buffer, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    return buffer;
  }
  static int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
  }

 private:
  int64_t start = 0;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

inline string normalize_path(const string& filename_) {
  auto filename = filename_;
  for (auto& c : filename)
    if (c == '\\') c = '/';
  if (filename.size() > 1 && filename[0] == '/' && filename[1] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  if (filename.size() > 3 && filename[1] == ':' && filename[2] == '/' &&
      filename[3] == '/') {
    throw std::invalid_argument("absolute paths are not supported");
    return filename_;
  }
  auto pos = (size_t)0;
  while ((pos = filename.find("//")) != filename.npos)
    filename = filename.substr(0, pos) + filename.substr(pos + 1);
  return filename;
}

// Get directory name (including '/').
inline string get_dirname(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return "";
  return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
inline string get_extension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return "";
  return filename.substr(pos + 1);
}

// Get filename without directory.
inline string get_filename(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('/');
  if (pos == string::npos) return filename;
  return filename.substr(pos + 1);
}

// Get extension.
inline string get_noextension(const string& filename_) {
  auto filename = normalize_path(filename_);
  auto pos      = filename.rfind('.');
  if (pos == string::npos) return filename;
  return filename.substr(0, pos);
}

// Get filename without directory and extension.
inline string get_basename(const string& filename) {
  return get_noextension(get_filename(filename));
}

// Return the preset type and the remaining filename
inline bool is_preset_filename(const string& filename) {
  return filename.find("::yocto::") == 0;
}
// Return the preset type and the filename. Call only if this is a preset.
inline pair<string, string> get_preset_type(const string& filename) {
  if (filename.find("::yocto::") == 0) {
    auto aux = filename.substr(string("::yocto::").size());
    auto pos = aux.find("::");
    if (pos == aux.npos) throw std::runtime_error("bad preset name" + filename);
    return {aux.substr(0, pos), aux.substr(pos + 2)};
  } else {
    return {"", filename};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
  // a simple concurrent queue that locks at every call
  concurrent_queue() {}
  concurrent_queue(const concurrent_queue<T>& other) {
    if (!other._queue.empty())
      throw std::invalid_argument("cannot copy full queue");
    clear();
  }
  concurrent_queue<T>& operator=(const concurrent_queue<T>& other) {
    if (!other._queue.empty())
      throw std::invalid_argument("cannot copy full queue");
    clear();
  }

  bool empty() {
    std::lock_guard<std::mutex> lock(_mutex);
    return _queue.empty();
  }
  void clear() {
    std::lock_guard<std::mutex> lock(_mutex);
    _queue.clear();
  }
  void push(const T& value) {
    std::lock_guard<std::mutex> lock(_mutex);
    _queue.push_back(value);
  }
  bool try_pop(T& value) {
    std::lock_guard<std::mutex> lock(_mutex);
    if (_queue.empty()) return false;
    value = _queue.front();
    _queue.pop_front();
    return true;
  }

 private:
  std::mutex _mutex;
  deque<T>   _queue;
};

// Runs a rask as an asycnrhonous operation.
template <typename Function>
inline auto async(Function&& function) {
  return std::async(std::launch::async, std::forward<Function>(function));
}
template <typename Function, typename... Args>
inline auto async(Function&& function, Args&&... args) {
  return std::async(std::launch::async, std::forward<Function>(function),
      std::forward<Args>(args)...);
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms.
template <typename Func>
inline void parallel_for(size_t begin, size_t end, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
  if (serial) {
    for (auto idx = begin; idx < end; idx++) {
      if (cancel && *cancel) break;
      func(idx);
    }
  } else {
    auto           futures  = vector<future<void>>{};
    auto           nthreads = std::thread::hardware_concurrency();
    atomic<size_t> next_idx(begin);
    for (auto thread_id = 0; thread_id < nthreads; thread_id++) {
      futures.emplace_back(async([&func, &next_idx, cancel, end]() {
        while (true) {
          if (cancel && *cancel) break;
          auto idx = next_idx.fetch_add(1);
          if (idx >= end) break;
          func(idx);
        }
      }));
    }
    for (auto& f : futures) f.get();
  }
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(size_t num, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(0, num, func, cancel, serial);
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms. `Func` takes a reference to a `T`.
template <typename T, typename Func>
inline void parallel_foreach(vector<T>& values, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
      cancel, serial);
}
template <typename T, typename Func>
inline void parallel_foreach(const vector<T>& values, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false) {
  parallel_for(
      0, (int)values.size(), [&func, &values](int idx) { func(values[idx]); },
      cancel, serial);
}

}  // namespace yocto

#endif
