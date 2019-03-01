//
// # Yocto/Utils: Tiny collection of utilities to support Yocto/GL
//
//
// Yocto/Utils is a collection of utilities used in writing other Yocto/GL
// libraries and example applications. We support printing builtin and
// Yocto/Math values, simple path manipulation, file lading/saving and basic
// concurrency utilities. These utilities are likely to change often and are to
// be considered internal to Yocto.
//
//
// ## Printing and parsing values
//
// Use `print_value()` to write a string in a stream or `println_values()`
// to print a line of values. Use `format_duraction()` and `format_num()`
// for pretty printing times and numbers. These will change once lib `fmt`
// is accepted in the standard.
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
// ## Path manipulation
//
// We define a few path manipulation utilities to split and join path
// components.
//
// 1. Get paths components with `get_dirname()`, `get_filename()` and
//   `get_extension()`
// 2. Replace the extension with `replace_path_extension()`
// 3. check if a file exists with `exists_file()`
//
//
// ## File IO
//
// 1. load and save text files with `load_text()` and `save_text()`
// 2. load and save binary files with `load_binary()` and `save_binary()`
// 3. use `input_file` and `output_file` as RIIA FILE* wrappers
// 4. use `read_XXX()` and `write_XXX()` to read/write to this files
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

#include <cctype>
#include <chrono>
#include <cstdio>
#include <deque>
#include <future>
#include <mutex>
#include <string>
#include <thread>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

using std::atomic;
using std::deque;
using std::future;
using std::lock_guard;
using std::mutex;
using std::thread;
using namespace std::chrono_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// APPLICATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Log info/error/fatal/trace message
inline void exit_error(const char* msg) {
    printf("%s\n", msg);
    exit(1);
}
inline void exit_error(const string& msg) {
    printf("%s\n", msg.c_str());
    exit(1);
}

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// io error
struct io_error : runtime_error {
    explicit io_error(const char* msg) : runtime_error{msg} {}
    explicit io_error(const std::string& msg) : runtime_error{msg} {}
};

// file inout stream
struct input_file {
    input_file(const string& filename, bool binary = false) {
        this->filename = filename;
        file           = fopen(filename.c_str(), binary ? "rb" : "rt");
        if (!file) throw io_error("could not open " + filename);
    }
    input_file(FILE* fs) {
        file  = fs;
        owned = false;
    }

    input_file(const input_file&) = delete;
    input_file& operator=(const input_file&) = delete;

    ~input_file() {
        if (file && owned) fclose(file);
    }

    string filename = "";
    FILE*  file     = nullptr;
    bool   owned    = true;
};

// file writer
struct output_file {
    output_file(const string& filename, bool binary = false) {
        this->filename = filename;
        file           = fopen(filename.c_str(), binary ? "wb" : "wt");
        if (!file) throw io_error("could not open " + filename);
    }
    output_file(FILE* fs) {
        file  = fs;
        owned = false;
    }

    output_file(const output_file&) = delete;
    output_file& operator=(const output_file&) = delete;

    ~output_file() {
        if (file && owned) fclose(file);
    }

    string filename = "";
    FILE*  file     = nullptr;
    bool   owned    = true;
};

// write a value to a file
template <typename T>
inline void write_value(const output_file& fs, const T& value) {
    if (fwrite(&value, sizeof(value), 1, fs.file) != 1) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// write values to a file
template <typename T>
inline void write_values(const output_file& fs, const vector<T>& values) {
    if (values.empty()) return;
    if (fwrite(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot write to " + fs.filename);
    }
}
template <typename T>
inline void write_values(const output_file& fs, const T* values, size_t count) {
    if (!count) return;
    if (fwrite(values, sizeof(values[0]), count, fs.file) != count) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// write text to a file
inline void write_text(const output_file& fs, const std::string& str) {
    if (fprintf(fs.file, "%s", str.c_str()) < 0) {
        throw io_error("cannot write to " + fs.filename);
    }
}

// read a value from a file
template <typename T>
inline void read_value(const input_file& fs, T& value) {
    if (fread(&value, sizeof(value), 1, fs.file) != 1) {
        throw io_error("cannot read from " + fs.filename);
    }
}

// read values from a file
template <typename T>
inline void read_values(const input_file& fs, T* values, size_t count) {
    if (!count) return;
    if (fread(values, sizeof(values[0]), count, fs.file) != count) {
        throw io_error("cannot read from " + fs.filename);
    }
}
template <typename T>
inline void read_values(const input_file& fs, vector<T>& values) {
    if (values.empty()) return;
    if (fread(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot read from " + fs.filename);
    }
}
// read characters from a file
inline void read_values(const input_file& fs, string& values) {
    if (values.empty()) return;
    if (fread(values.data(), sizeof(values[0]), values.size(), fs.file) !=
        values.size()) {
        throw io_error("cannot read from " + fs.filename);
    }
}

// read a line of text
inline bool read_line(const input_file& fs, string& str) {
    char buffer[4096];
    if (fgets(buffer, sizeof(buffer), fs.file) == nullptr) return false;
    str = buffer;
    return true;
}
inline bool read_line(const input_file& fs, char* buffer, size_t size) {
    if (fgets(buffer, size, fs.file) == nullptr) return false;
    return true;
}

// Printing values
inline void print_value(const output_file& fs, int value) {
    if (fprintf(fs.file, "%d", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, bool value) {
    if (fprintf(fs.file, "%d", (int)value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, float value) {
    if (fprintf(fs.file, "%g", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, char value) {
    if (fprintf(fs.file, "%c", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, const char* value) {
    if (fprintf(fs.file, "%s", value) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
inline void print_value(const output_file& fs, const string& value) {
    if (fprintf(fs.file, "%s", value.c_str()) < 0)
        throw io_error("cannot write to file " + fs.filename);
}
template <typename T, int N>
inline void print_value(const output_file& fs, const vec<T, N>& value) {
    for (auto i = 0; i < N; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}
template <typename T, int N, int M>
inline void print_value(const output_file& fs, const mat<T, N, M>& value) {
    for (auto i = 0; i < M; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}
template <typename T, int N>
inline void print_value(const output_file& fs, const frame<T, N>& value) {
    for (auto i = 0; i < N + 1; i++) {
        if (i) print_value(fs, ' ');
        print_value(fs, value[i]);
    }
}

// print values to file
template <typename Arg, typename... Args>
inline void println_values(
    const output_file& fs, const Arg& value, const Args&... values) {
    print_value(fs, value);
    if constexpr (sizeof...(values) > 0) {
        print_value(fs, ' ');
        println_values(fs, values...);
    } else {
        print_value(fs, '\n');
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Range helpper (this should not be used directly)
struct _range_helper {
    struct _iterator {
        int        _pos = 0;
        _iterator& operator++() {
            _pos++;
            return *this;
        }
        bool operator!=(const _iterator& other) const {
            return _pos != other._pos;
        }
        int operator*() const { return _pos; }
    };
    int       _start = 0, _end = 0;
    _iterator begin() const { return {_start}; }
    _iterator end() const { return {_end}; }
};

// Python `range()` equivalent. Construct an object to iterate over a sequence.
inline auto range(int max) { return _range_helper{0, max}; }
inline auto range(int min, int max) { return _range_helper{min, max}; }

// Enumerate helper (this should not be used directly)
template <typename T>
struct _enumerate_helper {
    struct _iterator {
        T*         _data = nullptr;
        int        _pos  = 0;
        _iterator& operator++() {
            _pos++;
            return *this;
        }
        bool operator!=(const _iterator& other) const {
            return _pos != other._pos;
        }
        pair<int&, T&> operator*() const { return {_pos, *(_data + _pos)}; }
    };
    T*        _data = nullptr;
    int       _size = 0;
    _iterator begin() const { return {_data, 0}; }
    _iterator end() const { return {_data, _size}; }
};

// Python `enumerate()` equivalent. Construct an object that iteraterates over a
// sequence of elements and numbers them.
template <typename T>
inline auto enumerate(const vector<T>& vals) {
    return _enumerate_helper<const T>{vals.data(), vals.size()};
};
template <typename T>
inline auto enumerate(vector<T>& vals) {
    return _enumerate_helper<T>{vals.data(), vals.size()};
};

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
// PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Normalize path delimiters.
inline string normalize_path(const string& filename);
// Get directory name (including '/').
inline string get_dirname(const string& filename);
// Get extension (not including '.').
inline string get_extension(const string& filename);
// Get filename without directory.
inline string get_filename(const string& filename);
// Replace extension.
inline string replace_extension(const string& filename, const string& ext);

// Check if a file can be opened for reading.
inline bool exists_file(const string& filename);

}  // namespace yocto

// -----------------------------------------------------------------------------
// FILE IO
// -----------------------------------------------------------------------------
namespace yocto {

// Load/save a text file
inline void load_text(const string& filename, string& str);
inline void save_text(const string& filename, const string& str);

// Load/save a binary file
inline void load_binary(const string& filename, vector<byte>& data);
inline void save_binary(const string& filename, const vector<byte>& data);

}  // namespace yocto

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
    concurrent_queue();
    concurrent_queue(const concurrent_queue& other);
    concurrent_queue& operator=(const concurrent_queue& other);

    bool empty();
    void clear();
    void push(const T& value);
    bool try_pop(T& value);

   private:
    mutex    _mutex;
    deque<T> _queue;
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
// parallel algorithms. `Func` takes the integer index.
template <typename Func>
inline void parallel_for(int begin, int end, const Func& func,
    atomic<bool>* cancel = nullptr, bool serial = false);
template <typename Func>
inline void parallel_for(int num, const Func& func,
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

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING FORMAT UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Format duration string from nanoseconds
inline string format_duration(int64_t duration) {
    auto elapsed = duration / 1000000;  // milliseconds
    auto hours   = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs  = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buffer[256];
    sprintf(buffer, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    return buffer;
}
// Format a large integer number in human readable form
inline string format_num(uint64_t num) {
    auto rem = num % 1000;
    auto div = num / 1000;
    if (div > 0) return format_num(div) + "," + std::to_string(rem);
    return std::to_string(rem);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF PATH UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

string normalize_path(const string& filename_) {
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
string get_dirname(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(0, pos + 1);
}

// Get extension (not including '.').
string get_extension(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('.');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Get filename without directory.
string get_filename(const string& filename_) {
    auto filename = normalize_path(filename_);
    auto pos      = filename.rfind('/');
    if (pos == string::npos) return "";
    return filename.substr(pos + 1);
}

// Replace extension.
string replace_extension(const string& filename_, const string& ext_) {
    auto filename = normalize_path(filename_);
    auto ext      = normalize_path(ext_);
    if (ext.at(0) == '.') ext = ext.substr(1);
    auto pos = filename.rfind('.');
    if (pos == string::npos) return filename;
    return filename.substr(0, pos) + "." + ext;
}

// Check if a file can be opened for reading.
bool exists_file(const string& filename) {
    auto f = fopen(filename.c_str(), "r");
    if (!f) return false;
    fclose(f);
    return true;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF FILE READING
// -----------------------------------------------------------------------------
namespace yocto {

// Load a text file
inline void load_text(const string& filename, string& str) {
    // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
    auto fs = input_file(filename);
    fseek(fs.file, 0, SEEK_END);
    auto length = ftell(fs.file);
    fseek(fs.file, 0, SEEK_SET);
    str.resize(length);
    read_values(fs, str);
}

// Save a text file
inline void save_text(const string& filename, const string& str) {
    auto fs = output_file(filename);
    write_text(fs, str);
}

// Load a binary file
inline void load_binary(const string& filename, vector<byte>& data) {
    // https://stackoverflow.com/questions/174531/how-to-read-the-content-of-a-file-to-a-string-in-c
    auto fs = input_file(filename, true);
    fseek(fs.file, 0, SEEK_END);
    auto length = ftell(fs.file);
    fseek(fs.file, 0, SEEK_SET);
    data.resize(length);
    read_values(fs, data);
}

// Save a binary file
inline void save_binary(const string& filename, const vector<byte>& data) {
    auto fs = output_file(filename, true);
    write_values(fs, data);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// a simple concurrent queue that locks at every call
template <typename T>
inline concurrent_queue<T>::concurrent_queue() {}
template <typename T>
inline concurrent_queue<T>::concurrent_queue(const concurrent_queue<T>& other) {
    if (!other._queue.empty())
        throw std::invalid_argument("cannot copy full queue");
    clear();
}
template <typename T>
inline concurrent_queue<T>& concurrent_queue<T>::operator=(
    const concurrent_queue<T>& other) {
    if (!other._queue.empty())
        throw std::invalid_argument("cannot copy full queue");
    clear();
}

template <typename T>
inline bool concurrent_queue<T>::empty() {
    lock_guard<mutex> lock(_mutex);
    return _queue.empty();
}
template <typename T>
inline void concurrent_queue<T>::clear() {
    lock_guard<mutex> lock(_mutex);
    _queue.clear();
}
template <typename T>
inline void concurrent_queue<T>::push(const T& value) {
    lock_guard<mutex> lock(_mutex);
    _queue.push_back(value);
}
template <typename T>
inline bool concurrent_queue<T>::try_pop(T& value) {
    lock_guard<mutex> lock(_mutex);
    if (_queue.empty()) return false;
    value = _queue.front();
    _queue.pop_front();
    return true;
}

// Simple parallel for used since our target platforms do not yet support
// parallel algorithms.
template <typename Func>
inline void parallel_for(
    int begin, int end, const Func& func, atomic<bool>* cancel, bool serial) {
    if (serial) {
        for (auto idx = begin; idx < end; idx++) {
            if (cancel && *cancel) break;
            func(idx);
        }
    } else {
        auto        futures  = vector<future<void>>{};
        auto        nthreads = thread::hardware_concurrency();
        atomic<int> next_idx(begin);
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

}  // namespace yocto

#endif
