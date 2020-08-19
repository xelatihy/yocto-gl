//
// # Yocto/Common: Common utilities
//
// Yocto/Common is a collection of utilities helpful in implementing other
// Yocto/GL libraries. Yocto/Common is implemented in `yocto_common.h`.
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

#ifndef _YOCTO_COMMON_H_
#define _YOCTO_COMMON_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <chrono>
#include <string>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::string;
using std::vector;

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

// Python `range()` equivalent. Construct an instance to iterate over a
// sequence.
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

#endif
