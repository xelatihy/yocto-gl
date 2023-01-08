//
// # Yocto/NDArray: N-dimensional arrays modeled after mdarray/mdspan
//
// Yocto/NDArray is an implementation of a multidimensional array,
// simplified to work well with other Yocto libraries, but mimicking
// the API of mdarray.
// Yocto/NDArray is implemented in `yocto_ndarray.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef YOCTO_NDARRAY_H_
#define YOCTO_NDARRAY_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <stdexcept>
#include <vector>

#include "yocto_math.h"
#include "yocto_views.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifndef kernel
#ifdef __CUDACC__
#define kernel __device__
#else
#define kernel
#endif
#endif

// -----------------------------------------------------------------------------
// ND-ARRAY
// -----------------------------------------------------------------------------
namespace yocto {

// N-dimensional array. We implement specialized versions for simplicity.
template <typename T, size_t N>
struct ndarray {
 public:
  // Constructors
  constexpr ndarray() : _extents{0}, _data{} {}
  constexpr explicit ndarray(const vec<size_t, N>& extents) :
      _extents{extents}, _data(_size(extents), T{}) {}
  template <typename I>
  constexpr explicit ndarray(const vec<I, N>& extents) :
      _extents{(vec<size_t, N>)extents},
      _data(_size((vec<size_t, N>)extents), T{}) {}
  constexpr ndarray(const T* data, const vec<size_t, N>& extents) :
      _extents{extents}, _data(data, data + _size(extents)) {}
  template <typename I>
  constexpr ndarray(const T* data, const vec<I, N>& extents) :
      _extents{(vec<size_t, N>)extents},
      _data(data, data + _size((vec<size_t, N>)extents)) {}
  constexpr ndarray(const ndarray& other) :
      _extents{other._extents}, _data{other._data} {}
  constexpr ndarray(ndarray&& other) : _extents{0}, _data{} {
    std::swap(_extents, other._extents);
    std::swap(_data, other._data);
  }

  // Assignments
  constexpr ndarray& operator=(const ndarray& other) {
    if (&other == this) return *this;
    _extents = other._extents;
    _data    = other._data;
    return *this;
  }
  constexpr ndarray& operator=(ndarray&& other) {
    if (&other == this) return *this;
    std::swap(_extents, other._extents);
    std::swap(_data, other._data);
    return *this;
  }

  // Spans
  constexpr operator ndspan<T, N>() { return {_data.data(), _extents}; }
  constexpr operator ndspan<const T, N>() const {
    return {_data.data(), _extents};
  }

  // Size
  constexpr bool           empty() const { return size() == 0; }
  constexpr size_t         size() const { return _size(_extents); }
  constexpr vec<size_t, N> extents() const { return _extents; }
  constexpr size_t         extent(size_t dimension) const {
    return _extents[dimension];
  }

  // Access
  constexpr T&       operator[](size_t idx) { return _data[idx]; }
  constexpr const T& operator[](size_t idx) const { return _data[idx]; }
  constexpr T&       operator[](const vec<size_t, N>& idx) {
    return _data[_index(idx, _extents)];
  }
  constexpr const T& operator[](const vec<size_t, N>& idx) const {
    return _data[_index(idx, _extents)];
  }
  template <typename I>
  constexpr T& operator[](const vec<I, N>& idx) {
    return _data[_index((vec2s)idx, _extents)];
  }
  template <typename I>
  constexpr const T& operator[](const vec<I, N>& idx) const {
    return _data[_index((vec2s)idx, _extents)];
  }

  // Iteration
  constexpr T*       begin() { return _data.data(); }
  constexpr T*       end() { return _data.data() + size(); }
  constexpr const T* begin() const { return _data.data(); }
  constexpr const T* end() const { return _data.data() + size(); }

  // Data access
  constexpr T*       data() { return _data.data(); }
  constexpr const T* data() const { return _data.data(); }

 private:
  vec<size_t, N> _extents = {0};
  vector<T>      _data    = {};

  static size_t _size(const vec<size_t, 1>& extents) { return extents[0]; }
  static size_t _size(const vec<size_t, 2>& extents) {
    return extents[0] * extents[1];
  }
  static size_t _size(const vec<size_t, 3>& extents) {
    return extents[0] * extents[1] * extents[2];
  }
  static size_t _index(
      const vec<size_t, 1>& index, const vec<size_t, 1>& extents) {
    return index[0];
  }
  static size_t _index(
      const vec<size_t, 2>& index, const vec<size_t, 2>& extents) {
    return index[1] * extents[0] + index[0];
  }
  static size_t _index(
      const vec<size_t, 3>& index, const vec<size_t, 3>& extents) {
    return (index[2] * extents[1] + index[1]) * extents[0] + index[0];
  }
};

// Using directives
template <typename T>
using array1d = ndarray<T, 1>;
template <typename T>
using array2d = ndarray<T, 2>;
template <typename T>
using array3d = ndarray<T, 3>;

// equality
template <typename T, size_t N>
constexpr bool operator==(const ndarray<T, N>& a, const ndarray<T, N>& b) {
  if (a.extents() != b.extents()) return false;
  for (auto idx : range(a.size()))
    if (a.data()[idx] != b.data()[idx]) return false;
  return true;
}
template <typename T, size_t N>
constexpr bool operator!=(const ndarray<T, N>& a, const ndarray<T, N>& b) {
  return !(a == b);
}

// Error handling
template <typename T1, typename T2, size_t N>
constexpr kernel void check_same_size(
    const ndarray<T1, N>& a, const ndarray<T2, N>& b) {
  if (a.extents() != b.extents())
    throw std::out_of_range{"arrays should have the same size"};
}

// Apply a function to each element of an array
template <typename T1, size_t N, typename Func, typename T = result_t<Func, T1>>
constexpr ndarray<T, N> fmap(const ndarray<T1, N>& a, Func&& func) {
  auto ret = ndarray<T, N>(a.extents());
  for (auto idx : range(a.size())) ret[idx] = func(a[idx]);
  return ret;
}
template <typename T, typename T1, size_t N, typename Func>
constexpr void fmap(ndarray<T, N>& ret, const ndarray<T1, N>& a, Func&& func) {
  for (auto idx : range(a.size())) ret[idx] = func(a[idx]);
}

// Python zip
template <typename T1, typename T2, size_t N>
constexpr kernel zip_view<span<const T1>, span<const T2>> zip(
    const ndarray<T1, N>& sequence1, const ndarray<T2, N>& sequence2) {
  check_same_size(sequence1, sequence2);
  return {span<const T1>{sequence1.data(), sequence1.size()},
      span<const T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2, size_t N>
constexpr kernel zip_view<span<const T1>, span<T2>> zip(
    const ndarray<T1, N>& sequence1, ndarray<T2, N>& sequence2) {
  check_same_size(sequence1, sequence2);
  return {span<const T1>{sequence1.data(), sequence1.size()},
      span<T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2, size_t N>
constexpr kernel zip_view<span<T1>, span<const T2>> zip(
    ndarray<T1, N>& sequence1, const ndarray<T2, N>& sequence2) {
  check_same_size(sequence1, sequence2);
  return {span<T1>{sequence1.data(), sequence1.size()},
      span<const T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2, size_t N>
constexpr kernel zip_view<span<T1>, span<T2>> zip(
    ndarray<T1, N>& sequence1, ndarray<T2, N>& sequence2) {
  check_same_size(sequence1, sequence2);
  return {span<T1>{sequence1.data(), sequence1.size()},
      span<T2>{sequence2.data(), sequence2.size()}};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
