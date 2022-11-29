//
// # Yocto/Views: Views and ranges
//
// Yocto/Views provides several views and ranges over data typical of graphics,
// applications.
// Yocto/Views is implemented in `yocto_views.h`.
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

#ifndef _YOCTO_VIEWS_H_
#define _YOCTO_VIEWS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <vector>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::vector;

}  // namespace yocto

// -----------------------------------------------------------------------------
// ONE-DIMENSIONAL SPAN
// -----------------------------------------------------------------------------
namespace yocto {

// Span similar to std::span. We'll switch to the standard library version when
// present in all compilers.
template <typename T>
struct span {
 public:
  // Constants
  static constexpr std::size_t extent = dynamic_extent;

  // Constructors
  constexpr span() noexcept : _data{nullptr}, _size{0} {}
  constexpr span(const span&) noexcept = default;
  constexpr span(span&&) noexcept      = default;
  constexpr span(T* data, size_type size) noexcept : _data{data}, _size{size} {}
  template <class It>
  constexpr span(T* begin, T* end) noexcept
      : _data{begin}, _size{end - begin} {}
  template <size_t N>
  constexpr span(std::array<T, N>& arr) noexcept
      : _data{arr.data()}, _size{N} {}

  // Assignments
  constexpr span& operator=(const span&) noexcept  = default;
  constexpr span& operator=(span&& other) noexcept = default;

  // Size
  constexpr bool   empty() const noexcept { return _size == 0; }
  constexpr size_t size() const noexcept { return _size; }

  // Access
  constexpr T& operator[](size_t idx) const noexcept { return _data[idx]; }

  // Iteration
  constexpr T* begin() const noexcept { return _data; }
  constexpr T* end() const noexcept { return _data + _size; }

  // Data access
  constexpr T* data() const noexcept { return _data; }

 private:
  T*     _data = nullptr;
  size_t _size = 0;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// N-DIMENSIONAL SPAN
// -----------------------------------------------------------------------------
namespace yocto {

// N-dimensional span similar to std::mdspan. We will switch to the standard
// version as it becomes available.
template <typename T, size_t N>
struct ndspan {
 public:
  // Constructors
  constexpr ndspan() noexcept: _extents{0}, _data{nullptr} {}
  template <typename... Indices>
  constexpr explicit ndspan(T* data, Indices... extents) noexcept
      : _extents{size_t(extents)...}
      , _data{data} {
    static_assert(N == sizeof...(Indices));
  }
  constexpr ndarray(T* data, const array<T, N>& extents) noexcept
      : _extents{extents}, _data{data} {}
  constexpr ndarray(const ndarray& other) noexcept = default;
  constexpr ndarray(ndarray&& other) noexcept = default;

  // Assignments
  constexpr ndarray& operator=(const ndarray& other) noexcept = default;
  constexpr ndarray& operator=(ndarray&& other) noexcept = default;

  // Size
  constexpr bool             empty() const noexcept { return size() == 0; }
  constexpr size_t           size() const noexcept { return _size(_extents); }
  constexpr array<size_t, N> extents() const noexcept { return _extents; }
  constexpr size_t extent(size_t dimension) const noexcept { return _extents[dimension]; }

  // Access
  constexpr T&       operator[](size_t idx) const noexcept { return _data[idx]; }
  constexpr T&       operator[](const array<size_t, N>& idx) const noexcept {
          return _data[_index(idx, _extents)];
  }

  // Iteration
  constexpr T*       begin()  const noexcept { return _data.data(); }
  constexpr T*       end()  const noexcept { return _data.data() + size(); }

  // Data access
  constexpr T*       data()  const noexcept { return _data.data(); }

 private:
  array<size_t, N> _extents = {0};
  T*        _data    = nullptr;

  static size_t _size(const array<size_t, 1>& extents) { return extents[0]; }
  static size_t _size(const array<size_t, 2>& extents) {
    return extents[0] * extents[1];
  }
  static size_t _size(const array<size_t, 3>& extents) {
    return extents[0] * extents[1] * extents[2];
  }
  static size_t _index(
      const array<size_t, 1>& index, const array<size_t, 1>& extents) {
    return index[0];
  }
  static size_t _index(
      const array<size_t, 2>& index, const array<size_t, 2>& extents) {
    return index[1] * extents[0] + index[0];
  }
  static size_t _index(
      const array<size_t, 3>& index, const array<size_t, 3>& extents) {
    return (index[2] * extents[1] + index[1]) * extents[0] + index[0];
  }
};

// Using directives
template <typename T>
using span1d = ndspan<T, 1>;
template <typename T>
using span2d = ndspan<T, 2>;
template <typename T>
using span3d = ndspan<T, 3>;

}

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

#endif
