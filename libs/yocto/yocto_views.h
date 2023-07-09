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

#ifndef YOCTO_VIEWS_H_
#define YOCTO_VIEWS_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <stdexcept>
#include <tuple>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::tuple;
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
// ONE-DIMENSIONAL SPAN
// -----------------------------------------------------------------------------
namespace yocto {

// Span similar to std::span. We'll switch to the standard library version when
// present in all compilers.
template <typename T>
struct span {
 public:
  // Constructors
  constexpr span() noexcept : _data{nullptr}, _size{0} {}
  constexpr span(const span&) noexcept = default;
  constexpr span(span&&) noexcept      = default;
  constexpr span(T* data, size_t size) noexcept : _data{data}, _size{size} {}
  constexpr span(T* begin, T* end) noexcept :
      _data{begin}, _size{end - begin} {}
  constexpr span(vector<T>& v) : _data{v.data()}, _size(v.size()) {}
  template <typename U>
  constexpr explicit span(const vector<U>& v) :
      _data{v.data()}, _size(v.size()) {}

  // Assignments
  constexpr span& operator=(const span&) noexcept  = default;
  constexpr span& operator=(span&& other) noexcept = default;

  // Size
  constexpr bool   empty() const noexcept { return _size == 0; }
  constexpr size_t size() const noexcept { return _size; }

  // Access
  constexpr T& operator[](size_t idx) const noexcept { return _data[idx]; }
  constexpr T& front() const noexcept { return _data[0]; }
  constexpr T& back() const noexcept { return _data[_size - 1]; }

  // Iteration
  constexpr T* begin() const noexcept { return _data; }
  constexpr T* end() const noexcept { return _data + _size; }

  // Data access
  constexpr T* data() const noexcept { return _data; }

 private:
  T*     _data = nullptr;
  size_t _size = 0;
};

// Constant span
template <typename T>
using cspan = span<const T>;

}  // namespace yocto

// -----------------------------------------------------------------------------
// N-DIMENSIONAL SPAN
// -----------------------------------------------------------------------------
namespace yocto {

// N-dimensional span similar to std::mdspan. We will switch to the standard
// version as it becomes available.
template <typename T, size_t N>
struct ndspan;

// N-dimensional span similar to std::mdspan. We will switch to the standard
// version as it becomes available.
template <typename T>
struct ndspan<T, 2> {
 public:
  // Constructors
  constexpr ndspan() noexcept : _extents{0}, _data{nullptr} {}
  constexpr ndspan(T* data, const vec2i& extents) noexcept :
      _data{data}, _extents{extents} {}
  constexpr ndspan(const ndspan& other) noexcept = default;
  constexpr ndspan(ndspan&& other) noexcept      = default;

  // Assignments
  constexpr ndspan& operator=(const ndspan& other) noexcept = default;
  constexpr ndspan& operator=(ndspan&& other) noexcept      = default;

  // Size
  constexpr bool    empty() const noexcept { return size() == 0; }
  constexpr size_t  size() const noexcept { return _size(_extents); }
  constexpr vec2i   extents() const noexcept { return _extents; }
  constexpr int32_t extent(int dimension) const noexcept {
    return _extents[dimension];
  }
  constexpr vec2i shape() const { return _extents; }
  constexpr int   rank() const { return 2; }

  // Access
  constexpr T& operator[](const vec2i& idx) const noexcept {
    return _data[_index(idx, _extents)];
  }

  // Iteration
  constexpr T* begin() const noexcept { return _data; }
  constexpr T* end() const noexcept { return _data + _size; }

  // Data access
  constexpr T* data() const noexcept { return _data; }

 private:
  T*    _data    = nullptr;
  vec2i _extents = {0, 0};

  static size_t _size(const vec2i& extents) {
    return (size_t)extents[0] * (size_t)extents[1];
  }
  static size_t _index(const vec2i& index, const vec2i& extents) {
    return (size_t)index[1] * (size_t)extents[0] + (size_t)index[0];
  }
};

// Using directives
template <typename T>
using span1d = ndspan<T, 1>;
template <typename T>
using span2d = ndspan<T, 2>;
template <typename T>
using span3d = ndspan<T, 3>;

}  // namespace yocto

// -----------------------------------------------------------------------------
// HELPERS FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Error handling
template <typename T1, typename T2>
constexpr kernel void check_same_size(span<T1> a, span<T2> b) {
  if (a.size() != b.size())
    throw std::out_of_range{"arrays should have the same size"};
}

// Error handling
template <typename T1, typename T2, size_t N>
constexpr kernel void check_same_size(
    const ndspan<T1, N>& a, const ndspan<T2, N>& b) {
  if (a.extents() != b.extents())
    throw std::out_of_range{"arrays should have the same size"};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ONE DIMENSIONAL VIEWS
// -----------------------------------------------------------------------------
namespace yocto {

// Python range: iterator and sequence
template <typename I>
struct range_view {
  struct iterator {
    constexpr kernel      iterator(I index_) : index{index_} {}
    constexpr kernel void operator++() { ++index; }
    constexpr kernel bool operator==(const iterator& other) const {
      return index == other.index;
    }
    constexpr kernel I operator*() const { return index; }

   private:
    I index;
  };
  using sentinel = iterator;

  constexpr kernel          range_view(I max_) : min{0}, max{max_} {}
  constexpr kernel          range_view(I min_, I max_) : min{min_}, max{max_} {}
  constexpr kernel iterator begin() const { return {min}; }
  constexpr kernel sentinel end() const { return {max}; }

 private:
  I min, max;
};

// Python range: iterator and sequence
template <typename I>
struct srange_view {
  struct iterator {
    constexpr kernel iterator(I index_, I step_) : index{index_}, step{step_} {}
    constexpr kernel void operator++() { index += step; }
    constexpr kernel bool operator==(const iterator& other) const {
      return index == other.index;
    }
    constexpr kernel I operator*() const { return index; }

   private:
    I index, step;
  };
  using sentinel = iterator;

  constexpr kernel srange_view(I min_, I max_, I step_) :
      min{min_}, max{max_}, step{step_} {}
  constexpr kernel iterator begin() const { return {min, step}; }
  constexpr kernel sentinel end() const {
    return {min + ((max - min) / step) * step, step};
  }

 private:
  I min, max, step;
};

// Python range. Construct an object that iterates over an integer sequence.
template <typename I>
constexpr kernel range_view<I> range(I max) {
  return range_view<I>(max);
}
template <typename I>
constexpr kernel range_view<I> range(I min, I max) {
  return range_view<I>(min, max);
}
template <typename I>
constexpr kernel srange_view<I> range(I min, I max, I step) {
  return srange_view<I>(min, max, step);
}

// Python enumerate view
template <typename View, typename I>
struct enumerate_view {
  using It = decltype(std::begin(std::declval<View>()));
  using Se = decltype(std::end(std::declval<View>()));
  using Rf = decltype(*std::begin(std::declval<View>()));

  struct iterator {
    constexpr kernel iterator(It cur_, I index_) : cur{cur_}, index{index_} {}
    constexpr kernel bool operator==(const iterator& other) const {
      return cur == other.cur;
    }
    constexpr kernel void operator++() {
      ++cur;
      ++index;
    }
    constexpr kernel tuple<I, Rf> operator*() const { return {index, *cur}; }

   private:
    It cur;
    I  index;
  };
  using sentinel = iterator;

  enumerate_view(View view_) : view{view_}, start{0} {}
  enumerate_view(View view_, I start_) : view{view_}, start{start_} {}
  constexpr kernel iterator begin() { return {std::begin(view), start}; }
  constexpr kernel sentinel end() {
    return {std::end(view), (I)(std::end(view) - std::begin(view))};
  }

 private:
  View view;
  I    start;
};

// Python enumerate over an array
template <typename T, typename I = size_t>
constexpr kernel enumerate_view<span<T>, I> enumerate(
    span<T> sequence, I start = 0) {
  return {sequence, start};
}
template <typename T, typename I = size_t>
constexpr kernel enumerate_view<span<T>, I> enumerate(
    vector<T>& sequence, I start = 0) {
  return {span<T>{sequence.data(), sequence.size()}, start};
}
template <typename T, typename I = size_t>
constexpr kernel enumerate_view<span<const T>, I> enumerate(
    const vector<T>& sequence, I start = 0) {
  return {span<const T>{sequence.data(), sequence.size()}, start};
}
template <typename T, size_t N, typename I = size_t>
constexpr kernel enumerate_view<span<T>, I> enumerate(
    array<T, N>& sequence, I start = 0) {
  return {span<T>{sequence.data(), N}, start};
}
template <typename T, size_t N, typename I = size_t>
constexpr kernel enumerate_view<span<const T>, I> enumerate(
    const array<T, N>& sequence, I start = 0) {
  return {span<const T>{sequence.data(), N}, start};
}

// Python zip: iterator and sequence
template <typename... Views>
struct zip_view {
  static constexpr auto N = sizeof...(Views);
  using ViewT             = tuple<Views...>;
  using ItT = tuple<decltype(std::begin(std::declval<Views>()))...>;
  using SeT = tuple<decltype(std::end(std::declval<Views>()))...>;
  using RfT = tuple<decltype(*std::begin(std::declval<Views>()))...>;

  struct iterator {
    constexpr kernel      iterator(ItT curs_) : curs{curs_} {}
    constexpr kernel bool operator==(const iterator& other) const {
      return curs == other.curs;
    }
    constexpr kernel void operator++() {
      if constexpr (N >= 1) ++get<0>(curs);
      if constexpr (N >= 2) ++get<1>(curs);
      if constexpr (N >= 3) ++get<2>(curs);
      if constexpr (N >= 4) ++get<3>(curs);
    }
    constexpr kernel RfT operator*() const {
      if constexpr (N == 1) return {*get<0>(curs)};
      if constexpr (N == 2) return {*get<0>(curs), *get<1>(curs)};
      if constexpr (N == 3)
        return {*get<0>(curs), *get<1>(curs), *get<2>(curs)};
      if constexpr (N == 4)
        return {*get<0>(curs), *get<1>(curs), *get<2>(curs), *get<3>(curs)};
    }

   private:
    ItT curs;
  };
  using sentinel = iterator;

  constexpr kernel          zip_view(Views... views_) : views{views_...} {}
  constexpr kernel iterator begin() {
    if constexpr (N == 2)
      return {{std::begin(get<0>(views)), std::begin(get<1>(views))}};
    if constexpr (N == 3)
      return {{std::begin(get<0>(views)), std::begin(get<1>(views)),
          std::begin(get<2>(views))}};
    if constexpr (N == 4)
      return {{std::begin(get<0>(views)), std::begin(get<1>(views)),
          std::begin(get<2>(views)), std::begin(get<3>(views))}};
  }
  constexpr kernel sentinel end() {
    if constexpr (N == 2)
      return {{std::end(get<0>(views)), std::end(get<1>(views))}};
    if constexpr (N == 3)
      return {{std::end(get<0>(views)), std::end(get<1>(views)),
          std::end(get<2>(views))}};
    if constexpr (N == 4)
      return {{std::end(get<0>(views)), std::end(get<1>(views)),
          std::end(get<2>(views)), std::end(get<3>(views))}};
  }

 private:
  ViewT views;
};

// Python zip
template <typename T1, typename T2>
constexpr kernel zip_view<span<T1>, span<T2>> zip(
    span<T1> sequence1, span<T2> sequence2) {
  return {sequence1, sequence2};
}
template <typename T1, typename T2>
constexpr kernel zip_view<span<const T1>, span<const T2>> zip(
    const vector<T1>& sequence1, const vector<T2>& sequence2) {
  return {span<const T1>{sequence1.data(), sequence1.size()},
      span<const T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2>
constexpr kernel zip_view<span<const T1>, span<T2>> zip(
    const vector<T1>& sequence1, vector<T2>& sequence2) {
  return {span<const T1>{sequence1.data(), sequence1.size()},
      span<T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2>
constexpr kernel zip_view<span<T1>, span<const T2>> zip(
    vector<T1>& sequence1, const vector<T2>& sequence2) {
  return {span<T1>{sequence1.data(), sequence1.size()},
      span<const T2>{sequence2.data(), sequence2.size()}};
}
template <typename T1, typename T2>
constexpr kernel zip_view<span<T1>, span<T2>> zip(
    vector<T1>& sequence1, vector<T2>& sequence2) {
  return {span<T1>{sequence1.data(), sequence1.size()},
      span<T2>{sequence2.data(), sequence2.size()}};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// ARRAY CREATION VIA VIEWS
// -----------------------------------------------------------------------------
namespace yocto {

// Error handling
template <typename T1, typename T2>
constexpr kernel void check_same_size(
    const vector<T1>& a, const vector<T2>& b) {
  if (a.size() != b.size())
    throw std::out_of_range{"arrays should have the same size"};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MULTI DIMENSIONAL VIEWS
// -----------------------------------------------------------------------------
namespace yocto {

// Range sequence in ND
template <typename I, size_t N, typename O = array<I, N>>
struct ndrange_view {
  struct iterator {
    constexpr kernel iterator(
        const array<I, N>& index_, const array<I, N>& end_) :
        index{index_}, end{end_} {}
    constexpr kernel void operator++() {
      ++index[0];
      if constexpr (N > 1) {
        if (index[0] >= end[0]) {
          index[0] = 0;
          index[1]++;
        }
      }
      if constexpr (N > 2) {
        if (index[1] >= end[1]) {
          index[1] = 0;
          index[2]++;
        }
      }
      if constexpr (N > 3) {
        if (index[2] >= end[2]) {
          index[2] = 0;
          index[3]++;
        }
      }
    }
    constexpr kernel bool operator==(const iterator& other) const {
      return index[N - 1] == other.index[N - 1];
    }
    constexpr kernel O operator*() const { return (O)index; }

   private:
    array<I, N> index, end;
  };
  using sentinel = iterator;

  constexpr kernel          ndrange_view(const array<I, N>& max_) : max{max_} {}
  constexpr kernel iterator begin() const { return {array<I, N>{0}, max}; }
  constexpr kernel sentinel end() const { return {max, max}; }

 private:
  array<I, N> max = {0};
};
// Python range in nd.
template <typename I, size_t N>
constexpr kernel ndrange_view<I, N> range(const array<I, N>& max) {
  return range_sequence<I, N>(max);
}
constexpr kernel ndrange_view<int, 2, vec2i> range(const vec2i& max) {
  return ndrange_view<int, 2, vec2i>(max);
}
constexpr kernel ndrange_view<int, 3, vec3i> range(const vec3i& max) {
  return ndrange_view<int, 3, vec3i>(max);
}
constexpr kernel ndrange_view<int, 4, vec4i> range(const vec4i& max) {
  return ndrange_view<int, 4, vec4i>(max);
}

// Enumerate sequence in ND
template <typename View, typename I, size_t N>
struct ndenumerate_view {
  using It = decltype(std::begin(std::declval<View>()));
  using Se = decltype(std::end(std::declval<View>()));
  using Rf = decltype(*std::begin(std::declval<View>()));

  struct iterator {
    constexpr kernel iterator(
        View view_, const array<I, N>& index_, const array<I, N>& end_) :
        index{index_}, end{end_} {}
    constexpr kernel void operator++() {
      ++cur;
      ++index.x;
      if constexpr (N > 1) {
        if (index[0] >= end[0]) {
          index[0] = 0;
          index[1]++;
        }
      }
      if constexpr (N > 2) {
        if (index[1] >= end[1]) {
          index[1] = 0;
          index[2]++;
        }
      }
      if constexpr (N > 3) {
        if (index[2] >= end[2]) {
          index[2] = 0;
          index[3]++;
        }
      }
    }
    constexpr kernel bool operator==(const iterator& other) const {
      return index[N - 1] == other.index[N - 1];
    }
    constexpr kernel tuple<array<I, N>, Rf> operator*() const {
      return {index, *cur};
    }

   private:
    It          cur;
    array<I, N> index, end;
  };
  using sentinel = iterator;

  constexpr kernel ndenumerate_view(View view_, const array<I, N>& max_) :
      view{view_}, max{max_} {}
  constexpr kernel iterator begin() const {
    return {view, array<I, N>{0}, max};
  }
  constexpr kernel sentinel end() const { return {view, max, max}; }

 private:
  View        view;
  array<I, N> max = {0};
};

// Python enumerate over an array
template <typename T, typename I = size_t, size_t N>
constexpr kernel ndenumerate_view<ndspan<T, N>, I, N> enumerate(
    ndspan<T, N> sequence) {
  return {sequence, sequence.extents()};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
