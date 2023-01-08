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

#include <algorithm>
#include <stdexcept>
#include <vector>

#include "yocto_math.h"

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
struct ndspan {
 public:
  // Constructors
  constexpr ndspan() noexcept : _extents{0}, _data{nullptr} {}
  constexpr ndspan(T* data, const vec<size_t, N>& extents) noexcept :
      _data{data}, _extents{extents} {}
  constexpr ndspan(const ndspan& other) noexcept = default;
  constexpr ndspan(ndspan&& other) noexcept      = default;

  // Assignments
  constexpr ndspan& operator=(const ndspan& other) noexcept = default;
  constexpr ndspan& operator=(ndspan&& other) noexcept      = default;

  // Size
  constexpr bool           empty() const noexcept { return size() == 0; }
  constexpr size_t         size() const noexcept { return _size(_extents); }
  constexpr vec<size_t, N> extents() const noexcept { return _extents; }
  constexpr size_t         extent(size_t dimension) const noexcept {
    return _extents[dimension];
  }

  // Access
  constexpr T& operator[](size_t idx) const noexcept { return _data[idx]; }
  constexpr T& operator[](const vec<size_t, N>& idx) const noexcept {
    return _data[_index(idx, _extents)];
  }

  // Iteration
  constexpr T* begin() const noexcept { return _data; }
  constexpr T* end() const noexcept { return _data + _size; }

  // Data access
  constexpr T* data() const noexcept { return _data; }

 private:
  T*             _data    = nullptr;
  vec<size_t, N> _extents = {0};

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

#ifndef __CUDACC__

// Make a vector from a view
template <typename R, typename T = rvalue_t<R>>
inline vector<T> to_vector(R&& range) {
  auto values = vector<T>{};
  for (auto value : range) values.push_back(value);
  return values;
}
template <typename R, typename Func, typename T = result_t<Func, rvalue_t<R>>>
inline vector<T> to_vector(R&& range, Func&& func) {
  auto values = vector<T>{};
  for (auto value : range) values.push_back(func(value));
  return values;
}

#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// ARRAY SEARCH AND SORT
// -----------------------------------------------------------------------------
namespace yocto {

// Find an element with linear search
template <typename T>
inline ptrdiff_t find_index(const vector<T>& values, const T& value) {
  auto pos = std::find(values.begin(), values.end(), value);
  if (pos == values.end()) return -1;
  return pos - values.begin();
}

// Find an element with binary search
template <typename T>
inline ptrdiff_t search_index(const vector<T>& values, const T& value) {
  auto pos = std::binary_search(values.begin(), values.end(), value);
  if (pos == values.end()) return -1;
  return pos - values.begin();
}

// Sort elements in an array
template <typename T>
inline void sort(vector<T>& values) {
  std::sort(values.begin(), values.end());
}
template <typename T>
inline vector<T> sorted(const vector<T>& values) {
  auto sorted = values;
  std::sort(sorted.begin(), sorted.end());
  return sorted;
}

// Sort and remove duplicates
template <typename T>
inline vector<T> remove_duplicates(const vector<T>& values_) {
  auto values = values_;
  std::sort(values.begin(), values.end());
  auto pos = std::unique(values.begin(), values.end());
  values.erase(pos, values.end());
  return values;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Python range: iterator and sequence
template <typename I>
struct range_view {
  struct range_iterator;
  struct range_sentinel {
    constexpr kernel   range_sentinel(I end_) : end{end_} {}
    constexpr kernel I sentinel() const { return end; }
    friend struct range_iterator;

   private:
    I end;
  };

  struct range_iterator {
    constexpr kernel      range_iterator(I index_) : current{index_} {}
    constexpr kernel I    index() const { return current; }
    constexpr kernel void operator++() { ++current; }
    constexpr kernel bool operator!=(const range_sentinel& other) const {
      return current != other.end;
    }
    constexpr kernel I operator*() const { return current; }
    friend struct range_sentinel;

   private:
    I current;
  };

  constexpr kernel range_view(I max_) : min{0}, max{max_} {}
  constexpr kernel range_view(I min_, I max_) : min{min_}, max{max_} {}
  constexpr kernel range_iterator begin() const { return {min}; }
  constexpr kernel range_sentinel end() const { return {max}; }

 private:
  I min, max;
};

// Python range: iterator and sequence
template <typename I>
struct srange_view {
  struct srange_iterator;
  struct srange_sentinel {
    constexpr kernel srange_sentinel(I end_) : end{end_} {}
    friend struct srange_iterator;

   private:
    I end;
  };
  struct srange_iterator {
    constexpr kernel srange_iterator(I index_, I step_) :
        index{index_}, step{step_} {}
    constexpr kernel void operator++() { index += step; }
    constexpr kernel bool operator!=(const srange_sentinel& other) const {
      return index != other.end;
    }
    constexpr kernel I operator*() const { return index; }

   private:
    I index, step;
  };

  constexpr kernel srange_view(I min_, I max_, I step_) :
      min{min_}, max{max_}, step{step_} {}
  constexpr kernel srange_iterator begin() const { return {min, step}; }
  constexpr kernel srange_sentinel end() const {
    return {min + ((max - min) / step) * step};
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

// Python range: iterator and sequence in 2D
template <typename I, size_t N>
struct ndrange_view {
  struct iterator;
  struct sentinel {
    constexpr kernel sentinel(const vec<I, N>& end_) : end{end_} {}
    friend struct iterator;

   private:
    vec<I, N> index, end;
  };
  struct iterator {
    constexpr kernel iterator(const vec<I, N>& cur_, const vec<I, N>& end_) :
        index{cur_}, end{end_} {}
    constexpr kernel void operator++() {
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
    constexpr kernel bool operator!=(const sentinel& other) const {
      return index[N - 1] != other.end[N - 1];
    }
    constexpr kernel vec<I, N> operator*() const { return index; }

   private:
    vec<I, N> index, end;
  };

  constexpr kernel          ndrange_view(const vec<I, N>& max_) : max{max_} {}
  constexpr kernel iterator begin() const { return {vec<I, N>{0}, max}; }
  constexpr kernel sentinel end() const { return {max}; }

 private:
  vec<I, N> max = {0};
};

// Python range in nd.
template <typename I, size_t N>
constexpr kernel ndrange_view<I, N> range(const vec<I, N>& max) {
  return ndrange_view<I, N>(max);
}
template <typename I, size_t N>
constexpr kernel ndrange_view<I, N> range(const array<I, N>& max) {
  return range_sequence<I, N>((vec<I, N>)max);
}

// Python enumerate view
template <typename View, typename I>
struct enumerate_view {
  using It = decltype(std::begin(std::declval<View>()));
  using Se = decltype(std::end(std::declval<View>()));
  using Rf = decltype(*std::begin(std::declval<View>()));

  struct iterator;
  struct sentinel {
    constexpr kernel sentinel(Se end_) : end{end_} {}
    friend struct iterator;

   private:
    Se end;
  };
  struct iterator {
    constexpr kernel iterator(It cur_, I index_) : cur{cur_}, index{index_} {}
    constexpr kernel bool operator!=(const sentinel& other) const {
      return cur != other.end;
    }
    constexpr kernel void operator++() {
      ++cur;
      ++index;
    }
    constexpr kernel pair<I, Rf> operator*() const { return {index, *cur}; }

   private:
    It cur;
    I  index;
  };

  enumerate_view(View view_) : view{view_}, start{0} {}
  enumerate_view(View view_, I start_) : view{view_}, start{start_} {}
  constexpr kernel iterator begin() { return {std::begin(view), start}; }
  constexpr kernel sentinel end() { return {std::end(view)}; }

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
template <typename T, size_t N, typename I = size_t>
constexpr kernel enumerate_view<span<const T>, I> enumerate(
    const vec<T, N>& sequence, I start = 0) {
  return {span<const T>{sequence.data(), N}, start};
}

// Python zip: iterator and sequence
template <typename View1, typename View2>
struct zip_view {
  using It1 = decltype(std::begin(std::declval<View1>()));
  using Se1 = decltype(std::end(std::declval<View1>()));
  using Rf1 = decltype(*std::begin(std::declval<View1>()));
  using It2 = decltype(std::begin(std::declval<View2>()));
  using Se2 = decltype(std::end(std::declval<View2>()));
  using Rf2 = decltype(*std::begin(std::declval<View2>()));

  struct iterator;
  struct sentinel {
    constexpr kernel sentinel(Se1 end1_, Se2 end2_) :
        end1{end1_}, end2{end2_} {}
    friend struct iterator;

   private:
    Se1 end1;
    Se2 end2;
  };
  struct iterator {
    constexpr kernel iterator(It1 cur1_, It2 cur2_) :
        cur1{cur1_}, cur2{cur2_} {}
    constexpr kernel bool operator!=(const sentinel& other) const {
      return cur1 != other.end1 && cur2 != other.end2;
    }
    constexpr kernel void operator++() {
      ++cur1;
      ++cur2;
    }
    constexpr kernel pair<Rf1, Rf2> operator*() const { return {*cur1, *cur2}; }

   private:
    It1 cur1;
    It2 cur2;
  };

  constexpr kernel zip_view(View1 view1_, View2 view2_) :
      view1{view1_}, view2{view2_} {}
  constexpr kernel iterator begin() {
    return {std::begin(view1), std::begin(view2)};
  }
  constexpr kernel sentinel end() { return {std::end(view1), std::end(view2)}; }

 private:
  View1 view1;
  View2 view2;
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
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
