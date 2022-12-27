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
  constexpr span(T* begin, T* end) noexcept
      : _data{begin}, _size{end - begin} {}
  constexpr span(std::vector<T>& arr) noexcept
      : _data{arr.data()}, _size{arr.size()} {}

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
  constexpr ndspan(T* data, const vec<size_t, N>& extents) noexcept
      : _data{data}, _extents{extents} {}
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
  constexpr T& operator[](const array<size_t, N>& idx) const noexcept {
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
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Python range: iterator and sequence
template <typename I>
struct range_sentinel {};
template <typename I>
struct range_iterator {
  constexpr kernel range_iterator(I index_, I end_)
      : index{index_}, end{end_} {}
  constexpr kernel void operator++() { ++index; }
  constexpr kernel bool operator!=(const range_sentinel<I>& other) const {
    return index != end;
  }
  constexpr kernel I operator*() const { return index; }

 private:
  I index, end;
};
template <typename I>
struct range_view {
  constexpr kernel range_view(I max_) : min{0}, max{max_} {}
  constexpr kernel range_view(I min_, I max_) : min{min_}, max{max_} {}
  constexpr kernel range_iterator<I> begin() const { return {min, max}; }
  constexpr kernel range_sentinel<I> end() const { return {}; }

 private:
  I min, max;
};

// Python range: iterator and sequence
template <typename I>
struct srange_sentinel {};
template <typename I>
struct srange_iterator {
  constexpr kernel srange_iterator(I index_, I end_, I step_)
      : index{index_}, end{end_}, step{step_} {}
  constexpr kernel void operator++() { index += step; }
  constexpr kernel bool operator!=(const srange_sentinel<I>& other) const {
    return index != end;
  }
  constexpr kernel I operator*() const { return index; }

 private:
  I index, end, step;
};
// Python range: iterator and sequence
template <typename I>
struct srange_view {
  constexpr kernel srange_view(I min_, I max_, I step_)
      : min{min_}, max{max_}, step{step_} {}
  constexpr kernel srange_iterator<I> begin() const {
    return {min, min + ((max - min) / step) * step, step};
  }
  constexpr kernel srange_sentinel<I> end() const { return {}; }

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
struct ndrange_sentinel {};
template <typename I, size_t N>
struct ndrange_iterator {
  constexpr kernel ndrange_iterator(
      const vec<I, N>& cur_, const vec<I, N>& end_)
      : index{cur_}, end{end_} {}
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
  }
  constexpr kernel bool operator!=(const ndrange_sentinel<I, N>&) const {
    return index[N - 1] != end[N - 1];
  }
  constexpr kernel vec<I, N> operator*() const { return index; }

 private:
  vec<I, N> index, end;
};
template <typename I, size_t N>
struct ndrange_view {
  constexpr kernel ndrange_view(const vec<I, N>& max_) : max{max_} {}
  constexpr kernel ndrange_iterator<I, N> begin() const {
    return {vec<I, N>{0}, max};
  }
  constexpr kernel ndrange_sentinel<I, N> end() const { return {}; }

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
template <typename T, typename I>
struct enumerate_sentinel {};
template <typename T, typename I>
struct enumerate_iterator {
  constexpr kernel enumerate_iterator(I index_, I end_, T* iterator_)
      : index{index_}, end{end_}, iterator{iterator_} {}
  constexpr kernel bool operator!=(
      const enumerate_sentinel<T, I>& other) const {
    return index != end;
  }
  constexpr kernel void operator++() {
    ++index;
    ++iterator;
  }
  constexpr kernel pair<I, T&> operator*() const { return {index, *iterator}; }

 private:
  I  index, end;
  T* iterator;
};
template <typename T, typename I>
struct enumerate_view {
  enumerate_view(span<T> sequence_) : sequence{sequence_}, start{0} {}
  enumerate_view(span<T> sequence_, I start_)
      : sequence{sequence_}, start{start_} {}
  constexpr kernel auto begin() {
    return enumerate_iterator{
        start, (I)sequence.size() + start, sequence.data()};
  }
  constexpr kernel enumerate_sentinel<T, I> end() { return {}; }

 private:
  span<T> sequence = {};
  I       start    = 0;
};

// Python enumerate over an array
template <typename T, typename I = size_t>
constexpr kernel enumerate_view<T, I> enumerate(span<T> sequence, I start = 0) {
  return enumerate_view<T, I>(sequence, start);
}
template <typename T, typename I = size_t>
constexpr kernel enumerate_view<const T, I> enumerate(
    const vector<T>& sequence, I start = 0) {
  return enumerate_view<const T, I>(
      span<const T>{sequence.data(), sequence.size()}, start);
}
template <typename T, size_t N, typename I = size_t>
constexpr kernel enumerate_view<const T, I> enumerate(
    const array<const T, N>& sequence, I start = 0) {
  return enumerate_view<const T, I>(span<const T>{sequence.data(), N}, start);
}
template <typename T, size_t N, typename I = size_t>
constexpr kernel enumerate_view<const T, I> enumerate(
    const vec<T, N>& sequence, I start = 0) {
  return enumerate_view<const T, I>(span<const T>{sequence.data(), N}, start);
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(
    const Sequence1& sequence1, const Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(const Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, const Sequence2& sequence2);

// Implementation of Python enumerate.
template <typename Sequence, typename I>
constexpr kernel auto enumerate(const Sequence& sequence, I start) {
  using Iterator  = typename Sequence::const_iterator;
  using Reference = typename Sequence::const_reference;
  struct enumerate_iterator {
    I                     index;
    Iterator              iterator;
    constexpr kernel bool operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel void operator++() {
      ++index;
      ++iterator;
    }
    constexpr kernel pair<const I&, Reference> operator*() const {
      return {index, *iterator};
    }
  };
  struct enumerate_helper {
    const Sequence&       sequence;
    I                     begin_, end_;
    constexpr kernel auto begin() {
      return enumerate_iterator{begin_, std::begin(sequence)};
    }
    constexpr kernel auto end() {
      return enumerate_iterator{end_, std::end(sequence)};
    }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python enumerate
template <typename Sequence, typename I>
constexpr kernel auto enumerate(Sequence& sequence, I start) {
  using Iterator  = typename Sequence::iterator;
  using Reference = typename Sequence::reference;
  struct enumerate_iterator {
    I                     index;
    Iterator              iterator;
    constexpr kernel bool operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel void operator++() {
      ++index;
      ++iterator;
    }
    constexpr kernel pair<I&, Reference> operator*() const {
      return {index, *iterator};
    }
  };
  struct enumerate_helper {
    Sequence&             sequence;
    I                     begin_, end_;
    constexpr kernel auto begin() {
      return enumerate_iterator{begin_, std::begin(sequence)};
    }
    constexpr kernel auto end() {
      return enumerate_iterator{end_, std::end(sequence)};
    }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(
    const Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1&      sequence1;
    const Sequence2&      sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1&            sequence1;
    Sequence2&            sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(const Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1&      sequence1;
    Sequence2&            sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1&            sequence1;
    const Sequence2&      sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
