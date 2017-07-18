# Yocto/Math

A collection of vector math functions and simple containers
used to implement YOCTO. Features include

- a few convenience math functions
- static length float vectors, with specialization for 2, 3, 4 length
- static length matrices, with specialization for 2x2, 3x3, 4x4
- static length rigid transforms (frames), specialized for 2d and 3d space
- linear algebra operations and transforms for fixed length matrices/vecs
- axis aligned bounding boxes
- rays
- ray-primitive intersection
- point-primitive distance and overlap tests
- normal amd tangent computation for meshes and lines
- generation of tesselated meshes
- random number generation via PCG32
- a few hash functions
- trivial image data structue and a few image operations

We developed our own library since we felt that all existing ones are either
complete, but unreadable or with lots of dependencies, or just as incomplete
and untested as ours.

This library has no dependencies.

This library includes code from the PCG random number generator,
boost hash_combine, Pixar multijittered sampling, code from "Real-Time
Collision Detection" by Christer Ericson and public domain code from
- https://github.com/sgorsten/linalg
- https://gist.github.com/badboy/6267743


## History

- v 0.18: bump to normal mapping convertion
- v 0.17: added example image geneation
- v 0.16: sampling
- v 0.15: enable specialization always
- v 0.14: move timer to Yocto/Utils
- v 0.13: more shape functions
- v 0.12: documentation update
- v 0.11: added more matrix and quaternion operations
- v 0.10: specialize some type and functions
- v 0.9: bbox containment tests
- v 0.8: remove std:array as base class for better control
- v 0.7: doxygen comments
- v 0.6: uniformed internal names
- v 0.5: simplification of constructors, raname bbox -> bbox
- v 0.4: overall type simplification
- v 0.3: internal C++ refactoring
- v 0.2: use of STL containers; removal of yocto containers
- v 0.1: C++ only implementation
- v 0.0: initial release in C99

## Namespace ym

Math types and utlities for 3D graphics and imaging

### Typedef byte

~~~ .cpp
using byte = unsigned char;
~~~

convenient typedef for bytes

### Typedef uint

~~~ .cpp
using uint = unsigned int;
~~~

convenient typedef for bytes

### Constant pif

~~~ .cpp
const float pif = 3.14159265f;
~~~

pi (float)

### Constant pi

~~~ .cpp
const double pi = 3.1415926535897932384626433832795;
~~~

pi (double)

### Constant flt_max

~~~ .cpp
constexpr const auto flt_max = std::numeric_limits<float>::max();
~~~

shortcat for float max value

### Constant flt_min

~~~ .cpp
constexpr const auto flt_min = std::numeric_limits<float>::lowest();
~~~

shortcat for float min value

### Constant int_max

~~~ .cpp
constexpr const auto int_max = std::numeric_limits<int>::max();
~~~

shortcat for int max value

### Constant int_min

~~~ .cpp
constexpr const auto int_min = std::numeric_limits<int>::min();
~~~

shortcat for int min value

### Function min()

~~~ .cpp
template <typename T>
constexpr inline T min(T x, T y);
~~~

Safe minimum value.

### Function max()

~~~ .cpp
template <typename T>
constexpr inline T max(T x, T y);
~~~

Safe maximum value.

### Function clamp()

~~~ .cpp
template <typename T>
constexpr inline T clamp(T x, T min_, T max_);
~~~

Clamp a value between a minimum and a maximum.

### Function lerp()

~~~ .cpp
template <typename T>
constexpr inline T lerp(T a, T b, T t);
~~~

Linear interpolation.

### Function pow2()

~~~ .cpp
constexpr inline int pow2(int x);
~~~

Integer power of two

### Function float_to_byte()

~~~ .cpp
constexpr inline byte float_to_byte(float x);
~~~

Safe float to byte conversion

### Function byte_to_float()

~~~ .cpp
constexpr inline float byte_to_float(byte x);
~~~

Safe byte to float conversion

### Function Alias sqrt()

~~~ .cpp
using std::sqrt;
~~~

sqrt

### Function Alias pow()

~~~ .cpp
using std::pow;
~~~

pow

### Function Alias sin()

~~~ .cpp
using std::sin;
~~~

sin

### Function Alias cos()

~~~ .cpp
using std::cos;
~~~

cos

### Function Alias tan()

~~~ .cpp
using std::tan;
~~~

tan

### Function Alias asin()

~~~ .cpp
using std::asin;
~~~

asin

### Function Alias acos()

~~~ .cpp
using std::acos;
~~~

acos

### Function Alias atan2()

~~~ .cpp
using std::atan2;
~~~

atan2

### Function Alias abs()

~~~ .cpp
using std::abs;
~~~

abs

### Struct vec

~~~ .cpp
template <typename T, int N>
struct vec {
    constexpr vec(); 
    constexpr vec(const std::initializer_list<T>& vv); 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    T v[N];
}
~~~

Vector of elements of compile time dimension with default initializer.

- Members:
    - vec():      default constructor
    - vec():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - v[N]:      element data


### Struct vec <T, 1 \>

~~~ .cpp
template <typename T>
struct vec<T, 1> {
    constexpr static const int N = 1;
    constexpr vec(); 
    constexpr vec(T x); 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    T x;
}
~~~

Specialization of vectors for 1 component and float coordinates.

- Members:
    - N:      size
    - vec():      default constructor
    - vec():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data


### Struct vec <T, 2 \>

~~~ .cpp
template <typename T>
struct vec<T, 2> {
    constexpr static const int N = 2;
    constexpr vec(); 
    constexpr vec(T x, T y); 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    T x;
    T y;
}
~~~

Specialization of vectors for 2 components and float coordinates.

- Members:
    - N:      size
    - vec():      default constructor
    - vec():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data


### Struct vec <T, 3 \>

~~~ .cpp
template <typename T>
struct vec<T, 3> {
    constexpr static const int N = 3;
    constexpr vec(); 
    constexpr vec(T x, T y, T z); 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    T x;
    T y;
    T z;
}
~~~

Specialization of vectors for 3 components and float coordinates.

- Members:
    - N:      size
    - vec():      default constructor
    - vec():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


### Struct vec <T, 4 \>

~~~ .cpp
template <typename T>
struct vec<T, 4> {
    constexpr static const int N = 4;
    constexpr vec(); 
    constexpr vec(T x, T y, T z, T w); 
    constexpr vec(const vec<T, 3>& xyz, T w); 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    constexpr vec<T, 3>& xyz(); 
    constexpr const vec<T, 3>& xyz() const; 
    T x;
    T y;
    T z;
    T w;
}
~~~

Specialization of vectors for 4 components and float coordinates.

- Members:
    - N:      size
    - vec():      default constructor
    - vec():      element constructor
    - vec():      constructor from smaller vector
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - xyz():      access xyz components
    - xyz():      access xyz components
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


### Typedef vec1f

~~~ .cpp
using vec1f = vec<float, 1>;
~~~

1-dimensional float vector

### Typedef vec2f

~~~ .cpp
using vec2f = vec<float, 2>;
~~~

2-dimensional float vector

### Typedef vec3f

~~~ .cpp
using vec3f = vec<float, 3>;
~~~

3-dimensional float vector

### Typedef vec4f

~~~ .cpp
using vec4f = vec<float, 4>;
~~~

4-dimensional float vector

### Typedef vec1i

~~~ .cpp
using vec1i = vec<int, 1>;
~~~

1-dimensional int vector

### Typedef vec2i

~~~ .cpp
using vec2i = vec<int, 2>;
~~~

2-dimensional int vector

### Typedef vec3i

~~~ .cpp
using vec3i = vec<int, 3>;
~~~

3-dimensional int vector

### Typedef vec4i

~~~ .cpp
using vec4i = vec<int, 4>;
~~~

4-dimensional int vector

### Typedef vec1b

~~~ .cpp
using vec1b = vec<byte, 1>;
~~~

1-dimensional byte vector

### Typedef vec2b

~~~ .cpp
using vec2b = vec<byte, 2>;
~~~

2-dimensional byte vector

### Typedef vec3b

~~~ .cpp
using vec3b = vec<byte, 3>;
~~~

3-dimensional byte vector

### Typedef vec4b

~~~ .cpp
using vec4b = vec<byte, 4>;
~~~

4-dimensional byte vector

### Function zero_vec()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> zero_vec();
~~~

Initialize a zero vector.

### Function zero_vec()

~~~ .cpp
template <>
constexpr inline vec3f zero_vec();
~~~

Sepcialization of Initialize a zero vector.

### Constant zero1f

~~~ .cpp
const auto zero1f = zero_vec<float, 1>();
~~~

1-dimensional float zero vector

### Constant zero2f

~~~ .cpp
const auto zero2f = zero_vec<float, 2>();
~~~

2-dimensional float zero vector

### Constant zero3f

~~~ .cpp
const auto zero3f = zero_vec<float, 3>();
~~~

3-dimensional float zero vector

### Constant zero4f

~~~ .cpp
const auto zero4f = zero_vec<float, 4>();
~~~

4-dimensional float zero vector

### Constant zero1i

~~~ .cpp
const auto zero1i = zero_vec<int, 1>();
~~~

1-dimensional int zero vector

### Constant zero2i

~~~ .cpp
const auto zero2i = zero_vec<int, 2>();
~~~

2-dimensional int zero vector

### Constant zero3i

~~~ .cpp
const auto zero3i = zero_vec<int, 3>();
~~~

3-dimensional int zero vector

### Constant zero4i

~~~ .cpp
const auto zero4i = zero_vec<int, 4>();
~~~

4-dimensional int zero vector

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline T* begin(vec<T, N>& a);
~~~

iteration support

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline const T* begin(const vec<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline T* end(vec<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline const T* end(const vec<T, N>& a);
~~~

iteration support

### Function operator==()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator==(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator !=

### Function operator <()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator<(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator < (lexicographic order - useful for std::map)

### Function operator==()

~~~ .cpp
template <>
constexpr inline bool operator==(const vec2f& a, const vec2f& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <>
constexpr inline bool operator!=(const vec2f& a, const vec2f& b);
~~~

vector operator !=

### Function operator==()

~~~ .cpp
template <>
constexpr inline bool operator==(const vec3f& a, const vec3f& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <>
constexpr inline bool operator!=(const vec3f& a, const vec3f& b);
~~~

vector operator !=

### Function operator==()

~~~ .cpp
template <>
constexpr inline bool operator==(const vec4f& a, const vec4f& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <>
constexpr inline bool operator!=(const vec4f& a, const vec4f& b);
~~~

vector operator !=

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, const T b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, const T b);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator+(T a, const vec<T, N>& b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator-(T a, const vec<T, N>& b);
~~~

vector operator -

### Function operator*()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, const T b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator*(const T a, const vec<T, N>& b);
~~~

vector operator *

### Function operator/()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, const T b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> operator/(const T a, const vec<T, N>& b);
~~~

vector operator /

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec2f operator+(const vec2f& a);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec2f operator-(const vec2f& a);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec2f operator+(const vec2f& a, const vec2f& b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec2f operator-(const vec2f& a, const vec2f& b);
~~~

vector operator -

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec2f operator*(const vec2f& a, const vec2f& b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec2f operator*(const vec2f& a, const float b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec2f operator*(const float a, const vec2f& b);
~~~

vector operator *

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec2f operator/(const vec2f& a, const vec2f& b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec2f operator/(const vec2f& a, const float b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec2f operator/(const float a, const vec2f& b);
~~~

vector operator /

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec3f operator+(const vec3f& a);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec3f operator-(const vec3f& a);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec3f operator+(const vec3f& a, const vec3f& b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec3f operator-(const vec3f& a, const vec3f& b);
~~~

vector operator -

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec3f operator*(const vec3f& a, const vec3f& b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec3f operator*(const vec3f& a, const float b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec3f operator*(const float a, const vec3f& b);
~~~

vector operator *

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec3f operator/(const vec3f& a, const vec3f& b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec3f operator/(const vec3f& a, const float b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec3f operator/(const float a, const vec3f& b);
~~~

vector operator /

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec4f operator+(const vec4f& a);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec4f operator-(const vec4f& a);
~~~

vector operator -

### Function operator+()

~~~ .cpp
template <>
constexpr inline vec4f operator+(const vec4f& a, const vec4f& b);
~~~

vector operator +

### Function operator-()

~~~ .cpp
template <>
constexpr inline vec4f operator-(const vec4f& a, const vec4f& b);
~~~

vector operator -

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec4f operator*(const vec4f& a, const vec4f& b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec4f operator*(const vec4f& a, const float b);
~~~

vector operator *

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec4f operator*(const float a, const vec4f& b);
~~~

vector operator *

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec4f operator/(const vec4f& a, const vec4f& b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec4f operator/(const vec4f& a, const float b);
~~~

vector operator /

### Function operator/()

~~~ .cpp
template <>
constexpr inline vec4f operator/(const float a, const vec4f& b);
~~~

vector operator /

### Function operator+=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator +=

### Function operator-=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator -=

### Function operator*=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator *=

### Function operator*=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, const T b);
~~~

vector operator *=

### Function operator/=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b);
~~~

vector operator /=

### Function operator/=()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, const T b);
~~~

vector operator /=

### Function dot()

~~~ .cpp
template <typename T, int N>
constexpr inline T dot(const vec<T, N>& a, const vec<T, N>& b);
~~~

vector dot product

### Function cross()

~~~ .cpp
template <typename T>
constexpr inline T cross(const vec<T, 2>& a, const vec<T, 2>& b);
~~~

vector cross product (2d)

### Function cross()

~~~ .cpp
template <typename T>
constexpr inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b);
~~~

vector cross product (3d)

### Function dot()

~~~ .cpp
template <>
constexpr inline float dot(const vec2f& a, const vec2f& b);
~~~

vector dot product

### Function dot()

~~~ .cpp
template <>
constexpr inline float dot(const vec3f& a, const vec3f& b);
~~~

vector dot product

### Function dot()

~~~ .cpp
template <>
constexpr inline float dot(const vec4f& a, const vec4f& b);
~~~

vector dot product

### Function cross()

~~~ .cpp
template <>
constexpr inline float cross(const vec2f& a, const vec2f& b);
~~~

vector cross product (2d)

### Function cross()

~~~ .cpp
template <>
constexpr inline vec3f cross(const vec3f& a, const vec3f& b);
~~~

vector cross product (3d)

### Function length()

~~~ .cpp
template <typename T, int N>
constexpr inline T length(const vec<T, N>& a);
~~~

vector length

### Function lengthsqr()

~~~ .cpp
template <typename T, int N>
constexpr inline T lengthsqr(const vec<T, N>& a);
~~~

vector length squared

### Function normalize()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> normalize(const vec<T, N>& a);
~~~

vector normalization

### Function dist()

~~~ .cpp
template <typename T, int N>
constexpr inline T dist(const vec<T, N>& a, const vec<T, N>& b);
~~~

point distance

### Function distsqr()

~~~ .cpp
template <typename T, int N>
constexpr inline T distsqr(const vec<T, N>& a, const vec<T, N>& b);
~~~

point distance squared

### Function uangle()

~~~ .cpp
template <typename T, int N>
constexpr inline T uangle(const vec<T, N>& a, const vec<T, N>& b);
~~~

angle between normalized vectors

### Function angle()

~~~ .cpp
template <typename T, int N>
constexpr inline T angle(const vec<T, N>& a, const vec<T, N>& b);
~~~

angle between vectors

### Function lerp()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T t);
~~~

vector linear interpolation

### Function nlerp()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> nlerp(const vec<T, N>& a, const vec<T, N>& b, T t);
~~~

vector normalized linear interpolation

### Function slerp()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T t);
~~~

vector spherical linear interpolation (vectors have to be normalized)

### Function orthogonal()

~~~ .cpp
// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
template <typename T>
constexpr inline vec<T, 3> orthogonal(const vec<T, 3>& v);
~~~

orthogonal vector

### Function orthonormalize()

~~~ .cpp
template <typename T>
constexpr inline vec<T, 3> orthonormalize(
    const vec<T, 3>& a, const vec<T, 3>& b);
~~~

orthonormalize two vectors

### Function clamp()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> clamp(
    const vec<T, N>& x, const T& min, const T& max);
~~~

vector component-wise clamp

### Function clamp()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> clamp(
    const vec<T, N>& x, const vec<T, N>& min, const vec<T, N>& max);
~~~

vector component-wise clamp

### Function clamplen()

~~~ .cpp
template <typename T, int N, typename T1>
constexpr inline vec<T, N> clamplen(const vec<T, N> x, T1 max);
~~~

clamp the length of a vector

### Function min_element_idx()

~~~ .cpp
template <typename T, int N>
constexpr inline int min_element_idx(const vec<T, N>& a);
~~~

index of the min vector element

### Function max_element_idx()

~~~ .cpp
template <typename T, int N>
constexpr inline int max_element_idx(const vec<T, N>& a);
~~~

index of the max vector element

### Function min_element_val()

~~~ .cpp
template <typename T, int N>
constexpr inline T min_element_val(const vec<T, N>& a);
~~~

index of the min vector element

### Function max_element_val()

~~~ .cpp
template <typename T, int N>
constexpr inline T max_element_val(const vec<T, N>& a);
~~~

index of the max vector element

### Function pow()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> pow(const vec<T, N>& a, const T b);
~~~

Element-wise pow

### Function float_to_byte()

~~~ .cpp
template <int N>
constexpr inline vec<byte, N> float_to_byte(const vec<float, N>& a);
~~~

Element-wise conversion

### Function byte_to_float()

~~~ .cpp
template <int N>
constexpr inline vec<float, N> byte_to_float(const vec<byte, N>& a);
~~~

Element-wise conversion

### Struct mat

~~~ .cpp
template <typename T, int N, int M>
struct mat {
    using V = vec<T, N>;
    constexpr mat(); 
    constexpr mat(const std::initializer_list<V>& vv); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    V v[M];
}
~~~

Matrix of elements of compile time dimensions, stored in column major
format, with default initializer.
Colums access via operator[].

- Members:
    - V:      column data type
    - mat():      default constructor
    - mat():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - v[M]:      element data


### Struct mat <float, 2, 2 \>

~~~ .cpp
template <>
struct mat<float, 2, 2> {
    constexpr static const int N = 2, M = 2;
    using T = float;
    using V = vec<T, N>;
    constexpr mat(); 
    constexpr mat(const V& x, const V& y); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    V x;
    V y;
}
~~~

Specialization for 2x2 float matrices.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - mat():      default constructor
    - mat():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data


### Struct mat <float, 3, 3 \>

~~~ .cpp
template <>
struct mat<float, 3, 3> {
    constexpr static const int N = 3, M = 3;
    using T = float;
    using V = vec<T, N>;
    constexpr mat(); 
    constexpr mat(const V& x, const V& y, const V& z); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    V x;
    V y;
    V z;
}
~~~

Specialization for 3x3 float matrices.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - mat():      default constructor
    - mat():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


### Struct mat <float, 4, 4 \>

~~~ .cpp
template <>
struct mat<float, 4, 4> {
    constexpr static const int N = 4, M = 4;
    using T = float;
    using V = vec<T, N>;
    constexpr mat(); 
    constexpr mat(const V& x, const V& y, const V& z, const V& w); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    V x;
    V y;
    V z;
    V w;
}
~~~

Specialization for 4x4 float matrices.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - mat():      default constructor
    - mat():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


### Typedef mat1f

~~~ .cpp
using mat1f = mat<float, 1, 1>;
~~~

1-dimensional float matrix

### Typedef mat2f

~~~ .cpp
using mat2f = mat<float, 2, 2>;
~~~

2-dimensional float matrix

### Typedef mat3f

~~~ .cpp
using mat3f = mat<float, 3, 3>;
~~~

3-dimensional float matrix

### Typedef mat4f

~~~ .cpp
using mat4f = mat<float, 4, 4>;
~~~

4-dimensional float matrix

### Function identity_mat()

~~~ .cpp
template <typename T, int N>
constexpr inline mat<T, N, N> identity_mat();
~~~

Initialize an identity matrix.

### Function identity_mat()

~~~ .cpp
template <>
constexpr inline mat3f identity_mat();
~~~

Specialization for Initialize an identity matrix.

### Function identity_mat()

~~~ .cpp
template <>
constexpr inline mat4f identity_mat();
~~~

Specialization for Initialize an identity matrix.

### Constant identity_mat1f

~~~ .cpp
const auto identity_mat1f = identity_mat<float, 1>();
~~~

1-dimensional float identity matrix

### Constant identity_mat2f

~~~ .cpp
const auto identity_mat2f = identity_mat<float, 2>();
~~~

2-dimensional float identity matrix

### Constant identity_mat3f

~~~ .cpp
const auto identity_mat3f = identity_mat<float, 3>();
~~~

3-dimensional float identity matrix

### Constant identity_mat4f

~~~ .cpp
const auto identity_mat4f = identity_mat<float, 4>();
~~~

4-dimensional float identity matrix

### Function begin()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline vec<T, N>* begin(mat<T, N, M>& a);
~~~

iteration support

### Function begin()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline const vec<T, N>* begin(const mat<T, N, M>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline vec<T, N>* end(mat<T, N, M>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline const vec<T, N>* end(const mat<T, N, M>& a);
~~~

iteration support

### Function operator==()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline bool operator==(const mat<T, N, M>& a, const mat<T, N, M>& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline bool operator!=(const mat<T, N, M>& a, const mat<T, N, M>& b);
~~~

vector operator !=

### Function operator-()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator-(const mat<T, N, M>& a);
~~~

matrix operator -

### Function operator+()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator+(
    const mat<T, N, M>& a, const mat<T, N, M>& b);
~~~

matrix operator +

### Function operator*()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator*(const mat<T, N, M>& a, T b);
~~~

matrix scalar multiply

### Function operator/()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator/(const mat<T, N, M>& a, T b);
~~~

matrix scalar division

### Function operator*()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline vec<T, N> operator*(
    const mat<T, N, M>& a, const vec<T, M>& b);
~~~

matrix-vector right multiply

### Function operator*()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline vec<T, M> operator*(
    const vec<T, N>& a, const mat<T, N, M>& b);
~~~

matrix-vector left multiply

### Function operator*()

~~~ .cpp
template <typename T, int N, int M, int K>
constexpr inline mat<T, N, M> operator*(
    const mat<T, N, K>& a, const mat<T, K, M>& b);
~~~

matrix-matrix multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec2f operator*(const mat2f& a, const vec2f& b);
~~~

matrix-vector right multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec2f operator*(const vec2f& a, const mat2f& b);
~~~

matrix-vector left multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline mat2f operator*(const mat2f& a, const mat2f& b);
~~~

matrix-matrix multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec3f operator*(const mat3f& a, const vec3f& b);
~~~

matrix-vector right multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec3f operator*(const vec3f& a, const mat3f& b);
~~~

matrix-vector left multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline mat3f operator*(const mat3f& a, const mat3f& b);
~~~

matrix-matrix multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec4f operator*(const mat4f& a, const vec4f& b);
~~~

matrix-vector right multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline vec4f operator*(const vec4f& a, const mat4f& b);
~~~

matrix-vector left multiply

### Function operator*()

~~~ .cpp
template <>
constexpr inline mat4f operator*(const mat4f& a, const mat4f& b);
~~~

matrix-matrix multiply

### Function operator+=()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator+=(
    mat<T, N, M>& a, const mat<T, N, M>& b);
~~~

matrix sum assignment

### Function operator*=()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator*=(
    mat<T, N, M>& a, const mat<T, N, M>& b);
~~~

matrix-matrix multiply assignment

### Function operator*=()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator*=(mat<T, N, M>& a, const T& b);
~~~

matrix scaling assignment

### Function operator/=()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator/=(mat<T, N, M>& a, const T& b);
~~~

matrix scaling assignment

### Function mat_diagonal()

~~~ .cpp
template <typename T, int N>
constexpr vec<T, N> mat_diagonal(const mat<T, N, N>& a);
~~~

matrix diagonal

### Function transpose()

~~~ .cpp
template <typename T, int N, int M>
constexpr inline mat<T, M, N> transpose(const mat<T, N, M>& a);
~~~

matrix transpose

### Function adjugate()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a);
~~~

matrix adjugate (2x2)

### Function adjugate()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a);
~~~

matrix adjugate (3x3)

### Function adjugate()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a);
~~~

matrix adjugate (4x4)

### Function determinant()

~~~ .cpp
template <typename T>
constexpr inline T determinant(const mat<T, 2, 2>& a);
~~~

matrix determinant (2x2)

### Function determinant()

~~~ .cpp
template <typename T>
constexpr inline T determinant(const mat<T, 3, 3>& a);
~~~

matrix determinant (3x3)

### Function determinant()

~~~ .cpp
template <typename T>
constexpr inline T determinant(const mat<T, 4, 4>& a);
~~~

matrix determinant (4x4)

### Function inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline mat<T, N, N> inverse(const mat<T, N, N>& a);
~~~

matrix inverse (uses adjugate and determinant)

### Struct frame

~~~ .cpp
template <typename T, int N>
struct frame {
    using V = vec<T, N>;
    using M = mat<T, N, N>;
    constexpr frame(); 
    constexpr frame(const std::initializer_list<vec<T, N>>& vv); 
    constexpr frame(const M& m, const V& t); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    constexpr V& pos(); 
    constexpr const V& pos() const; 
    constexpr M& rot(); 
    constexpr const M& rot() const; 
    V v[N + 1];
}
~~~

Rigid transforms stored as a column-major affine matrix Nx(N+1).
In memory, this representation is equivalent to storing an NxN rotation
followed by a Nx1 translation. Viewed this way, the representation allows
also to retrive the axis of the coordinate frame as the first N column and
the translation as the N+1 column.
Colums access via operator[]. Access rotation and position with pos() and
rot().

- Members:
    - V:      column data type
    - M:      rotation data type
    - frame():      default constructor
    - frame():      element constructor
    - frame():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - pos():      access position
    - pos():      access position
    - rot():      access rotation
    - rot():      access rotation
    - 1]:      element data


### Struct frame <float, 2 \>

~~~ .cpp
template <>
struct frame<float, 2> {
    constexpr static const int N = 2;
    using T = float;
    using V = vec<T, N>;
    using M = mat<T, N, N>;
    constexpr frame(); 
    constexpr frame(const V& x, const V& y, const V& o); 
    constexpr frame(const M& m, const V& t); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    constexpr V& pos(); 
    constexpr const V& pos() const; 
    constexpr M& rot(); 
    constexpr const M& rot() const; 
    V x;
    V y;
    V o;
}
~~~

Specialization for 3D float frames.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - M:      rotation data type
    - frame():      default constructor
    - frame():      element constructor
    - frame():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - pos():      access position
    - pos():      access position
    - rot():      access rotation
    - rot():      access rotation
    - x:      element data
    - y:      element data
    - o:      element data


### Struct frame <float, 3 \>

~~~ .cpp
template <>
struct frame<float, 3> {
    constexpr static const int N = 3;
    using T = float;
    using V = vec<T, N>;
    using M = mat<T, N, N>;
    constexpr frame(); 
    constexpr frame(const V& x, const V& y, const V& z, const V& o); 
    constexpr frame(const M& m, const V& t); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    constexpr V* data(); 
    constexpr const V* data() const; 
    constexpr V& pos(); 
    constexpr const V& pos() const; 
    constexpr M& rot(); 
    constexpr const M& rot() const; 
    V x;
    V y;
    V z;
    V o;
}
~~~

Specialization for 3D float frames.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - M:      rotation data type
    - frame():      default constructor
    - frame():      element constructor
    - frame():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - pos():      access position
    - pos():      access position
    - rot():      access rotation
    - rot():      access rotation
    - x:      element data
    - y:      element data
    - z:      element data
    - o:      element data


### Typedef frame1f

~~~ .cpp
using frame1f = frame<float, 1>;
~~~

1-dimensional float frame

### Typedef frame2f

~~~ .cpp
using frame2f = frame<float, 2>;
~~~

2-dimensional float frame

### Typedef frame3f

~~~ .cpp
using frame3f = frame<float, 3>;
~~~

3-dimensional float frame

### Typedef frame4f

~~~ .cpp
using frame4f = frame<float, 4>;
~~~

4-dimensional float frame

### Function identity_frame()

~~~ .cpp
template <typename T, int N>
constexpr inline frame<T, N> identity_frame();
~~~

Initialize an identity frame.

### Function identity_frame()

~~~ .cpp
template <>
constexpr inline frame2f identity_frame();
~~~

Initialize an identity frame.

### Function identity_frame()

~~~ .cpp
template <>
constexpr inline frame3f identity_frame();
~~~

Initialize an identity frame.

### Constant identity_frame1f

~~~ .cpp
const auto identity_frame1f = identity_frame<float, 1>();
~~~

1-dimensional float identity frame

### Constant identity_frame2f

~~~ .cpp
const auto identity_frame2f = identity_frame<float, 2>();
~~~

2-dimensional float identity frame

### Constant identity_frame3f

~~~ .cpp
const auto identity_frame3f = identity_frame<float, 3>();
~~~

3-dimensional float identity frame

### Constant identity_frame4f

~~~ .cpp
const auto identity_frame4f = identity_frame<float, 4>();
~~~

4-dimensional float identity frame

### Function pos()

~~~ .cpp
template <typename T, int N>
constexpr inline const vec<T, N>& pos(const frame<T, N>& f);
~~~

frame position const access

### Function rot()

~~~ .cpp
template <typename T, int N>
constexpr inline const mat<T, N, N>& rot(const frame<T, N>& f);
~~~

frame rotation const access

### Function pos()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>& pos(frame<T, N>& f);
~~~

frame position reference

### Function rot()

~~~ .cpp
template <typename T, int N>
constexpr inline mat<T, N, N>& rot(frame<T, N>& f);
~~~

frame rotation reference

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>* begin(frame<T, N>& a);
~~~

iteration support

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline const vec<T, N>* begin(const frame<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>* end(frame<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline const vec<T, N>* end(const frame<T, N>& a);
~~~

iteration support

### Function to_mat()

~~~ .cpp
template <typename T, int N>
constexpr inline mat<T, N + 1, N + 1> to_mat(const frame<T, N>& a);
~~~

frame to matrix conversion

### Function to_frame()

~~~ .cpp
template <typename T, int N>
constexpr inline frame<T, N - 1> to_frame(const mat<T, N, N>& a);
~~~

matrix to frame conversion

### Function operator==()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator==(const frame<T, N>& a, const frame<T, N>& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b);
~~~

vector operator !=

### Function operator*()

~~~ .cpp
template <typename T, int N>
constexpr inline frame<T, N> operator*(
    const frame<T, N>& a, const frame<T, N>& b);
~~~

frame composition (equivalent to affine matrix multiply)

### Function inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline frame<T, N> inverse(const frame<T, N>& a);
~~~

frame inverse (equivalent to rigid affine inverse)

### Function operator*()

~~~ .cpp
template <>
constexpr inline frame3f operator*(const frame3f& a, const frame3f& b);
~~~

frame composition (equivalent to affine matrix multiply)

### Function inverse()

~~~ .cpp
template <>
constexpr inline frame3f inverse(const frame3f& a);
~~~

frame inverse (equivalent to rigid affine inverse)

### Struct quat

~~~ .cpp
template <typename T, int N>
struct quat;
~~~

Quaternion placeholder. Only helpful in the specialization.

### Struct quat <T, 4 \>

~~~ .cpp
template <typename T>
struct quat<T, 4> {
    constexpr static const int N = 4;
    constexpr quat(); 
    constexpr explicit quat(const vec<T, N>& vv); 
    constexpr explicit operator vec<T, N>() const; 
    constexpr T& operator[](int i); 
    constexpr const T& operator[](int i) const; 
    constexpr T* data(); 
    constexpr const T* data() const; 
    T x;
    T y;
    T z;
    T w;
}
~~~

Quaternions implemented as a vec<T,4>. Data access via operator[].
Quaterions are xi + yj + zk + w.

- Members:
    - N:      size
    - quat():      default constructor
    - quat():      conversion from vec
    - operator N>():      conversion to vec
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      data
    - y:      data
    - z:      data
    - w:      data


### Typedef quat4f

~~~ .cpp
using quat4f = quat<float, 4>;
~~~

float quaterion

### Constant identity_quat4f

~~~ .cpp
const auto identity_quat4f = quat<float, 4>{0, 0, 0, 1};
~~~

float identity quaterion

### Function operator==()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator==(const quat<T, N>& a, const quat<T, N>& b);
~~~

vector operator ==

### Function operator!=()

~~~ .cpp
template <typename T, int N>
constexpr inline bool operator!=(const quat<T, N>& a, const quat<T, N>& b);
~~~

vector operator !=

### Function operator*()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b);
~~~

quaterion multiply

### Function conjugate()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> conjugate(const quat<T, 4>& v);
~~~

quaterion conjugate

### Function inverse()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> inverse(const quat<T, 4>& v);
~~~

quaterion inverse

### Function normalize()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> normalize(const quat<T, 4>& v);
~~~

quaterion inverse

### Function nlerp()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T t);
~~~

quaterion normalized linear interpolation

### Function slerp()

~~~ .cpp
template <typename T>
constexpr quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T t);
~~~

quaterion spherical linear interpolation

### Struct bbox

~~~ .cpp
template <typename T, int N>
struct bbox {
    using V = vec<T, N>;
    constexpr bbox(); 
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    V min;
    V max;
}
~~~

Axis aligned bounding box represented as a min/max vector pair.
Access min/max with operator[].

- Members:
    - V:      column data type
    - bbox():      initializes an invalid bbox
    - bbox():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


### Struct bbox <float, 1 \>

~~~ .cpp
template <>
struct bbox<float, 1> {
    constexpr static const int N = 1;
    using T = float;
    using V = vec<T, N>;
    constexpr bbox(); 
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    V min;
    V max;
}
~~~

Specialization for float 3D bounding boxes.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - bbox():      initializes an invalid bbox
    - bbox():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


### Struct bbox <float, 2 \>

~~~ .cpp
template <>
struct bbox<float, 2> {
    constexpr static const int N = 2;
    using T = float;
    using V = vec<T, N>;
    constexpr bbox(); 
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    V min;
    V max;
}
~~~

Specialization for float 3D bounding boxes.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - bbox():      initializes an invalid bbox
    - bbox():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


### Struct bbox <float, 3 \>

~~~ .cpp
template <>
struct bbox<float, 3> {
    constexpr static const int N = 3;
    using T = float;
    using V = vec<T, N>;
    constexpr bbox(); 
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    V min;
    V max;
}
~~~

Specialization for float 3D bounding boxes.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - bbox():      initializes an invalid bbox
    - bbox():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


### Struct bbox <float, 4 \>

~~~ .cpp
template <>
struct bbox<float, 4> {
    constexpr static const int N = 4;
    using T = float;
    using V = vec<T, N>;
    constexpr bbox(); 
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M); 
    constexpr V& operator[](int i); 
    constexpr const V& operator[](int i) const; 
    V min;
    V max;
}
~~~

Specialization for float 3D bounding boxes.

- Members:
    - N:      size
    - T:      type
    - V:      column data type
    - bbox():      initializes an invalid bbox
    - bbox():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


### Typedef bbox1f

~~~ .cpp
using bbox1f = bbox<float, 1>;
~~~

1-dimensional float bbox

### Typedef bbox2f

~~~ .cpp
using bbox2f = bbox<float, 2>;
~~~

2-dimensional float bbox

### Typedef bbox3f

~~~ .cpp
using bbox3f = bbox<float, 3>;
~~~

3-dimensional float bbox

### Typedef bbox4f

~~~ .cpp
using bbox4f = bbox<float, 4>;
~~~

4-dimensional float bbox

### Function invalid_bbox()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> invalid_bbox();
~~~

initializes an empty bbox

### Function make_bbox()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> make_bbox(int count, const vec<T, N>* v);
~~~

initialize a bonding box from a list of points

### Function make_bbox()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> make_bbox(
    const std::initializer_list<vec<T, N>>& v);
~~~

initialize a bonding box from a list of points

### Constant invalid_bbox1f

~~~ .cpp
const auto invalid_bbox1f = bbox1f();
~~~

1-dimensional float empty bbox

### Constant invalid_bbox2f

~~~ .cpp
const auto invalid_bbox2f = bbox2f();
~~~

2-dimensional float empty bbox

### Constant invalid_bbox3f

~~~ .cpp
const auto invalid_bbox3f = bbox3f();
~~~

3-dimensional float empty bbox

### Constant invalid_bbox4f

~~~ .cpp
const auto invalid_bbox4f = bbox4f();
~~~

4-dimensional float empty bbox

### Function center()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> center(const bbox<T, N>& a);
~~~

computes the center of a bbox

### Function diagonal()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> diagonal(const bbox<T, N>& a);
~~~

computes the diagonal of a bbox

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>* begin(bbox<T, N>& a);
~~~

iteration support

### Function begin()

~~~ .cpp
template <typename T, int N>
constexpr inline const vec<T, N>* begin(const bbox<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N>* end(bbox<T, N>& a);
~~~

iteration support

### Function end()

~~~ .cpp
template <typename T, int N>
constexpr inline const vec<T, N>* end(const bbox<T, N>& a);
~~~

iteration support

### Function expand()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> expand(const bbox<T, N>& a, const vec<T, N>& b);
~~~

expands a bounding box with a point

### Function expand()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> expand(const bbox<T, N>& a, const bbox<T, N>& b);
~~~

expands a bounding box with a bounding box

### Function contains()

~~~ .cpp
template <typename T, int N>
constexpr inline bool contains(const bbox<T, N>& a, const vec<T, N>& b);
~~~

check if a bounding box contains a point

### Function contains()

~~~ .cpp
template <typename T, int N>
constexpr inline bool contains(const bbox<T, N>& a, const bbox<T, N>& b);
~~~

check if a bounding box contains a bounding box

### Function expand()

~~~ .cpp
template <>
constexpr inline bbox3f expand(const bbox3f& a, const vec3f& b);
~~~

expands a bounding box with a point

### Function expand()

~~~ .cpp
template <>
constexpr inline bbox3f expand(const bbox3f& a, const bbox3f& b);
~~~

expands a bounding box with a bounding box

### Function contains()

~~~ .cpp
template <>
constexpr inline bool contains(const bbox3f& a, const vec3f& b);
~~~

check if a bounding box contains a point

### Function contains()

~~~ .cpp
template <>
constexpr inline bool contains(const bbox3f& a, const bbox3f& b);
~~~

check if a bounding box contains a bounding box

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> operator+(const bbox<T, N>& a, const T& b);
~~~

same as expand()

### Function operator+()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N> operator+(
    const bbox<T, N>& a, const bbox<T, N>& b);
~~~

same as expand()

### Function operator+=()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const vec<T, N>& b);
~~~

assign to expand()

### Function operator+=()

~~~ .cpp
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const bbox<T, N>& b);
~~~

assign to expand()

### Struct ray

~~~ .cpp
template <typename T, int N>
struct ray {
    vec<T, N> o;
    vec<T, N> d;
    T tmin;
    T tmax;
    constexpr ray(); 
    constexpr ray(const vec<T, N>& o, const vec<T, N>& d, T tmin = 0, T tmax = std::numeric_limits<T>::max()); 
}
~~~

Rays with origin, direction and min/max t value.

- Members:
    - o:      origin
    - d:      direction
    - tmin:      minimum distance
    - tmax:      maximum distance
    - ray():      default constructor
    - ray():      initializes a ray from its elements


### Struct ray <float, 3 \>

~~~ .cpp
template <>
struct ray<float, 3> {
    constexpr static const int N = 3;
    using T = float;
    vec<T, N> o;
    vec<T, N> d;
    T tmin;
    T tmax;
    constexpr ray(); 
    constexpr ray( const vec<T, N>& o, const vec<T, N>& d, T tmin = 0, T tmax = flt_max); 
}
~~~

Sepcialization for 3D float rays.

- Members:
    - N:      size
    - T:      type
    - o:      origin
    - d:      direction
    - tmin:      minimum distance
    - tmax:      maximum distance
    - ray():      default constructor
    - ray():      initializes a ray from its elements


### Typedef ray1f

~~~ .cpp
using ray1f = ray<float, 1>;
~~~

1-dimensional float ray

### Typedef ray2f

~~~ .cpp
using ray2f = ray<float, 2>;
~~~

2-dimensional float ray

### Typedef ray3f

~~~ .cpp
using ray3f = ray<float, 3>;
~~~

3-dimensional float ray

### Typedef ray4f

~~~ .cpp
using ray4f = ray<float, 4>;
~~~

4-dimensional float ray

### Function eval()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> eval(const ray<T, N>& ray, T t);
~~~

evalutes the position along the ray

### Function transform_point()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b);
~~~

transforms a point by a matrix

### Function transform_vector()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b);
~~~

transforms a vector by a matrix

### Function transform_direction()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b);
~~~

transforms a direction by a matrix

### Function transform_point()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

transforms a point by a frame (rigid affine transform)

### Function transform_vector()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

transforms a vector by a frame (rigid affine transform)

### Function transform_direction()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

transforms a direction by a frame (rigid affine transform)

### Function transform_frame()

~~~ .cpp
template <typename T, int N>
constexpr inline frame<T, N> transform_frame(
    const frame<T, N>& a, const frame<T, N>& b);
~~~

transforms a frame by a frame (rigid affine transform)

### Function transform_point_inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

inverse transforms a point by a frame (rigid affine transform)

### Function transform_vector_inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

inverse transforms a vector by a frame (rigid affine transform)

### Function transform_direction_inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b);
~~~

inverse transforms a direction by a frame (rigid affine transform)

### Function transform_point()

~~~ .cpp
template <>
constexpr inline vec3f transform_point(const mat4f& a, const vec3f& b);
~~~

transforms a point by a matrix

### Function transform_vector()

~~~ .cpp
template <>
constexpr inline vec3f transform_vector(const mat4f& a, const vec3f& b);
~~~

transforms a vector by a matrix

### Function transform_direction()

~~~ .cpp
template <>
constexpr inline vec3f transform_direction(const mat4f& a, const vec3f& b);
~~~

transforms a direction by a matrix

### Function transform_point()

~~~ .cpp
template <>
constexpr inline vec3f transform_point(const frame3f& a, const vec3f& b);
~~~

transforms a point by a frame (rigid affine transform)

### Function transform_vector()

~~~ .cpp
template <>
constexpr inline vec3f transform_vector(const frame3f& a, const vec3f& b);
~~~

transforms a vector by a frame (rigid affine transform)

### Function transform_direction()

~~~ .cpp
template <>
constexpr inline vec3f transform_direction(const frame3f& a, const vec3f& b);
~~~

transforms a direction by a frame (rigid affine transform)

### Function transform_frame()

~~~ .cpp
template <>
constexpr inline frame3f transform_frame(const frame3f& a, const frame3f& b);
~~~

transforms a frame by a frame (rigid affine transform)

### Function transform_point_inverse()

~~~ .cpp
template <>
constexpr inline vec3f transform_point_inverse(
    const frame3f& a, const vec3f& b);
~~~

inverse transforms a point by a frame (rigid affine transform)

### Function transform_vector_inverse()

~~~ .cpp
template <>
constexpr inline vec3f transform_vector_inverse(
    const frame3f& a, const vec3f& b);
~~~

inverse transforms a vector by a frame (rigid affine transform)

### Function transform_direction_inverse()

~~~ .cpp
template <>
constexpr inline vec3f transform_direction_inverse(
    const frame3f& a, const vec3f& b);
~~~

inverse transforms a direction by a frame (rigid affine transform)

### Function transform_ray()

~~~ .cpp
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const mat<T, N + 1, N + 1>& a, const ray<T, N>& b);
~~~

transforms a ray by a matrix

### Function transform_bbox()

~~~ .cpp
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const mat<T, 4, 4>& a, const bbox<T, 3>& b);
~~~

transforms a bbox by a matrix

### Function transform_ray()

~~~ .cpp
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const frame<T, N>& a, const ray<T, N>& b);
~~~

transforms a ray by a frame (rigid affine transform)

### Function transform_bbox()

~~~ .cpp
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const frame<T, 3>& a, const bbox<T, 3>& b);
~~~

transforms a bbox by a frame (rigid affine transform)

### Function transform_ray_inverse()

~~~ .cpp
template <typename T, int N>
constexpr inline ray<T, N> transform_ray_inverse(
    const frame<T, N>& a, const ray<T, N>& b);
~~~

inverse transforms a ray by a frame (rigid affine transform)

### Function transform_bbox_inverse()

~~~ .cpp
template <typename T>
constexpr inline bbox<T, 3> transform_bbox_inverse(
    const frame<T, 3>& a, const bbox<T, 3>& b);
~~~

inverse transforms a bbox by a frame (rigid affine transform)

### Function rotation_mat3()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 3, 3> rotation_mat3(const vec<T, 3>& axis, T angle);
~~~

rotation matrix from axis-angle

### Function translation_frame3()

~~~ .cpp
template <typename T>
constexpr inline frame<T, 3> translation_frame3(const vec<T, 3>& a);
~~~

translation frame

### Function translation_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> translation_mat4(const vec<T, 3>& a);
~~~

translation matrix

### Function scaling_frame3()

~~~ .cpp
template <typename T>
constexpr inline frame<T, 3> scaling_frame3(const vec<T, 3>& a);
~~~

scaling frame (this is not rigid and it is only here for symmatry of
API)

### Function scaling_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> scaling_mat4(const vec<T, 3>& a);
~~~

scaling matrix

### Function rotation_frame3()

~~~ .cpp
template <typename T>
constexpr inline frame<T, 3> rotation_frame3(const vec<T, 3>& axis, T angle);
~~~

rotation frame

### Function rotation_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const mat<T, 3, 3>& rot);
~~~

rotation matrix

### Function rotation_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const vec<T, 3>& axis, T angle);
~~~

rotation matrix

### Function rotation_axisangle4()

~~~ .cpp
template <typename T>
constexpr inline vec<T, 4> rotation_axisangle4(const quat<T, 4>& a);
~~~

quaternion axis-angle conversion

### Function rotation_quat4()

~~~ .cpp
template <typename T>
constexpr inline quat<T, 4> rotation_quat4(const vec<T, 4>& axis_angle);
~~~

axis-angle to quaternion

### Function rotation_mat3()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 3, 3> rotation_mat3(const quat<T, 4>& v);
~~~

quaterion to matrix conversion

### Function rotation_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const quat<T, 4>& v);
~~~

rotation matrix

### Function rotation_quat4()

~~~ .cpp
template <typename T>
constexpr inline quat<T, 4> rotation_quat4(const mat<T, 3, 3>& m_);
~~~

matrix to quaternion

### Function lookat_frame3()

~~~ .cpp
template <typename T>
constexpr inline frame<T, 3> lookat_frame3(
    const vec<T, 3>& eye, const vec<T, 3>& center, const vec<T, 3>& up);
~~~

OpenGL lookat frame

### Function lookat_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> lookat_mat4(
    const vec<T, 3>& eye, const vec<T, 3>& center, const vec<T, 3>& up);
~~~

OpenGL lookat matrix

### Function frustum_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> frustum_mat4(T l, T r, T b, T t, T n, T f);
~~~

OpenGL frustum matrix

### Function ortho_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> ortho_mat4(T l, T r, T b, T t, T n, T f);
~~~

OpenGL orthographic matrix

### Function ortho2d_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> ortho2d_mat4(T left, T right, T bottom, T top);
~~~

OpenGL orthographic 2D matrix

### Function ortho_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> ortho_mat4(T xmag, T ymag, T near, T far);
~~~

OpenGL/GLTF orthographic matrix

### Function perspective_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat4(
    T fovy, T aspect, T near, T far);
~~~

OpenGL/GLTF perspective matrix

### Function perspective_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat4(T fovy, T aspect, T near);
~~~

OpenGL/GLTF infinite perspective matrix

### Function decompose_mat4()

~~~ .cpp
template <typename T>
constexpr inline void decompose_mat4(const mat<T, 4, 4>& m,
    vec<T, 3>& translation, mat<T, 3, 3>& rotation, vec<T, 3>& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

### Function decompose_mat4()

~~~ .cpp
template <typename T>
constexpr inline void decompose_mat4(const mat<T, 4, 4>& m,
    vec<T, 3>& translation, quat<T, 4>& rotation, vec<T, 3>& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

### Function compose_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> compose_mat4(const vec<T, 3>& translation,
    const mat<T, 3, 3>& rotation, const vec<T, 3>& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

### Function compose_mat4()

~~~ .cpp
template <typename T>
constexpr inline mat<T, 4, 4> compose_mat4(const vec<T, 3>& translation,
    const quat<T, 4>& rotation, const vec<T, 3>& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

### Struct rng_pcg32

~~~ .cpp
struct rng_pcg32 {
~~~

PCG random numbers. A family of random number generators that supports
multiple sequences. In our code, we allocate one sequence for each sample.
PCG32 from http://www.pcg-random.org/

### Function next()

~~~ .cpp
constexpr inline uint32_t next(rng_pcg32* rng);
~~~

Next random number

### Function init()

~~~ .cpp
constexpr inline void init(rng_pcg32* rng, uint64_t state, uint64_t seq);
~~~

Init a random number generator with a state state from the sequence seq.

### Function next1f()

~~~ .cpp
inline float next1f(rng_pcg32* rng);
~~~

Next random float in [0,1).

### Function next2f()

~~~ .cpp
inline vec2f next2f(rng_pcg32* rng);
~~~

Next random float in [0,1)x[0,1).

### Function hash_permute()

~~~ .cpp
constexpr inline uint32_t hash_permute(uint32_t i, uint32_t n, uint32_t key);
~~~

Computes the i-th term of a permutation of l values keyed by p.
From Correlated Multi-Jittered Sampling by Kensler @ Pixar

### Function hash_randfloat()

~~~ .cpp
constexpr inline float hash_randfloat(uint32_t i, uint32_t key);
~~~

Computes a float value by hashing i with a key p.
From Correlated Multi-Jittered Sampling by Kensler @ Pixar

### Function hash_uint64()

~~~ .cpp
constexpr inline uint64_t hash_uint64(uint64_t a);
~~~

64 bit integer hash. Public domain code.

### Function hash_uint64_32()

~~~ .cpp
constexpr inline uint32_t hash_uint64_32(uint64_t a);
~~~

64-to-32 bit integer hash. Public domain code.

### Function hash_combine()

~~~ .cpp
constexpr inline int hash_combine(int a, int b);
~~~

Combines two 64 bit hashes as in boost::hash_combine

### Function hash_vec()

~~~ .cpp
template <typename T, int N>
constexpr inline int hash_vec(const vec<T, N>& v);
~~~

Hash a vector with hash_combine() and std::hash

### Function triangle_normal()

~~~ .cpp
template <typename T>
constexpr inline vec<T, 3> triangle_normal(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2);
~~~

triangle normal

### Function triangle_area()

~~~ .cpp
template <typename T>
constexpr inline T triangle_area(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2);
~~~

triangle area

### Function tetrahedron_volume()

~~~ .cpp
template <typename T>
constexpr inline T tetrahedron_volume(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3);
~~~

tetrahedron volume

### Function blerp()

~~~ .cpp
template <typename T, typename T1>
constexpr inline T blerp(const T& a, const T& b, const T& c, const T1& w);
~~~

triangle baricentric interpolation

### Function compute_tangents()

~~~ .cpp
inline void compute_tangents(int nlines, const vec2i* lines, int nverts,
    const vec3f* pos, vec3f* tang, bool weighted = true);
~~~

Compute smoothed tangents (for lines).

Parameters:
- nverts/pos: array pf vertex positions
- npoints/points: array of point indices
- nlines/lines: array of point indices
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- tang: preallocated array of computed normals

### Function compute_tangents()

~~~ .cpp
inline void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    bool weighted = true);
~~~

Compute smoothed tangents.

Parameters:
- nverts/pos: array pf vertex positions
- npoints/points: array of point indices
- nlines/lines: array of point indices
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- tang: array of computed tangents

### Function compute_normals()

~~~ .cpp
inline void compute_normals(int ntriangles, const vec3i* triangles, int nverts,
    const vec3f* pos, vec3f* norm, bool weighted = true);
~~~

Compute smoothed normals.

Parameters:
- nverts/pos: array pf vertex positions
- npoints/points: array of point indices
- nlines/lines: array of point indices
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- norm: preallocated array of computed normals

### Function compute_normals()

~~~ .cpp
inline void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true);
~~~

Compute smoothed normals.

Parameters:
- nverts/pos: array pf vertex positions
- npoints/points: array of point indices
- nlines/lines: array of point indices
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- norm: array of computed normals

### Function compute_tangent_frame()

~~~ .cpp
inline void compute_tangent_frame(int ntriangles, const vec3i* triangles,
    int nverts, const vec3f* pos, const vec3f* norm, const vec2f* texcoord,
    vec4f* tangsp, bool weighted = true);
~~~

Compute tangent frame for triangle mesh. Tangent space is defined by
a four component vector. The first three components are the tangent
with respect to the U texcoord. The fourth component is the sign of the
tangent wrt the V texcoord. Tangent frame is useful in normal mapping.

Parameters:
- nverts/pos: array pf vertex positions
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- tangsp: preallocated array of computed tangent space

### Function compute_tangent_frame()

~~~ .cpp
inline void compute_tangent_frame(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted = true);
~~~

Compute tangent frame for triangle mesh. Tangent space is defined by
a four component vector. The first three components are the tangent
with respect to the U texcoord. The fourth component is the sign of the
tangent wrt the V texcoord. Tangent frame is useful in normal mapping.

Parameters:
- nverts/pos: array pf vertex positions
- ntriangles/triangles: array of point indices
- weighted: whether to use area weighting (typically true)

Out Parameters:
- tangsp: array of computed tangent space

### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(int nverts, const vec3f* pos, const vec3f* norm,
    const vec4f* weights, const vec4i* joints, const mat4f* xforms,
    vec3f* skinned_pos, vec3f* skinned_norm);
~~~

Apply skinning

### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(int nverts, const vec3f* pos, const vec3f* norm,
    const vec4f* weights, const vec4i* joints, const frame3f* xforms,
    vec3f* skinned_pos, vec3f* skinned_norm);
~~~

Apply skinning

### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
~~~

Apply skinning

### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
~~~

Apply skinning

### Function compute_matrix_skinning()

~~~ .cpp
inline void compute_matrix_skinning(int nverts, const vec3f* pos,
    const vec3f* norm, const vec4f* weights, const vec4i* joints,
    const mat4f* xforms, vec3f* skinned_pos, vec3f* skinned_norm);
~~~

Apply skinning as specified in Khronos glTF

### Function compute_matrix_skinning()

~~~ .cpp
inline void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
~~~

Apply skinning as specified in Khronos glTF

### Function make_triangles()

~~~ .cpp
template <typename PosFunc, typename NormFunc, typename TexcoordFunc>
inline void make_triangles(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const PosFunc& pos_fn, const NormFunc& norm_fn,
    const TexcoordFunc& texcoord_fn);
~~~

Generate a parametric surface with callbacks.

Parameters:
- usteps: subdivisions in u
- vsteps: subdivisions in v
- pos_fn: pos callbacks (vec2f -> vec3f)
- norm_fn: norm callbacks (vec2f -> vec3f)
- texcoord_fn: texcoord callbacks (vec2f -> vec2f)

Out Parameters:
- triangles: element array
- pos/norm/texcoord: vertex position/normal/texcoords

### Function make_lines()

~~~ .cpp
template <typename PosFunc, typename TangFunc, typename TexcoordFunc,
    typename RadiusFunc>
inline void make_lines(int num, int usteps, std::vector<vec2i>& lines,
    std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    std::vector<vec2f>& texcoord, std::vector<float>& radius,
    const PosFunc& pos_fn, const TangFunc& tang_fn,
    const TexcoordFunc& texcoord_fn, const RadiusFunc& radius_fn);
~~~

Generate parametric lines with callbacks.

Parameters:
- usteps: subdivisions in u
- num: number of lines
- pos_fn: pos callbacks ((int, float) -> vec3f)
- tang_fn: tangent callbacks ((int, float) -> vec3f)
- texcoord_fn: texcoord callbacks ((int, float) -> vec2f)
- radius_fn: radius callbacks ((int, float) -> float)

Out Parameters:
- lines: element array
- pos/tang/texcoord/radius: vertex position/tangent/texcoords/radius

### Function make_points()

~~~ .cpp
template <typename PosFunc, typename NormFunc, typename TexcoordFunc,
    typename RadiusFunc>
inline void make_points(int num, std::vector<int>& points,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<float>& radius,
    const PosFunc& pos_fn, const NormFunc& norm_fn,
    const TexcoordFunc& texcoord_fn, const RadiusFunc& radius_fn);
~~~

Generate a parametric point set. Mostly here for completeness.

Parameters:
- num: number of points
- pos_fn: pos callbacks (int -> vec3f)
- norm_fn: norm callbacks (int -> vec3f)
- texcoord_fn: texcoord callbacks (int -> vec2f)
- radius_fn: radius callbacks (int -> float)

Out Parameters:
- points: element array
- pos/norm/texcoord/radius: vertex position/normal/texcoords/radius

### Function merge_triangles()

~~~ .cpp
inline void merge_triangles(std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, const std::vector<vec3i>& mtriangles,
    const std::vector<vec3f>& mpos, const std::vector<vec3f>& mnorm,
    const std::vector<vec2f>& mtexcoord);
~~~

Merge a triangle mesh into another.

### Function sample_lines_cdf()

~~~ .cpp
inline void sample_lines_cdf(
    int nlines, const vec2i* lines, const vec3f* pos, float* cdf);
~~~

Compute a distribution for sampling lines uniformly

### Function sample_lines_cdf()

~~~ .cpp
inline std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos);
~~~

Compute a distribution for sampling lines uniformly

### Function sample_lines()

~~~ .cpp
inline std::pair<int, vec2f> sample_lines(
    int nlines, const float* cdf, float re, float ruv);
~~~

Pick a point on lines

### Function sample_lines()

~~~ .cpp
inline std::pair<int, vec2f> sample_lines(
    const std::vector<float>& cdf, float re, float ruv);
~~~

Pick a point on lines

### Function sample_triangles_cdf()

~~~ .cpp
inline void sample_triangles_cdf(
    int ntriangles, const vec3i* triangles, const vec3f* pos, float* cdf);
~~~

Compute a distribution for sampling triangle meshes uniformly

### Function sample_triangles_cdf()

~~~ .cpp
inline std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos);
~~~

Pick a point on lines

### Function sample_triangles()

~~~ .cpp
inline std::pair<int, vec3f> sample_triangles(
    int ntriangles, const float* cdf, float re, const vec2f& ruv);
~~~

Pick a point on a triangle mesh

### Function sample_triangles()

~~~ .cpp
inline std::pair<int, vec3f> sample_triangles(
    const std::vector<float>& cdf, float re, const vec2f& ruv);
~~~

Pick a point on a triangle mesh

### Function sample_triangles_points()

~~~ .cpp
inline void sample_triangles_points(int ntriangles, const vec3i* triangles,
    const vec3f* pos, const vec3f* norm, const vec2f* texcoord, int npoints,
    vec3f* sampled_pos, vec3f* sampled_norm, vec2f* sampled_texcoord,
    uint64_t seed);
~~~

Samples a set of points over a triangle mesh uniformly. The rng function
takes the point index and returns vec3f numbers uniform directibuted in
[0,1]^3. ù norm and texcoord are optional.

### Function sample_triangles_points()

~~~ .cpp
inline void sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints,
    std::vector<vec3f>& sampled_pos, std::vector<vec3f>& sampled_norm,
    std::vector<vec2f>& sampled_texcoord, uint64_t seed);
~~~

Samples a set of points over a triangle mesh uniformly.
Wrapper to the above function.

### Function make_uvsphere()

~~~ .cpp
inline void make_uvsphere(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord);
~~~

Make a sphere.

### Function make_uvhemisphere()

~~~ .cpp
inline void make_uvhemisphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a sphere.

### Function make_uvflippedsphere()

~~~ .cpp
inline void make_uvflippedsphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make an inside-out sphere.

### Function make_uvflippedhemisphere()

~~~ .cpp
inline void make_uvflippedhemisphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make an inside-out hemisphere

### Function make_uvquad()

~~~ .cpp
inline void make_uvquad(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvcube()

~~~ .cpp
inline void make_uvcube(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvspherecube()

~~~ .cpp
inline void make_uvspherecube(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvspherizedcube()

~~~ .cpp
inline void make_uvspherizedcube(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvflipcapsphere()

~~~ .cpp
inline void make_uvflipcapsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvcutsphere()

~~~ .cpp
inline void make_uvcutsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function make_uvflippedcutsphere()

~~~ .cpp
inline void make_uvflippedcutsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord);
~~~

Make a quad.

### Function intersect_point()

~~~ .cpp
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& ray_t);
~~~

Intersect a ray with a point (approximate)

Parameters:
- ray: ray origin and direction, parameter min, max range
- p: point position
- r: point radius

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: primitive uv ( {0,0} for points )

Returns:
- whether the intersection occurred

Iplementation Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm finds the closest point on the ray segment to the point and
   test their distance with the point radius
- based on http://geomalgorithms.com/a02-lines.html.

### Function intersect_line()

~~~ .cpp
inline bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv);
~~~

Intersect a ray with a line

Parameters:
- ray: ray origin and direction, parameter min, max range
- v0, v1: line segment points
- r0, r1: line segment radia

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: euv[0] is the line parameter at the intersection ( euv[1] is zero )

Returns:
- whether the intersection occurred

Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm find the closest points on line and ray segment and test
  their distance with the line radius at that location
- based on http://geomalgorithms.com/a05-intersect-1.html
- based on http://geomalgorithms.com/a07-distance.html#
    dist3D_Segment_to_Segment

### Function intersect_triangle()

~~~ .cpp
inline bool intersect_triangle(const ray3f& ray, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float& ray_t, vec3f& euv);
~~~

Intersect a ray with a triangle

Parameters:
- ray: ray origin and direction, parameter min, max range
- v0, v1, v2: triangle vertices

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: baricentric coordinates of the intersection

Returns:
- whether the intersection occurred

Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm based on Muller-Trombone intersection test

### Function intersect_tetrahedron()

~~~ .cpp
inline bool intersect_tetrahedron(const ray3f& ray_, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float& ray_t,
    vec4f& euv);
~~~

Intersect a ray with a tetrahedron. Note that we consider only intersection
wiht the tetrahedra surface and discount intersction with the interior.

Parameters:
- ray: ray to intersect with
- v0, v1, v2: triangle vertices

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: baricentric coordinates of the intersection

Returns:
- whether the intersection occurred

TODO: check order
TODO: uv

### Function intersect_check_bbox()

~~~ .cpp
inline bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox);
~~~

Intersect a ray with a axis-aligned bounding box

Parameters:
- ray: ray to intersect with
- bbox: bounding box min/max bounds

Returns:
- whether the intersection occurred

### Function _safemin()

~~~ .cpp
template <typename T>
static inline const T& _safemin(const T& a, const T& b);
~~~

Min/max used in BVH traversal. Copied here since the traversal code relies
on the specific behaviour wrt NaNs.

### Function _safemax()

~~~ .cpp
template <typename T>
static inline const T& _safemax(const T& a, const T& b);
~~~

Min/max used in BVH traversal. Copied here since the traversal code relies
on the specific behaviour wrt NaNs.

### Function intersect_check_bbox()

~~~ .cpp
inline bool intersect_check_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);
~~~

Intersect a ray with a axis-aligned bounding box

Parameters:
- ray_o, ray_d: ray origin and direction
- ray_tmin, ray_tmax: ray parameter min, max range
- ray_dinv: ray inverse direction
- ray_dsign: ray direction sign
- bbox_min, bbox_max: bounding box min/max bounds

Returns:
- whether the intersection occurred

Implementation Notes:
- based on "Robust BVH Ray Traversal" by T. Ize published at
http://jcgt.org/published/0002/02/02/paper.pdf

### Function point_bbox()

~~~ .cpp
inline bbox3f point_bbox(const vec3f& p, float r = 0);
~~~

Point bounds

### Function line_bbox()

~~~ .cpp
inline bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0);
~~~

Line bounds

### Function triangle_bbox()

~~~ .cpp
inline bbox3f triangle_bbox(const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

Triangle bounds

### Function tetrahedron_bbox()

~~~ .cpp
inline bbox3f tetrahedron_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3);
~~~

Tetrahedron bounds

### Function turntable()

~~~ .cpp
template <typename T>
constexpr inline void turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);
~~~

Turntable for UI navigation from a from/to/up parametrization of the camera.

### Function turntable()

~~~ .cpp
template <typename T>
constexpr inline void turntable(frame<T, 3>& frame, float& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);
~~~

Turntable for UI navigation for a frame/distance parametrization of the
camera.

### Struct image

~~~ .cpp
template <typename T>
struct image {
    constexpr image(); 
    constexpr image(int w, int h, const T& v =; 
    constexpr image(int w, int h, const T* v); 
    int width() const; 
    int height() const; 
    vec2i size() const; 
    void resize(int w, int h, const T& v =; 
    void assign(int w, int h, const T& v); 
    void set(const T& v); 
    T& operator[](const vec2i& ij); 
    const T& operator[](const vec2i& ij) const; 
    T& at(const vec2i& ij); 
    const T& at(const vec2i& ij) const; 
    T& at(int i, int j); 
    const T& at(int i, int j) const; 
    T* data(); 
    const T* data() const; 
}
~~~

Image of a specified type

- Members:
    - image():      empty image constructor
    - image():      image constructor
    - image():      image constructor
    - width():      width
    - height():      height
    - size():      size
    - resize():      reallocate memory
    - assign():      reallocate memory
    - set():      set values
    - operator[]():      element access
    - operator[]():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - data():      data access
    - data():      data access


### Typedef image1f

~~~ .cpp
using image1f = image<vec<float, 1>>;
~~~

1-dimensional float image

### Typedef image2f

~~~ .cpp
using image2f = image<vec<float, 2>>;
~~~

2-dimensional float image

### Typedef image3f

~~~ .cpp
using image3f = image<vec<float, 3>>;
~~~

3-dimensional float image

### Typedef image4f

~~~ .cpp
using image4f = image<vec<float, 4>>;
~~~

4-dimensional float image

### Typedef image4b

~~~ .cpp
using image4b = image<vec<byte, 4>>;
~~~

4-dimensional byte image

### Typedef imagef

~~~ .cpp
using imagef = image<float>;
~~~

float image

### Function image_lookup()

~~~ .cpp
template <typename T>
constexpr inline vec<T, 4> image_lookup(
    int width, int height, int ncomp, const T* img, int x, int y, T alpha = 0);
~~~

Lookup an image value from a generic image

### Function image_set()

~~~ .cpp
template <typename T>
constexpr inline void image_set(int width, int height, int ncomp, T* img, int x,
    int y, const vec<T, 4>& vv);
~~~

Set an image value for a generic image

### Function srgb_to_linear()

~~~ .cpp
inline vec3f srgb_to_linear(const vec3b& srgb);
~~~

Approximate conversion from srgb.

### Function srgb_to_linear()

~~~ .cpp
inline vec4f srgb_to_linear(const vec4b& srgb);
~~~

Conversion from srgb.

### Function tonemap_filmic()

~~~ .cpp
inline vec3f tonemap_filmic(const vec3f& hdr);
~~~

Tone map with a fitted filmic curve.

Implementation from
https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/

### Function tonemap_image()

~~~ .cpp
inline void tonemap_image(int width, int height, int ncomp, const float* hdr,
    byte* ldr, tonemap_type tm, float exposure, float gamma);
~~~

Tone mapping HDR to LDR images.

### Function image_over()

~~~ .cpp
inline void image_over(
    vec4f* img, int width, int height, int nlayers, vec4f** layers);
~~~

Image over operator

### Function image_over()

~~~ .cpp
inline void image_over(
    vec4b* img, int width, int height, int nlayers, vec4b** layers);
~~~

Image over operator

### Function hsv_to_rgb()

~~~ .cpp
inline vec4b hsv_to_rgb(const vec4b& hsv);
~~~

Convert HSV to RGB

Implementatkion from
http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both

### Function make_grid_image()

~~~ .cpp
inline image<vec4b> make_grid_image(int size, int tile = 64,
    const vec4b& c0 =;
~~~

Make a grid image

### Function make_checker_image()

~~~ .cpp
inline image<vec4b> make_checker_image(int size, int tile = 64,
    const vec4b& c0 =;
~~~

Make a checkerboard image

### Function make_gammaramp_image()

~~~ .cpp
inline image<vec4b> make_gammaramp_image(int size);
~~~

Make a gamma ramp image

### Function make_gammaramp_imagef()

~~~ .cpp
inline image<vec4f> make_gammaramp_imagef(int size);
~~~

Make a gamma ramp image

### Function make_uvgrid_image()

~~~ .cpp
inline image<vec4b> make_uvgrid_image(
    int size, int tile = 64, bool colored = true);
~~~

Make a uv colored grid

### Function make_recuvgrid_image()

~~~ .cpp
inline image<vec4b> make_recuvgrid_image(
    int size, int tile = 64, bool colored = true);
~~~

Make a uv recusive colored grid

### Function bump_to_normal_map()

~~~ .cpp
inline image<vec4b> bump_to_normal_map(
    const image<vec4b>& img, float scale = 1);
~~~

Comvert a bump map to a normal map.

