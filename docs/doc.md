## Basic math constants and functions

#### Typedef byte

    using ygl::byte = unsigned char;

Convenient typedef for bytes.

#### Typedef uint

    using ygl::uint = unsigned int;

Convenient typedef for unsigned ints.

#### Constant pif

    const auto ygl::pif = = 3.14159265f;

Pi (float).

#### Constant pi

    const auto ygl::pi = = 3.1415926535897932384626433832795;

Pi (double).

#### Constant flt_max

    const auto ygl::flt_max = = std::numeric_limits<float>::max();

Shortcat for float max value.

#### Constant flt_min

    const auto ygl::flt_min = = std::numeric_limits<float>::lowest();

Shortcat for float min value.

#### Constant flt_eps

    const auto ygl::flt_eps = = std::numeric_limits<float>::epsilon();

Shortcat for float epsilon.

#### Constant int_max

    const auto ygl::int_max = = std::numeric_limits<int>::max();

Shortcat for int max value.

#### Constant int_min

    const auto ygl::int_min = = std::numeric_limits<int>::min();

Shortcat for int min value.

#### Function sqrt()

    template<typename T>
    T ygl::sqrt(T a);

Square root.

#### Function pow()

    template<typename T, typename T1>
    auto ygl::pow(T a, T1 b);

Power.

#### Function exp()

    template<typename T>
    T ygl::exp(T a);

Exponential.

#### Function log()

    template<typename T>
    T ygl::log(T a);

Logarithm.

#### Function sin()

    template<typename T>
    T ygl::sin(T a);

Sine.

#### Function cos()

    template<typename T>
    T ygl::cos(T a);

Cosine.

#### Function tan()

    template<typename T>
    T ygl::tan(T a);

Tangent.

#### Function asin()

    template<typename T>
    T ygl::asin(T a);

Arc sine.

#### Function acos()

    template<typename T>
    T ygl::acos(T a);

Arc cosine.

#### Function atan()

    template<typename T>
    T ygl::atan(T a);

Arc tangent.

#### Function atan2()

    template<typename T, typename T1>
    auto ygl::atan2(T a, T1 b);

Arc tangent.

#### Function abs()

    template<typename T>
    T ygl::abs(T a);

Absolute value.

#### Function floor()

    template<typename T>
    T ygl::floor(T a);

Floor.

#### Function round()

    template<typename T>
    T ygl::round(T a);

Round.

#### Function min()

    template<typename T>
    T ygl::min(T x, T y);

Safe minimum value.

#### Function min()

    template<typename T>
    T ygl::min(std::initializer_list<T> vs);

Safe minimum value.

#### Function max()

    template<typename T>
    T ygl::max(T x, T y);

Safe maximum value.

#### Function max()

    template<typename T>
    T ygl::max(std::initializer_list<T> vs);

Safe maximum value.

#### Function clamp()

    template<typename T>
    T ygl::clamp(T x, T min_, T max_);

Clamp a value between a minimum and a maximum.

#### Function lerp()

    template<typename T, typename T1>
    T ygl::lerp(const T& a, const T& b, T1 u);

Linear interpolation.

#### Function bilerp()

    template<typename T, typename T1>
    float ygl::bilerp(const T& a, const T& b, const T& c, const T& d, T1 u, T1
    v);

Bilinear interpolation. Order is specified like quads counter-
clockwise, so a,b,c,d correspond to parameters (0,0), (0,1), (1,1),
(0,1).

#### Function pow2()

    int ygl::pow2(int x);

Integer power of two.

#### Function fastfloor()

    int ygl::fastfloor(float x);

Fast floor.

#### Function float_to_byte()

    byte ygl::float_to_byte(float x);

Safe float to byte conversion.

#### Function byte_to_float()

    float ygl::byte_to_float(byte x);

Safe byte to float conversion.

## Fixed-size vectors

#### Typedef vec1f

    using ygl::vec1f = vec<float, 1>;

1-dimensional float vector.

#### Typedef vec2f

    using ygl::vec2f = vec<float, 2>;

2-dimensional float vector.

#### Typedef vec3f

    using ygl::vec3f = vec<float, 3>;

3-dimensional float vector

#### Typedef vec4f

    using ygl::vec4f = vec<float, 4>;

4-dimensional float vector

#### Typedef vec1i

    using ygl::vec1i = vec<int, 1>;

1-dimensional int vector.

#### Typedef vec2i

    using ygl::vec2i = vec<int, 2>;

2-dimensional int vector.

#### Typedef vec3i

    using ygl::vec3i = vec<int, 3>;

3-dimensional int vector.

#### Typedef vec4i

    using ygl::vec4i = vec<int, 4>;

4-dimensional int vector.

#### Typedef vec4b

    using ygl::vec4b = vec<byte, 4>;

4-dimensional byte vector.

#### Constant zero1f

    const auto ygl::zero1f = = vec1f();

1-dimensional float zero vector.

#### Constant zero2f

    const auto ygl::zero2f = = vec2f();

2-dimensional float zero vector.

#### Constant zero3f

    const auto ygl::zero3f = = vec3f();

3-dimensional float zero vector.

#### Constant zero4f

    const auto ygl::zero4f = = vec4f();

4-dimensional float zero vector.

#### Constant zero1i

    const auto ygl::zero1i = = vec1i();

1-dimensional int zero vector.

#### Constant zero2i

    const auto ygl::zero2i = = vec2i();

2-dimensional int zero vector.

#### Constant zero3i

    const auto ygl::zero3i = = vec3i();

3-dimensional int zero vector.

#### Constant zero4i

    const auto ygl::zero4i = = vec4i();

4-dimensional int zero vector.

#### Constant zero4b

    const auto ygl::zero4b = = vec4b();

4-dimensional byte zero vector.

#### Function begin()

    template<typename T, int>
    T* ygl::begin(vec<T, N>& a);

Element iteration.

#### Function begin()

    template<typename T, int>
    const T* ygl::begin(const vec<T, N>& a);

Element iteration.

#### Function end()

    template<typename T, int>
    T* ygl::end(vec<T, N>& a);

Element iteration.

#### Function end()

    template<typename T, int>
    const T* ygl::end(const vec<T, N>& a);

Element iteration.

#### Function data()

    template<typename T, int>
    T* ygl::data(vec<T, N>& a);

Element access.

#### Function data()

    template<typename T, int>
    const T* ygl::data(const vec<T, N>& a);

Element access.

#### Function size()

    template<typename T, int>
    int ygl::size(vec<T, N>& a);

Number of elements.

#### Function empty()

    template<typename T, int>
    bool ygl::empty(vec<T, N>& a);

Empty check (always false for useful for templated code).

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const vec<T, 1>& a, const vec<T, 1>& b);

Vector equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const vec<T, 1>& a, const vec<T, 1>& b);

Vector inequality.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const vec<T, 2>& a, const vec<T, 2>& b);

Vector equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const vec<T, 2>& a, const vec<T, 2>& b);

Vector inequality.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const vec<T, 3>& a, const vec<T, 3>& b);

Vector equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const vec<T, 3>& a, const vec<T, 3>& b);

Vector inequality.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const vec<T, 4>& a, const vec<T, 4>& b);

Vector equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const vec<T, 4>& a, const vec<T, 4>& b);

Vector inequality.

#### Function operator<()

    template<typename T>
    bool ygl::operator<(const vec<T, 2>& a, const vec<T, 2>& b);

Vector comparison using lexicographic order, useful for map.

#### Function operator<()

    template<typename T>
    bool ygl::operator<(const vec<T, 3>& a, const vec<T, 3>& b);

Vector comparison using lexicographic order, useful for map.

#### Function operator<()

    template<typename T>
    bool ygl::operator<(const vec<T, 4>& a, const vec<T, 4>& b);

Vector comparison using lexicographic order, useful for map.

#### Function operator+()

    template<typename T>
    vec<T, 2> ygl::operator+(const vec<T, 2>& a);

Vector unary plus (for completeness).

#### Function operator-()

    template<typename T>
    vec<T, 2> ygl::operator-(const vec<T, 2>& a);

Vector negation.

#### Function operator+()

    template<typename T>
    vec<T, 2> ygl::operator+(const vec<T, 2>& a, const vec<T, 2>& b);

Vector sum.

#### Function operator-()

    template<typename T>
    vec<T, 2> ygl::operator-(const vec<T, 2>& a, const vec<T, 2>& b);

Vector difference.

#### Function operator*()

    template<typename T>
    vec<T, 2> ygl::operator*(const vec<T, 2>& a, const vec<T, 2>& b);

Vector scalar product.

#### Function operator*()

    template<typename T, typename T1>
    vec<T, 2> ygl::operator*(const vec<T, 2>& a, T1 b);

Vector scalar product.

#### Function operator*()

    template<typename T>
    vec<T, 2> ygl::operator*(float a, const vec<T, 2>& b);

Vector scalar product.

#### Function operator/()

    template<typename T>
    vec<T, 2> ygl::operator/(const vec<T, 2>& a, const vec<T, 2>& b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 2> ygl::operator/(const vec<T, 2>& a, T1 b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 2> ygl::operator/(T1 a, const vec<T, 2>& b);

Vector scalar division.

#### Function operator+()

    template<typename T>
    vec<T, 3> ygl::operator+(const vec<T, 3>& a);

Vector unary plus (for completeness).

#### Function operator-()

    template<typename T>
    vec<T, 3> ygl::operator-(const vec<T, 3>& a);

Vector negation.

#### Function operator+()

    template<typename T>
    vec<T, 3> ygl::operator+(const vec<T, 3>& a, const vec<T, 3>& b);

Vector sum.

#### Function operator-()

    template<typename T>
    vec<T, 3> ygl::operator-(const vec<T, 3>& a, const vec<T, 3>& b);

Vector operator -.

#### Function operator*()

    template<typename T>
    vec<T, 3> ygl::operator*(const vec<T, 3>& a, const vec<T, 3>& b);

Vector scalar product.

#### Function operator*()

    template<typename T, typename T1>
    vec<T, 3> ygl::operator*(const vec<T, 3>& a, T1 b);

Vector scalar product.

#### Function operator*()

    template<typename T, typename T1>
    vec<T, 3> ygl::operator*(T1 a, const vec<T, 3>& b);

Vector scalar product.

#### Function operator/()

    template<typename T>
    vec<T, 3> ygl::operator/(const vec<T, 3>& a, const vec<T, 3>& b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 3> ygl::operator/(const vec<T, 3>& a, T1 b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 3> ygl::operator/(T1 a, const vec<T, 3>& b);

Vector scalar division.

#### Function operator+()

    template<typename T>
    vec<T, 4> ygl::operator+(const vec<T, 4>& a);

Vector unary plus (for completeness).

#### Function operator-()

    template<typename T>
    vec<T, 4> ygl::operator-(const vec<T, 4>& a);

Vector negation.

#### Function operator+()

    template<typename T>
    vec<T, 4> ygl::operator+(const vec<T, 4>& a, const vec<T, 4>& b);

Vector sum.

#### Function operator-()

    template<typename T>
    vec<T, 4> ygl::operator-(const vec<T, 4>& a, const vec<T, 4>& b);

Vector difference.

#### Function operator*()

    template<typename T>
    vec<T, 4> ygl::operator*(const vec<T, 4>& a, const vec<T, 4>& b);

Vector scalar product.

#### Function operator*()

    template<typename T>
    vec<T, 4> ygl::operator*(const vec<T, 4>& a, float b);

Vector scalar product.

#### Function operator*()

    template<typename T>
    vec<T, 4> ygl::operator*(float a, const vec<T, 4>& b);

Vector scalar product.

#### Function operator/()

    template<typename T>
    vec<T, 4> ygl::operator/(const vec<T, 4>& a, const vec<T, 4>& b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 4> ygl::operator/(const vec<T, 4>& a, T1 b);

Vector scalar division.

#### Function operator/()

    template<typename T, typename T1>
    vec<T, 4> ygl::operator/(T1 a, const vec<T, 4>& b);

Vector scalar division.

#### Function operator+=()

    template<typename T, int>
    vec<T, N>& ygl::operator+=(vec<T, N>& a, const vec<T, N>& b);

Vector assignment.

#### Function operator-=()

    template<typename T, int>
    vec<T, N>& ygl::operator-=(vec<T, 2>& a, const vec<T, N>& b);

Vector assignment.

#### Function operator*=()

    template<typename T, int>
    vec<T, N>& ygl::operator*=(vec<T, N>& a, const vec<T, N>& b);

Vector assignment.

#### Function operator*=()

    template<typename T, int, typename T1>
    vec<T, N>& ygl::operator*=(vec<T, N>& a, T1 b);

Vector assignment.

#### Function operator/=()

    template<typename T, int>
    vec<T, N>& ygl::operator/=(vec<T, N>& a, const vec<T, N>& b);

Vector assignment.

#### Function operator/=()

    template<typename T, int, typename T1>
    vec<T, N>& ygl::operator/=(vec<T, N>& a, T1 b);

Vector assignment.

#### Function dot()

    template<typename T>
    T ygl::dot(const vec<T, 2>& a, const vec<T, 2>& b);

Vector dot product.

#### Function dot()

    template<typename T>
    T ygl::dot(const vec<T, 3>& a, const vec<T, 3>& b);

Vector dot product.

#### Function dot()

    template<typename T>
    T ygl::dot(const vec<T, 4>& a, const vec<T, 4>& b);

Vector dot product.

#### Function cross()

    template<typename T>
    T ygl::cross(const vec<T, 2>& a, const vec<T, 2>& b);

Vector cross product.

#### Function cross()

    template<typename T>
    vec<T, 3> ygl::cross(const vec<T, 3>& a, const vec<T, 3>& b);

Vector cross product.

#### Function length()

    template<typename T, int>
    T ygl::length(const vec<T, N>& a);

Vector length.

#### Function normalize()

    template<typename T, int>
    vec<T, N> ygl::normalize(const vec<T, N>& a);

Vector normalization.

#### Function angle()

    template<typename T, int>
    T ygl::angle(const vec<T, N>& a, const vec<T, N>& b);

Angle between vectors.

#### Function slerp()

    template<typename T, int, typename T1>
    vec<T, N> ygl::slerp(const vec<T, N>& a, const vec<T, N>& b, T1 u);

Vector spherical linear interpolation (vectors have to be normalized).

#### Function orthogonal()

    template<typename T>
    vec<T, 3> ygl::orthogonal(const vec<T, 3>& v);

Orthogonal vector.

#### Function orthonormalize()

    template<typename T>
    vec<T, 3> ygl::orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b);

Orthonormalize two vectors.

#### Function reflect()

    template<typename T>
    vec<T, 3> ygl::reflect(const vec<T, 3>& w, const vec<T, 3>& n);

Reflected vector.

#### Function refract()

    template<typename T>
    vec<T, 3> ygl::refract(const vec<T, 3>& w, const vec<T, 3>& n, T eta);

Refracted vector.

#### Function clamp()

    template<typename T, typename T1>
    vec<T, 2> ygl::clamp(const vec<T, 2>& x, T1 min, T1 max);

Component-wise clamp.

#### Function clamp()

    template<typename T, typename T1>
    vec<T, 3> ygl::clamp(const vec<T, 3>& x, T1 min, T1 max);

Component-wise clamp.

#### Function clamp()

    template<typename T, typename T1>
    vec<T, 4> ygl::clamp(const vec<T, 4>& x, T1 min, T1 max);

Component-wise clamp.

#### Function clamplen()

    template<typename T, int, typename T1>
    vec<T, N> ygl::clamplen(const vec<T, N>& x, T1 max);

Clamp a vector to a maximum length.

#### Function min_element()

    template<typename T, int>
    int ygl::min_element(const vec<T, N>& a);

Index of minimum element.

#### Function min_element_value()

    template<typename T, int>
    T ygl::min_element_value(const vec<T, N>& a);

Value of minimum element.

#### Function max_element()

    template<typename T, int>
    int ygl::max_element(const vec<T, N>& a);

Index of maximum element.

#### Function max_element_value()

    template<typename T, int>
    T ygl::max_element_value(const vec<T, N>& a);

Value of maximum element.

#### Function float_to_byte()

    vec4b ygl::float_to_byte(const vec4f& a);

Element-wise float to byte conversion.

#### Function byte_to_float()

    vec4f ygl::byte_to_float(const vec4b& a);

Element-wise byte to float conversion.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const vec<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, vec<T, N>& a);

Stream read.

## Fixed-size matrices

#### Typedef mat2f

    using ygl::mat2f = mat<float, 2>;

2-dimensional float matrix.

#### Typedef mat3f

    using ygl::mat3f = mat<float, 3>;

3-dimensional float matrix.

#### Typedef mat4f

    using ygl::mat4f = mat<float, 4>;

4-dimensional float matrix.

#### Constant identity_mat2f

    const auto ygl::identity_mat2f = = mat2f();

2-dimensional float identity matrix.

#### Constant identity_mat3f

    const auto ygl::identity_mat3f = = mat3f();

3-dimensional float identity matrix.

#### Constant identity_mat4f

    const auto ygl::identity_mat4f = = mat4f();

4-dimensional float identity matrix.

#### Function begin()

    template<typename T, int>
    vec<T, N>* ygl::begin(mat<T, N>& m);

Column iteration.

#### Function end()

    template<typename T, int>
    vec<T, N>* ygl::end(mat<T, N>& m);

Column iteration.

#### Function begin()

    template<typename T, int>
    const vec<T, N>* ygl::begin(const mat<T, N>& m);

Column iteration.

#### Function end()

    template<typename T, int>
    const vec<T, N>* ygl::end(const mat<T, N>& m);

Column iteration.

#### Function data()

    template<typename T, int>
    vec<T, N>* ygl::data(mat<T, N>& m);

Column access.

#### Function data()

    template<typename T, int>
    const vec<T, N>* ygl::data(const mat<T, N>& m);

Column access.

#### Function size()

    template<typename T, int>
    int ygl::size(mat<T, N>& a);

Number of columns.

#### Function empty()

    template<typename T, int>
    bool ygl::empty(mat<T, N>& a);

Empty check (always false for useful for templated code).

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const mat<T, 2>& a, const mat<T, 2>& b);

Matrix equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const mat<T, 2>& a, const mat<T, 2>& b);

Matrix inequality.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const mat<T, 3>& a, const mat<T, 3>& b);

Matrix equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const mat<T, 3>& a, const mat<T, 3>& b);

Matrix inequality.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const mat<T, 4>& a, const mat<T, 4>& b);

Matrix equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const mat<T, 4>& a, const mat<T, 4>& b);

Matrix inequality.

#### Function operator+()

    template<typename T>
    mat<T, 2> ygl::operator+(const mat<T, 2>& a, const mat<T, 2>& b);

Matrix sum.

#### Function operator*()

    template<typename T, typename T1>
    mat<T, 2> ygl::operator*(const mat<T, 2>& a, T1 b);

Matrix scalar product.

#### Function operator/()

    template<typename T, typename T1>
    mat<T, 2> ygl::operator/(const mat<T, 2>& a, T1 b);

matrix scalar division.

#### Function operator*()

    template<typename T>
    vec<T, 2> ygl::operator*(const mat<T, 2>& a, const vec<T, 2>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    vec<T, 2> ygl::operator*(const vec<T, 2>& a, const mat<T, 2>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    mat<T, 2> ygl::operator*(const mat<T, 2>& a, const mat<T, 2>& b);

Matrix-matrix product.

#### Function operator+()

    template<typename T>
    mat<T, 3> ygl::operator+(const mat<T, 3>& a, const mat<T, 3>& b);

Matrix sum.

#### Function operator*()

    template<typename T, typename T1>
    mat<T, 3> ygl::operator*(const mat<T, 3>& a, T1 b);

Matrix scalar product.

#### Function operator/()

    template<typename T, typename T1>
    mat<T, 3> ygl::operator/(const mat<T, 3>& a, T1 b);

Matrix scalar division.

#### Function operator*()

    template<typename T>
    vec<T, 3> ygl::operator*(const mat<T, 3>& a, const vec<T, 3>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    vec<T, 3> ygl::operator*(const vec<T, 3>& a, const mat<T, 3>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    mat<T, 3> ygl::operator*(const mat<T, 3>& a, const mat<T, 3>& b);

Matrix-matrix product.

#### Function operator+()

    template<typename T>
    mat<T, 4> ygl::operator+(const mat<T, 4>& a, const mat<T, 4>& b);

Matrix sum.

#### Function operator*()

    template<typename T, typename T1>
    mat<T, 4> ygl::operator*(const mat<T, 4>& a, T1 b);

Matrix scalar product.

#### Function operator/()

    template<typename T, typename T1>
    mat<T, 4> ygl::operator/(const mat<T, 4>& a, T1 b);

Matrix scalar division.

#### Function operator*()

    template<typename T>
    vec<T, 4> ygl::operator*(const mat<T, 4>& a, const vec<T, 4>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    vec<T, 4> ygl::operator*(const vec<T, 4>& a, const mat<T, 4>& b);

Matrix-vector product.

#### Function operator*()

    template<typename T>
    mat<T, 4> ygl::operator*(const mat<T, 4>& a, const mat<T, 4>& b);

Matrix-matrix product.

#### Function operator+=()

    template<typename T, int>
    mat<T, N>& ygl::operator+=(mat<T, N>& a, const mat<T, N>& b);

Matrix assignment.

#### Function operator*=()

    template<typename T, int>
    mat<T, N>& ygl::operator*=(mat<T, N>& a, const mat<T, N>& b);

Matrix assignment.

#### Function operator*=()

    template<typename T, int, typename T1>
    mat<T, N>& ygl::operator*=(mat<T, N>& a, T1 b);

Matrix assignment.

#### Function operator/=()

    template<typename T, int, typename T1>
    mat<T, N>& ygl::operator/=(mat<T, N>& a, T1 b);

Matrix assignment.

#### Function mat_diagonal()

    template<typename T>
    vec<T, 2> ygl::mat_diagonal(const mat<T, 2>& a);

Matrix diagonal.

#### Function mat_diagonal()

    template<typename T>
    vec<T, 3> ygl::mat_diagonal(const mat<T, 3>& a);

Matrix diagonal.

#### Function mat_diagonal()

    template<typename T>
    vec<T, 4> ygl::mat_diagonal(const mat<T, 4>& a);

Matrix diagonal.

#### Function transpose()

    template<typename T>
    mat<T, 2> ygl::transpose(const mat<T, 2>& a);

Matrix transpose.

#### Function transpose()

    template<typename T>
    mat<T, 3> ygl::transpose(const mat<T, 3>& a);

Matrix transpose.

#### Function transpose()

    template<typename T>
    mat<T, 4> ygl::transpose(const mat<T, 4>& a);

Matrix transpose.

#### Function adjugate()

    template<typename T>
    mat<T, 2> ygl::adjugate(const mat<T, 2>& a);

Matrix adjugate.

#### Function adjugate()

    template<typename T>
    mat<T, 3> ygl::adjugate(const mat<T, 3>& a);

Matrix adjugate.

#### Function adjugate()

    template<typename T>
    mat<T, 4> ygl::adjugate(const mat<T, 4>& a);

Matrix adjugate.

#### Function determinant()

    template<typename T>
    T ygl::determinant(const mat<T, 2>& a);

Matrix determinant.

#### Function determinant()

    template<typename T>
    T ygl::determinant(const mat<T, 3>& a);

Matrix determinant.

#### Function determinant()

    template<typename T>
    T ygl::determinant(const mat<T, 4>& a);

Matrix determinant.

#### Function inverse()

    template<typename T, int>
    mat<T, N> ygl::inverse(const mat<T, N>& a);

Matrix inverse.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const mat<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, mat<T, N>& a);

Stream read.

## Rigid-body frames

#### Typedef frame3f

    using ygl::frame3f = frame<float, 3>;

3-dimensional float frame.

#### Constant identity_frame3f

    const auto ygl::identity_frame3f = =     frame3f{{1, 0, 0}, {0, 1, 0}, {0,
    0, 1}, {0, 0, 0}};

Indentity frame.

#### Function begin()

    template<typename T, int>
    vec<T, N>* ygl::begin(frame<T, N>& a);

Element/column iteration.

#### Function begin()

    template<typename T, int>
    const vec<T, N>* ygl::begin(const frame<T, N>& a);

Element/column iteration.

#### Function end()

    template<typename T, int>
    vec<T, N>* ygl::end(frame<T, N>& a);

Element/column iteration.

#### Function end()

    template<typename T, int>
    const vec<T, N>* ygl::end(const frame<T, N>& a);

Element/column iteration.

#### Function data()

    template<typename T, int>
    vec<T, N>* ygl::data(frame<T, N>& a);

Element/column access.

#### Function data()

    template<typename T, int>
    const vec<T, N>* ygl::data(const frame<T, N>& a);

Element/column access.

#### Function size()

    template<typename T, int>
    int ygl::size(frame<T, N>& a);

Number of columns in the underlying affine matrix.

#### Function empty()

    template<typename T, int>
    bool ygl::empty(frame<T, N>& a);

Empty check (always false for useful for templated code).

#### Function make_frame_fromz()

    template<typename T>
    frame<T, 3> ygl::make_frame_fromz(const vec<T, 3>& o, const vec<T, 3>&
    z_);



#### Function make_frame_fromzx()

    template<typename T>
    frame<T, 3> ygl::make_frame_fromzx(const vec<T, 3>& o, const vec<T, 3>&
    z_, const vec<T, 3>& x_);



#### Function frame_to_mat()

    template<typename T>
    mat<T, 4> ygl::frame_to_mat(const frame<T, 3>& a);

Frame to matrix conversion.

#### Function mat_to_frame()

    template<typename T>
    frame<T, 3> ygl::mat_to_frame(const mat<T, 4>& a);

Matrix to frame conversion.

#### Function frame_pos()

    template<typename T, int>
    vec<T, N>& ygl::frame_pos(frame<T, N>& a);

Frame origin.

#### Function frame_pos()

    template<typename T, int>
    const vec<T, N>& ygl::frame_pos(const frame<T, N>& a);

Frame origin.

#### Function frame_rot()

    template<typename T, int>
    mat<T, 3>& ygl::frame_rot(frame<T, N>& a);

Frame rotation.

#### Function frame_rot()

    template<typename T, int>
    const mat<T, 3>& ygl::frame_rot(const frame<T, N>& a);

Frame rotation.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const frame<T, 3>& a, const frame<T, 3>& b);

Frame equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const frame<T, 3>& a, const frame<T, 3>& b);

Frame inequality.

#### Function operator*()

    template<typename T>
    frame<T, 3> ygl::operator*(const frame<T, 3>& a, const frame<T, 3>& b);

Frame composition, equivalent to affine matrix product.

#### Function inverse()

    template<typename T>
    frame<T, 3> ygl::inverse(const frame<T, 3>& a);

Frame inverse, equivalent to rigid affine inverse.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const frame<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, frame<T, N>& a);

Stream read.

## Quaternions

#### Typedef quat4f

    using ygl::quat4f = quat<float, 4>;

4-dimensional float quaternion.

#### Constant identity_quat4f

    const auto ygl::identity_quat4f = = quat4f{0, 0, 0, 1};

Float identity quaternion.

#### Function begin()

    template<typename T, int>
    T* ygl::begin(quat<T, N>& a);

Element iteration.

#### Function begin()

    template<typename T, int>
    const T* ygl::begin(const quat<T, N>& a);

Element iteration.

#### Function end()

    template<typename T, int>
    T* ygl::end(quat<T, N>& a);

Element iteration.

#### Function end()

    template<typename T, int>
    const T* ygl::end(const quat<T, N>& a);

Element iteration.

#### Function data()

    template<typename T, int>
    T* ygl::data(quat<T, N>& a);

Element access.

#### Function data()

    template<typename T, int>
    const T* ygl::data(const quat<T, N>& a);

Element access.

#### Function size()

    template<typename T, int>
    int ygl::size(quat<T, N>& a);

Number of elements.

#### Function empty()

    template<typename T, int>
    bool ygl::empty(quat<T, N>& a);

Empty check (always false for useful for templated code).

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const quat<T, 4>& a, const quat<T, 4>& b);

Quaternion equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const quat<T, 4>& a, const quat<T, 4>& b);

Quaternion inequality.

#### Function operator+()

    template<typename T>
    quat<T, 4> ygl::operator+(const quat<T, 4>& a, const quat<T, 4>& b);

Quaternion sum.

#### Function operator*()

    template<typename T, typename T1>
    quat<T, 4> ygl::operator*(const quat<T, 4>& a, T1 b);

Quaternion scalar product.

#### Function operator/()

    template<typename T, typename T1>
    quat<T, 4> ygl::operator/(const quat<T, 4>& a, T1 b);

Quaternion scalar division.

#### Function operator*()

    template<typename T>
    quat<T, 4> ygl::operator*(const quat<T, 4>& a, const quat<T, 4>& b);

Quaternion product.

#### Function conjugate()

    template<typename T>
    quat<T, 4> ygl::conjugate(const quat<T, 4>& v);

Quaternion conjugate.

#### Function inverse()

    template<typename T>
    quat<T, 4> ygl::inverse(const quat<T, 4>& v);

Quaternion inverse.

#### Function normalize()

    template<typename T>
    quat<T, 4> ygl::normalize(const quat<T, 4>& v);

Quaternion normalization.

#### Function slerp()

    template<typename T, typename T1>
    quat<T, 4> ygl::slerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t);

Quaternion spherical linear interpolation.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const quat<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, quat<T, N>& a);

Stream read.

## Axis-aligned bounding boxes

#### Typedef bbox1f

    using ygl::bbox1f = bbox<float, 1>;

1-dimensional float bounding box.

#### Typedef bbox2f

    using ygl::bbox2f = bbox<float, 2>;

2-dimensional float bounding box.

#### Typedef bbox3f

    using ygl::bbox3f = bbox<float, 3>;

3-dimensional float bounding box.

#### Typedef bbox4f

    using ygl::bbox4f = bbox<float, 4>;

4-dimensional float bounding box.

#### Constant invalid_bbox1f

    const auto ygl::invalid_bbox1f = = bbox1f();

1-dimensional float empty bbox.

#### Constant invalid_bbox2f

    const auto ygl::invalid_bbox2f = = bbox2f();

2-dimensional float empty bbox.

#### Constant invalid_bbox3f

    const auto ygl::invalid_bbox3f = = bbox3f();

3-dimensional float empty bbox.

#### Constant invalid_bbox4f

    const auto ygl::invalid_bbox4f = = bbox4f();

4-dimensional float empty bbox.

#### Function operator==()

    template<typename T, int>
    bool ygl::operator==(const bbox<T, N>& a, const bbox<T, N>& b);

Bounding box equality.

#### Function operator!=()

    template<typename T, int>
    bool ygl::operator!=(const bbox<T, N>& a, const bbox<T, N>& b);

Bounding box inequality.

#### Function bbox_center()

    template<typename T, int>
    vec<T, N> ygl::bbox_center(const bbox<T, N>& a);

Bounding box center.

#### Function bbox_diagonal()

    template<typename T, int>
    vec<T, N> ygl::bbox_diagonal(const bbox<T, N>& a);

Bounding box diagonal.

#### Function expand()

    template<typename T>
    bbox<T, 1> ygl::expand(const bbox<T, 1>& a, T b);

Expands a bounding box with a point.

#### Function expand()

    template<typename T>
    bbox<T, 1> ygl::expand(const bbox<T, 1>& a, const vec<T, 1>& b);

Expands a bounding box with a point.

#### Function expand()

    template<typename T>
    bbox<T, 2> ygl::expand(const bbox<T, 2>& a, const vec<T, 2>& b);

Expands a bounding box with a point.

#### Function expand()

    template<typename T>
    bbox<T, 3> ygl::expand(const bbox<T, 3>& a, const vec<T, 3>& b);

Expands a bounding box with a point.

#### Function expand()

    template<typename T>
    bbox<T, 4> ygl::expand(const bbox<T, 4>& a, const vec<T, 4>& b);

Expands a bounding box with a point.

#### Function expand()

    template<typename T>
    bbox<T, 1> ygl::expand(const bbox<T, 1>& a, const bbox<T, 1>& b);

Expands a bounding box with a bounding box.

#### Function expand()

    template<typename T>
    bbox<T, 2> ygl::expand(const bbox<T, 2>& a, const bbox<T, 2>& b);

Expands a bounding box with a bounding box.

#### Function expand()

    template<typename T>
    bbox<T, 3> ygl::expand(const bbox<T, 3>& a, const bbox<T, 3>& b);

Expands a bounding box with a bounding box.

#### Function expand()

    template<typename T>
    bbox<T, 4> ygl::expand(const bbox<T, 4>& a, const bbox<T, 4>& b);

Expands a bounding box with a bounding box.

#### Function contains()

    template<typename T, int>
    bool ygl::contains(const bbox<T, N>& a, const vec<T, N>& b);

Check if a bounding box contains a point.

#### Function contains()

    template<typename T, int>
    bool ygl::contains(const bbox<T, 3>& a, const bbox<T, 3>& b);

Check if a bounding box contains a bounding box.

#### Function operator+=()

    template<typename T, int>
    bbox<T, N>& ygl::operator+=(bbox<T, N>& a, const vec<T, N>& b);

Expands a bounding box with a point.

#### Function operator+=()

    template<typename T, int>
    bbox<T, N>& ygl::operator+=(bbox<T, N>& a, const bbox<T, N>& b);

Expands a bounding box with a bounding box.

#### Function make_bbox()

    template<typename T, int>
    bbox<T, N> ygl::make_bbox(int count, const vec<T, N>* v);

Initialize a bonding box from a list of points.

#### Function make_bbox()

    template<typename T, int>
    bbox<T, N> ygl::make_bbox(const std::initializer_list<vec<T, N>>& v);

Initialize a bonding box from a list of points.

#### Function bbox_corners()

    template<typename T>
    std::array<vec<T, 2>, 4> ygl::bbox_corners(const bbox<T, 2>& a);

Computes the corners of a bounding boxes.

#### Function bbox_corners()

    template<typename T>
    std::array<vec<T, 3>, 8> ygl::bbox_corners(const bbox<T, 3>& a);

Computes the corners of a bounding boxes.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const bbox<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, bbox<T, N>& a);

Stream read.

## Primitive bounding boxes

#### Function point_bbox()

    template<typename T, typename T1>
    bbox<T, 3> ygl::point_bbox(const vec<T, 3>& p, T1 r=0);

Point bounds.

#### Function line_bbox()

    template<typename T, typename T1>
    bbox<T, 3> ygl::line_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1, T1
    r0=0, T1 r1=0);

Line bounds.

#### Function triangle_bbox()

    template<typename T>
    bbox<T, 3> ygl::triangle_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2);

Triangle bounds.

#### Function quad_bbox()

    template<typename T>
    bbox<T, 3> ygl::quad_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1, const
    vec<T, 3>& v2, const vec<T, 3>& v3);

Quad bounds.

#### Function tetrahedron_bbox()

    template<typename T, typename T1>
    bbox<T, 3> ygl::tetrahedron_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3);

Tetrahedron bounds.

## Rays

#### Typedef ray2f

    using ygl::ray2f = ray<float, 2>;

2-dimensional float ray.

#### Typedef ray3f

    using ygl::ray3f = ray<float, 3>;

3-dimensional float ray.

#### Typedef ray4f

    using ygl::ray4f = ray<float, 4>;

4-dimensional float ray.

#### Function make_ray()

    template<typename T, int>
    ray<T, N> ygl::make_ray(const vec<T, N>& o, const vec<T, N>& d, T
    eps=1e-4f);

Construct a ray using a default epsilon.

#### Function make_segment()

    template<typename T, int>
    ray<T, N> ygl::make_segment(const vec<T, N>& p1, const vec<T, N>& p2, T
    eps=1e-4f);

Construct a ray segment using a default epsilon.

#### Function operator<<()

    template<typename T, int>
    std::ostream& ygl::operator<<(std::ostream& os, const ray<T, N>& a);

Stream write.

#### Function operator>>()

    template<typename T, int>
    std::istream& ygl::operator>>(std::istream& is, ray<T, N>& a);

Stream read.

## Transforms

#### Function transform_point()

    template<typename T>
    vec<T, 2> ygl::transform_point(const mat<T, 3>& a, const vec<T, 2>& b);

Transforms a point by a matrix.

#### Function transform_point()

    template<typename T>
    vec<T, 3> ygl::transform_point(const mat<T, 4>& a, const vec<T, 3>& b);

Transforms a point by a matrix.

#### Function transform_vector()

    template<typename T>
    vec<T, 2> ygl::transform_vector(const mat<T, 3>& a, const vec<T, 2>& b);

Transforms a vector by a matrix.

#### Function transform_vector()

    template<typename T>
    vec<T, 3> ygl::transform_vector(const mat<T, 4>& a, const vec<T, 3>& b);

Transforms a vector by a matrix.

#### Function transform_direction()

    template<typename T, int>
    vec<T, N> ygl::transform_direction(const mat<T, N+1>& a, const vec<T, N>&
    b);

Transforms a direction by a matrix.

#### Function transform_ray()

    template<typename T, int>
    ray<T, N> ygl::transform_ray(const mat<T, N+1>& a, const ray<T, N>& b);

Transforms a ray by a matrix, leaving the direction not normalized.

#### Function transform_bbox()

    template<typename T, int>
    bbox<T, N> ygl::transform_bbox(const mat<T, N+1>& a, const bbox<T, N>& b);

transforms a bbox by a matrix

#### Function transform_point()

    template<typename T>
    vec<T, 2> ygl::transform_point(const frame<T, 2>& a, const vec<T, 2>& b);

Transforms a point by a frame, i.e. an affine transform.

#### Function transform_point()

    template<typename T>
    vec<T, 3> ygl::transform_point(const frame<T, 3>& a, const vec<T, 3>& b);

Transforms a point by a frame, i.e. an affine transform.

#### Function transform_vector()

    template<typename T>
    vec<T, 2> ygl::transform_vector(const frame<T, 2>& a, const vec<T, 2>& b);

Transforms a vector by a frame, i.e. an affine transform.

#### Function transform_vector()

    template<typename T>
    vec<T, 3> ygl::transform_vector(const frame<T, 3>& a, const vec<T, 3>& b);

Transforms a vector by a frame, i.e. an affine transform.

#### Function transform_direction()

    template<typename T, int>
    vec<T, N> ygl::transform_direction(const frame<T, N>& a, const vec<T, N>&
    b);

Transforms a direction by a frame, i.e. an affine transform.

#### Function transform_frame()

    template<typename T, int>
    frame<T, N> ygl::transform_frame(const frame<T, N>& a, const frame<T, N>&
    b);

Transforms a frame by a frame, i.e. an affine transform.

#### Function transform_ray()

    template<typename T, int>
    ray<T, N> ygl::transform_ray(const frame<T, 3>& a, const ray<T, N>& b);

Transforms a ray by a frame, i.e. an affine transform.

#### Function transform_bbox()

    template<typename T, int>
    bbox<T, N> ygl::transform_bbox(const frame<T, N>& a, const bbox<T, N>& b);

Transforms a bbox by a frame, i.e. an affine transform.

#### Function transform_bbox()

    template<typename T>
    bbox<T, 3> ygl::transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b);

Transforms a bbox by a frame, i.e. an affine transform.

#### Function transform_point_inverse()

    template<typename T>
    vec<T, 2> ygl::transform_point_inverse(const frame<T, 2>& a, const vec<T,
    2>& b);

Inverse transforms a point by a frame, assuming a rigid transform.

#### Function transform_point_inverse()

    template<typename T>
    vec<T, 3> ygl::transform_point_inverse(const frame<T, 3>& a, const vec<T,
    3>& b);

Inverse transforms a point by a frame, assuming a rigid transform.

#### Function transform_vector_inverse()

    template<typename T>
    vec<T, 2> ygl::transform_vector_inverse(const frame<T, 2>& a, const vec<T,
    2>& b);

Inverse transforms a vector by a frame, assuming a rigid transform.

#### Function transform_vector_inverse()

    template<typename T>
    vec<T, 3> ygl::transform_vector_inverse(const frame<T, 3>& a, const vec<T,
    3>& b);

Inverse transforms a vector by a frame, assuming a rigid transform.

#### Function transform_direction_inverse()

    template<typename T, int>
    vec<T, N> ygl::transform_direction_inverse(const frame<T, N>& a, const
    vec<T, N>& b);

Inverse transforms a direction by a frame, assuming a rigid transform.

#### Function transform_ray_inverse()

    template<typename T, int>
    ray<T, N> ygl::transform_ray_inverse(const frame<T, N>& a, const ray<T,
    N>& b);

Inverse transforms a direction by a frame, assuming a rigid transform.

#### Function transform_bbox_inverse()

    template<typename T, int>
    bbox<T, N> ygl::transform_bbox_inverse(const frame<T, N>& a, const bbox<T,
    N>& b);

Inverse transforms a bbox by a frame, assuming a rigid transform.

#### Function translation_frame()

    template<typename T>
    frame<T, 3> ygl::translation_frame(const vec<T, 3>& a);

Translation affine transform.

#### Function scaling_frame()

    template<typename T>
    frame<T, 3> ygl::scaling_frame(const vec<T, 3>& a);

Scaling affine transform; this is not rigid and here for symmatry of
API.

#### Function rotation_frame()

    template<typename T, typename T1>
    frame<T, 3> ygl::rotation_frame(const vec<T, 3>& axis, T1 angle);

Rotation affine transform.

#### Function rotation_frame()

    template<typename T>
    frame<T, 3> ygl::rotation_frame(const quat<T, 4>& v);

Rotation affine transform.

#### Function lookat_frame()

    template<typename T>
    frame<T, 3> ygl::lookat_frame(const vec<T, 3>& eye, const vec<T, 3>&
    center, const vec<T, 3>& up, bool inv_xz=false);

OpenGL lookat frame. Z-axis can be inverted with inv_xz.

#### Function frustum_mat()

    template<typename T>
    mat<T, 4> ygl::frustum_mat(T l, T r, T b, T t, T n, T f);

OpenGL frustum matrix.

#### Function ortho_mat()

    template<typename T>
    mat<T, 4> ygl::ortho_mat(T l, T r, T b, T t, T n, T f);

OpenGL orthographic matrix.

#### Function ortho2d_mat()

    template<typename T>
    mat<T, 4> ygl::ortho2d_mat(T left, T right, T bottom, T top);

OpenGL orthographic 2D matrix.

#### Function ortho_mat()

    template<typename T>
    mat<T, 4> ygl::ortho_mat(T xmag, T ymag, T near, T far);

OpenGL orthographic matrix.

#### Function perspective_mat()

    template<typename T>
    mat<T, 4> ygl::perspective_mat(T fovy, T aspect, T near, T far);

OpenGL perspective matrix.

#### Function perspective_mat()

    template<typename T>
    mat<T, 4> ygl::perspective_mat(T fovy, T aspect, T near);

OpenGL infinite perspective matrix.

#### Function rotation_axisangle()

    template<typename T>
    std::pair<vec<T, 4>, T> ygl::rotation_axisangle(const quat<T, 4>& a);

Rotation affine transform.

#### Function rotation_quat()

    template<typename T>
    quat<T, 4> ygl::rotation_quat(const vec<T, 3>& axis, T angle);

Axis-angle to quaternion conversion.

#### Function rotation_quat()

    template<typename T>
    quat<T, 4> ygl::rotation_quat(const mat<T, 3>& m_);

Rotation matrix to quaternion conversion.

#### Function decompose_frame()

    template<typename T>
    std::tuple<vec<T, 3>, mat<T, 3>, vec<T, 3>> ygl::decompose_frame(const
    frame<T, 3>& m);

Decompose an affine matrix into translation, rotation, scale. Assumes
there is no shear.

#### Function compose_frame()

    template<typename T>
    frame<T, 4> ygl::compose_frame(const vec<T, 3>& translation, const mat<T,
    3>& rotation, const vec<T, 3>& scale);

Decompose an affine matrix into translation, rotation, scale. Assumes
there is no shear and the matrix is affine.

#### Function compose_frame()

    template<typename T>
    frame<T, 4> ygl::compose_frame(const vec<T, 3>& translation, const quat<T,
    4>& rotation, const vec<T, 3>& scale);

Decompose an affine matrix into translation, rotation, scale. Assumes
there is no shear and the matrix is affine.

## User interface utilities

#### Function camera_turntable()

    void ygl::camera_turntable(vec3f& from, vec3f& to, vec3f& up, const vec2f&
    rotate, float dolly, const vec2f& pan);

Turntable for UI navigation.

#### Function camera_turntable()

    void ygl::camera_turntable(frame3f& frame, float& focus, const vec2f&
    rotate, float dolly, const vec2f& pan);

Turntable for UI navigation.

#### Function camera_fps()

    void ygl::camera_fps(frame3f& frame, const vec3f& transl, const vec2f&
    rotate);

FPS camera for UI navigation.

## Random number generation

#### Function advance_rng()

    uint32_t ygl::advance_rng(rng_pcg32& rng);

Next random number.

#### Function advance_rng()

    void ygl::advance_rng(rng_pcg32& rng, uint64_t delta);

Multi-step advance function (jump-ahead, jump-back).

#### Function advance_rng()

    void ygl::advance_rng(rng_pcg32& rng, int64_t delta);

Multi-step advance function (jump-ahead, jump-back).

#### Function seed_rng()

    void ygl::seed_rng(rng_pcg32& rng, uint64_t state, uint64_t seq=1);

Seeds a random number generator with a state state from the sequence
seq.

#### Function init_rng()

    rng_pcg32 ygl::init_rng(uint64_t state, uint64_t seq=1);

Init a random number generator with a state state from the sequence
seq.

#### Function next_rand1i()

    uint32_t ygl::next_rand1i(rng_pcg32& rng, uint32_t n);

Next random uint in [0,n) range with proper weighting.

#### Function next_rand1f()

    float ygl::next_rand1f(rng_pcg32& rng);

Next random float in [0,1).

#### Function next_rand1f()

    float ygl::next_rand1f(rng_pcg32& rng, float a, float b);

Next random float in [a,b).

#### Function next_rand2f()

    vec2f ygl::next_rand2f(rng_pcg32& rng);

Next random float2 in [0,1)x[0,1).

#### Function next_rand2f()

    vec2f ygl::next_rand2f(rng_pcg32& rng, const vec2f& a, const vec2f& b);

Next random float in [a.x,b.x)x[a.y,b.y).

#### Function next_rand3f()

    vec3f ygl::next_rand3f(rng_pcg32& rng);

Next random float3 in [0,1)x[0,1)x[0,1).

#### Function next_rand2f()

    vec3f ygl::next_rand2f(rng_pcg32& rng, const vec3f& a, const vec3f& b);

Next random float in [a.x,b.x)x[a.y,b.y)x[a.z,b.z).

#### Function next_rand1d()

    double ygl::next_rand1d(rng_pcg32& rng);

Next random double in [0, 1). Only 32 mantissa bits are filled, but
still better than float that uses 23.

#### Function rng_distance()

    int64_t ygl::rng_distance(const rng_pcg32& a, const rng_pcg32& b);

Distance between random number generators.

#### Function rng_shuffle()

    template<typename Iterator>
    void ygl::rng_shuffle(rng_pcg32& rng, Iterator begin, Iterator end);

Random shuffle of a sequence.

#### Function rng_shuffle()

    template<typename T>
    void ygl::rng_shuffle(rng_pcg32& rng, T* vals, int num);

Random shuffle of a sequence.

#### Function rng_shuffle()

    template<typename T>
    void ygl::rng_shuffle(rng_pcg32& rng, std::vector<T>& vals);

Random shuffle of a sequence.

#### Function operator==()

    bool ygl::operator==(const rng_pcg32& a, const rng_pcg32& b);

Random number generator equality.

#### Function operator!=()

    bool ygl::operator!=(const rng_pcg32& a, const rng_pcg32& b);

Random number generator inequality.

## Monte Carlo sampling

#### Function sample_hemisphere()

    vec3f ygl::sample_hemisphere(const vec2f& ruv);

Sample an hemispherical direction with uniform distribution.

#### Function sample_hemisphere_pdf()

    float ygl::sample_hemisphere_pdf(const vec3f& w);

Pdf for uniform hemispherical direction.

#### Function sample_sphere()

    vec3f ygl::sample_sphere(const vec2f& ruv);

Sample a spherical direction with uniform distribution.

#### Function sample_sphere_pdf()

    float ygl::sample_sphere_pdf(const vec3f& w);

Pdf for uniform spherical direction.

#### Function sample_hemisphere_cosine()

    vec3f ygl::sample_hemisphere_cosine(const vec2f& ruv);

Sample an hemispherical direction with cosine distribution.

#### Function sample_hemisphere_cosine_pdf()

    float ygl::sample_hemisphere_cosine_pdf(const vec3f& w);

Pdf for cosine hemispherical direction.

#### Function sample_hemisphere_cospower()

    vec3f ygl::sample_hemisphere_cospower(float n, const vec2f& ruv);

Sample an hemispherical direction with cosine power distribution.

#### Function sample_hemisphere_cospower_pdf()

    float ygl::sample_hemisphere_cospower_pdf(float n, const vec3f& w);

Pdf for cosine power hemispherical direction.

#### Function sample_disk()

    vec3f ygl::sample_disk(const vec2f& ruv);

Sample a point uniformly on a disk.

#### Function sample_disk_pdf()

    float ygl::sample_disk_pdf();

Pdf for uniform disk sampling.

#### Function sample_cylinder()

    vec3f ygl::sample_cylinder(const vec2f& ruv);

Sample a point uniformly on a cylinder, without caps.

#### Function sample_cylinder_pdf()

    float ygl::sample_cylinder_pdf();

Pdf for uniform cylinder sampling.

#### Function sample_triangle()

    vec2f ygl::sample_triangle(const vec2f& ruv);

Sample a point uniformly on a triangle.

#### Function sample_triangle()

    vec3f ygl::sample_triangle(const vec3f& v0, const vec3f& v1, const vec3f&
    v2, const vec2f& ruv);

Sample a point uniformly on a triangle.

#### Function sample_triangle_pdf()

    float ygl::sample_triangle_pdf(const vec3f& v0, const vec3f& v1, const
    vec3f& v2);

Pdf for uniform triangle sampling, i.e. triangle area.

#### Function sample_index()

    int ygl::sample_index(int size, float r);

Sample an index with uniform distribution.

#### Function sample_index_pdf()

    float ygl::sample_index_pdf(int size);

Pdf for uniform index sampling.

#### Function sample_discrete()

    int ygl::sample_discrete(const std::vector<float>& cdf, float r);

Sample a discrete distribution represented by its cdf.

#### Function sample_discrete_pdf()

    float ygl::sample_discrete_pdf(const std::vector<float>& cdf, int idx);

Pdf for uniform discrete distribution sampling.

## Hashing

#### Function cmjs_permute()

    uint32_t ygl::cmjs_permute(uint32_t i, uint32_t n, uint32_t key);

Computes the i-th term of a permutation of l values keyed by p. From
Correlated Multi-Jittered Sampling by Kensler @ Pixar.

#### Function cmjs_randfloat()

    float ygl::cmjs_randfloat(uint32_t i, uint32_t key);

Computes a float value by hashing i with a key p. From Correlated
Multi-Jittered Sampling by Kensler @ Pixar.

#### Function hash_uint32()

    uint32_t ygl::hash_uint32(uint64_t a);

32 bit integer hash. Public domain code.

#### Function hash_uint64()

    uint64_t ygl::hash_uint64(uint64_t a);

64 bit integer hash. Public domain code.

#### Function hash_uint64_32()

    uint32_t ygl::hash_uint64_32(uint64_t a);

64-to-32 bit integer hash. Public domain code.

#### Function hash_combine()

    size_t ygl::hash_combine(size_t a, size_t b);

Combines two 64 bit hashes as in boost::hash_combine.

## Perlin noise

#### Function perlin_noise()

    float ygl::perlin_noise(const vec3f& p, const vec3i& wrap=zero3i);

Compute the revised Pelin noise function. Wrap provides a wrapping
noise but must be power of two (wraps at 256 anyway). For octave based
noise, good values are obtained with octaves=6 (numerber of noise
calls), lacunarity=~2.0 (spacing between successive octaves: 2.0 for
warpping output), gain=0.5 (relative weighting applied to each
successive octave), offset=1.0 (used to invert the ridges).

#### Function perlin_ridge_noise()

    float ygl::perlin_ridge_noise(const vec3f& p, float lacunarity=2.0f, float
    gain=0.5f, float offset=1.0f, int octaves=6, const vec3i& wrap=zero3i);

Ridge noise function. See perlin_noise() for params.

#### Function perlin_fbm_noise()

    float ygl::perlin_fbm_noise(const vec3f& p, float lacunarity=2.0f, float
    gain=0.5f, int octaves=6, const vec3i& wrap=zero3i);

Fractal brownian motion noise. See perlin_noise() for params.

#### Function perlin_turbulence_noise()

    float ygl::perlin_turbulence_noise(const vec3f& p, float lacunarity=2.0f,
    float gain=0.5f, int octaves=6, const vec3i& wrap=zero3i);

Fractal turbulence noise. See perlin_noise() for params.

## Python-like iterators

#### Function range()

    range_generator ygl::range(int max);

Python-like range ierator.

#### Function range()

    range_generator ygl::range(int min, int max, int step=1);

Python-like range ierator.

#### Function enumerate()

    template<typename T>
    enumerate_generator<const T> ygl::enumerate(const std::vector<T>& vv);

Python-like enumerate.

#### Function enumerate()

    template<typename T>
    enumerate_generator<T> ygl::enumerate(std::vector<T>& vv);

Python-like enumerate.

## Container operations

#### Enum basic_variant_type

    enum struct basic_variant_type {

        none,
        boolean,
        integer,
        number,
        string,
    }

The type of the basic variant.

Members:
    - none: None variant to indicate uninitialized type.
    - boolean: Boolean variant.
    - integer: Integer variant.
    - number: Number variant.
    - string: String variant.

#### Function operator==()

    template<typename T>
    bool ygl::operator==(const optional<T>& a, const optional<T>& b);

Optional equality.

#### Function operator!=()

    template<typename T>
    bool ygl::operator!=(const optional<T>& a, const optional<T>& b);

Optional inequality.

#### Function make_optional()

    template<typename T>
    optional<T> ygl::make_optional(const T& v);

Creates an optional.

#### Function make_optional()

    template<typename T, typename...>
    optional<T> ygl::make_optional(Args& &... args);

Creates an optional.

#### Function operator<<()

    std::ostream& ygl::operator<<(std::ostream& os, const basic_variant& a);

Stream write.

#### Function operator>>()

    std::istream& ygl::operator>>(std::istream& is, basic_variant& a);

Stream read.

#### Function append()

    template<typename T>
    void ygl::append(std::vector<T>& v, const std::vector<T>& vv);

Append a vector to a vector.

#### Function join()

    template<typename T>
    std::vector<T> ygl::join(const std::vector<T>& a, const std::vector<T>&
    b);

Append two vectors.

#### Function get_key()

    template<typename K, typename V>
    K ygl::get_key(const std::vector<std::pair<K, V>>& kvs, const V& v);

Get a key from a list of key-value pairs and its value.

#### Function get_value()

    template<typename K, typename V>
    V ygl::get_value(const std::vector<std::pair<K, V>>& kvs, const K& k);

Get a value from a list of key-value pairs and its key.

#### Function find()

    template<typename T>
    int ygl::find(const std::vector<T>& v, const T& vv);

Find the position of a value in an array. Returns -1 if not found.
Wrapper for std::find().

#### Function upper_bound()

    template<typename T>
    int ygl::upper_bound(const std::vector<T>& v, const T& vv);

Find the first array value that is greater than the argument. Assumes
that the array is sorted. Wrapper for std::upper_bound().

#### Function lower_bound()

    template<typename T>
    int ygl::lower_bound(const std::vector<T>& v, const T& vv);

Find the first array value smaller that is greater or equal to the
argument. Assumes that the array is sorted. Wrapper for
std::lower_bound().

#### Function contains()

    template<typename T>
    bool ygl::contains(const std::vector<T>& v, const T& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V>
    bool ygl::contains(const std::map<K, V>& v, const K& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V>
    bool ygl::contains(const std::unordered_map<K, V>& v, const K& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V>
    bool ygl::contains(const std::set<K, V>& v, const K& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V>
    bool ygl::contains(const std::unordered_set<K, V>& v, const K& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V, typename K1>
    bool ygl::contains(const std::unordered_map<K, V>& v, const K1& vv);

Checks if a containers contains a value.

#### Function contains()

    template<typename K, typename V, typename K1>
    bool ygl::contains(const std::unordered_set<K, V>& v, const K1& vv);

Checks if a containers contains a value.

## Type support

#### Enum visit_var_type

    enum struct visit_var_type {

        value,
        name,
        path,
        object,
        reference,
        color,
        noneditable,
    }

Types of variable semantic for visit().

Members:
    - value: Generic value.
    - name: Name.
    - path: Path.
    - object: Object.
    - reference: Reference.
    - color: Color.
    - noneditable: Generic value not editable.

#### Function enum_names()

    template<typename T>
    const std::vector<std::pair<std::string, T>>& ygl::enum_names();

Names of enum values. Specialized by enums that support reflection.

#### Function enum_names()

    template<typename T>
    const std::vector<std::pair<std::string, T>>& ygl::enum_names(T v);

Names of enum values.

#### Function operator<<()

    template<typename T, typename std::enable_if<std::is_enum<T>::value,
    int>::type>
    std::ostream&  ygl::operator<<(std::ostream& os, const T& a);

Stream write.

#### Function operator>>()

    template<typename T, typename std::enable_if<std::is_enum<T>::value,
    int>::type>
    std::istream& ygl::operator>>(std::istream& is, T& a);

Stream read.

#### Function visit()

    template<typename T, typename Visitor>
    void ygl::visit(T& val, Visitor& &visitor);

Visit struct elements. Calls visitor(name,val.var,sem) for each
variable of a structure, where name is the name of the variable, var
is the variable and sem is one a visit_sem value. Implemented by
structures that support reflection.

#### Function visit()

    template<typename T, typename Visitor>
    void ygl::visit(T*& val, Visitor& &visitor);

Visit pointer elements.

## Geometry utilities

#### Function line_tangent()

    vec3f ygl::line_tangent(const vec3f& v0, const vec3f& v1);

Line tangent.

#### Function line_length()

    float ygl::line_length(const vec3f& v0, const vec3f& v1);

Line length.

#### Function triangle_normal()

    vec3f ygl::triangle_normal(const vec3f& v0, const vec3f& v1, const vec3f&
    v2);

Triangle normal.

#### Function triangle_area()

    float ygl::triangle_area(const vec3f& v0, const vec3f& v1, const vec3f&
    v2);

Triangle area.

#### Function quad_normal()

    vec3f ygl::quad_normal(const vec3f& v0, const vec3f& v1, const vec3f& v2,
    const vec3f& v3);

Quad normal.

#### Function quad_area()

    float ygl::quad_area(const vec3f& v0, const vec3f& v1, const vec3f& v2,
    const vec3f& v3);

Quad area.

#### Function tetrahedron_volume()

    float ygl::tetrahedron_volume(const vec3f& v0, const vec3f& v1, const
    vec3f& v2, const vec3f& v3);

tetrahedron volume

#### Function triangle_tangents_fromuv()

    std::pair<vec3f, vec3f> ygl::triangle_tangents_fromuv(const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2);

Triangle tangent and bitangent from uv (not othornormalized with
themselfves not the normal). Follows the definition in
http://www.terathon.com/code/tangent.html and
https://gist.github.com/aras-p/2843984.

#### Function interpolate_point()

    template<typename T>
    T ygl::interpolate_point(const std::vector<T>& vals, int p);

Copies of point value. Here only for completeness.

#### Function interpolate_line()

    template<typename T, typename T1>
    T ygl::interpolate_line(const T& v0, const T& v1, const T1 u);

Interpolates values over a line parametrized from a to b by u. Same as
lerp.

#### Function interpolate_line()

    template<typename T, typename T1>
    T ygl::interpolate_line(const std::vector<T>& vals, const vec2i& l, T1 u);

Interpolates values over lines parametrized from a to b by u. Same as
lerp.

#### Function interpolate_triangle()

    template<typename T, typename T1>
    T ygl::interpolate_triangle(const T& v0, const T& v1, const T& v2, const
    vec<T1, 2>& uv);

Interpolates values over a triangle parametrized by u and v along the
(v1-v0) and (v2-v0) directions. Same as barycentric interpolation.

#### Function interpolate_triangle()

    template<typename T, typename T1>
    T ygl::interpolate_triangle(const std::vector<T>& vals, const vec3i& t,
    const vec<T1, 2>& uv);

Interpolates values over triangles parametrized by u and v along the
(v1-v0) and (v2-v0) directions. Same as barycentric interpolation.

#### Function interpolate_triangle()

    template<typename T, typename T1>
    T ygl::interpolate_triangle(const T& v0, const T& v1, const T& v2, const
    T& v3, const vec<T1, 2>& uv);

Interpolates values over a quad parametrized by u and v along the
(v1-v0) and (v2-v1) directions. Same as bilear interpolation.

#### Function interpolate_quad()

    template<typename T, typename T1>
    T ygl::interpolate_quad(const std::vector<T>& vals, const vec4i& t, const
    vec<T1, 2>& uv);

Interpolates values over quads parametrized by u and v along the
(v1-v0) and (v2-v1) direction. Same as bilear interpolation.

#### Function eval_bernstein()

    template<typename T>
    T ygl::eval_bernstein(T u, int i, int degree);

Evaluates the i-th Bernstein polynomial of degree degree at u.

#### Function eval_bernstein_derivative()

    template<typename T>
    T ygl::eval_bernstein_derivative(T u, int i, int degree);

Evaluates the derivative of i-th Bernstein polynomial of degree degree
at u.

#### Function interpolate_bezier()

    template<typename T, typename T1>
    T ygl::interpolate_bezier(const T& v0, const T& v1, const T& v2, const T&
    v3, T1 u);

Interpolates values along a cubic Bezier segment parametrized by u.

#### Function interpolate_bezier()

    template<typename T, typename T1>
    T ygl::interpolate_bezier(const std::vector<T>& vals, const vec4i& b, T1
    u);

Interpolates values along cubic Bezier segments parametrized by u.

#### Function interpolate_bezier_derivative()

    template<typename T, typename T1>
    T ygl::interpolate_bezier_derivative(const T& v0, const T& v1, const T&
    v2, const T& v3, T1 u);

Computes the derivative of a cubic Bezier segment parametrized by u.

#### Function interpolate_bezier_derivative()

    template<typename T, typename T1>
    T ygl::interpolate_bezier_derivative(const std::vector<T>& vals, const
    vec4i& b, T1 u);

Computes the derivative of a cubic Bezier segments parametrized by u.

## Animation utilities

#### Function eval_keyframed_step()

    template<typename T>
    T ygl::eval_keyframed_step(const std::vector<float>& times, const
    std::vector<T>& vals, float time);

Evalautes a keyframed value using step interpolation.

#### Function eval_keyframed_lerp()

    template<typename T>
    T ygl::eval_keyframed_lerp(const T& a, const T& b, float t);



#### Function eval_keyframed_lerp()

    template<typename T>
    quat<T, 4> ygl::eval_keyframed_lerp(const quat<T, 4>& a, const quat<T, 4>&
    b, float t);



#### Function eval_keyframed_linear()

    template<typename T>
    T ygl::eval_keyframed_linear(const std::vector<float>& times, const
    std::vector<T>& vals, float time);

Evalautes a keyframed value using linear interpolation.

#### Function eval_keyframed_cubic()

    template<typename T>
    T ygl::eval_keyframed_cubic(const T& a, const T& b, const T& c, const T&
    d, float t);



#### Function eval_keyframed_cubic()

    template<typename T>
    quat<T, 4> ygl::eval_keyframed_cubic(const quat<T, 4>& a, const quat<T,
    4>& b, const quat<T, 4>& c, const quat<T, 4>& d, float t);



#### Function eval_keyframed_bezier()

    template<typename T>
    T ygl::eval_keyframed_bezier(const std::vector<float>& times, const
    std::vector<T>& vals, float time);

Evalautes a keyframed value using Bezier interpolation.

## Shape utilities

#### Function compute_normals()

    std::vector<vec3f> ygl::compute_normals(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, bool weighted=true);

Compute per-vertex normals/tangents for lines, triangles and quads
with positions pos. Weighted indicated whether the normals/tangents
are weighted by length/area.

#### Function compute_tangent_frames()

    std::vector<vec4f> ygl::compute_tangent_frames(const std::vector<vec3i>&
    triangles, const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, bool weighted=true);

Compute per-vertex tangent frames for triangle meshes. Tangent space
is defined by a four component vector. The first three components are
the tangent with respect to the u texcoord. The fourth component is
the sign of the tangent wrt the v texcoord. Tangent frame is useful in
normal mapping.

#### Function compute_skinning()

    void ygl::compute_skinning(const std::vector<vec3f>& pos, const
    std::vector<vec3f>& norm, const std::vector<vec4f>& weights, const
    std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

Apply skinning to vertex position and normals.

#### Function compute_skinning()

    void ygl::compute_skinning(const std::vector<vec3f>& pos, const
    std::vector<vec3f>& norm, const std::vector<vec4f>& weights, const
    std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

Apply skinning.

#### Function compute_matrix_skinning()

    void ygl::compute_matrix_skinning(const std::vector<vec3f>& pos, const
    std::vector<vec3f>& norm, const std::vector<vec4f>& weights, const
    std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

Apply skinning as specified in Khronos glTF.

#### Function get_edges()

    std::vector<vec2i> ygl::get_edges(const std::vector<vec2i>& lines, const
    std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

Create an array of edges.

#### Function get_boundary_edges()

    std::vector<vec2i> ygl::get_boundary_edges(const std::vector<vec2i>&
    lines, const std::vector<vec3i>& triangles, const std::vector<vec4i>&
    quads);

Create an array of boundary edges. Lines are always considered
boundaries.

#### Function get_verts()

    std::vector<int> ygl::get_verts(const std::vector<vec2i>& lines, const
    std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

Get a list of all unique vertices.

#### Function get_boundary_verts()

    std::vector<int> ygl::get_boundary_verts(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

Create an array of boundary vertices. Lines are always considered
boundaries.

#### Function convert_quads_to_triangles()

    std::vector<vec3i> ygl::convert_quads_to_triangles(const
    std::vector<vec4i>& quads);

Convert quads to triangles.

#### Function convert_quads_to_triangles()

    std::vector<vec3i> ygl::convert_quads_to_triangles(const
    std::vector<vec4i>& quads, int row_length);

Convert quads to triangles with a diamond-like topology. Quads have to
be consecutive one row after another.

#### Function convert_bezier_to_lines()

    std::vector<vec2i> ygl::convert_bezier_to_lines(const std::vector<vec4i>&
    beziers);

Convert beziers to lines using 3 lines for each bezier.

#### Function convert_face_varying()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::convert_face_varying(const std::vector<vec4i>&
    quads_pos, const std::vector<vec4i>& quads_norm, const std::vector<vec4i>&
    quads_texcoord, const std::vector<vec3f>& pos, const std::vector<vec3f>&
    norm, const std::vector<vec2f>& texcoord);

Convert face-varying data to single primitives. Returns the quads
indices and filled vectors for pos, norm and texcoord.

#### Function subdivide_elems_linear()

    std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>,
    std::vector<vec2i>, std::vector<vec4i>> ygl::subdivide_elems_linear(const
    std::vector<vec2i>& lines, const std::vector<vec3i>& triangles, const
    std::vector<vec4i>& quads, int nverts);

Tesselate lines, triangles and quads by splitting edges and faces for
quads. Returns the tesselated elements and edge and face array for
vertex calculations.

#### Function subdivide_vert_linear()

    template<typename T>
    std::vector<T> ygl::subdivide_vert_linear(const std::vector<T>& vert,
    const std::vector<vec2i>& edges, const std::vector<vec4i>& faces, bool
    normalized=false);

Subdivide vertex properties for edges and faces. This is implemented
for vecs and floats.

#### Function subdivide_vert_catmullclark()

    template<typename T>
    std::vector<T> ygl::subdivide_vert_catmullclark(const std::vector<vec4i>&
    quads, const std::vector<T>& vert, const std::vector<vec2i>&
    crease_tlines, const std::vector<int>& crease_tpoints, bool
    normalized=false);

Performs the smoothing step of Catmull-Clark. Start with a tesselate
quad mesh obtained with subdivide_elems_linear() and
subdivide_vert_linear(). To handle open meshes with boundary, get the
boundary from make_boundary_edge() and pass it as crease_lines. To fix
the boundary entirely, just get the boundary vertices and pass it as
creases. This is implemented for vecs and floats.

#### Function subdivide_bezier_recursive()

    std::tuple<std::vector<vec4i>, std::vector<int>, std::vector<vec4i>>
    ygl::subdivide_bezier_recursive(const std::vector<vec4i>& beziers, int
    nverts);

Subdivide Bezier segments recursively by splitting each segment into
two in the middle. Returns the tesselated elements and dictionaries
for vertex calculations.

#### Function subdivide_vert_bezier()

    template<typename T>
    std::vector<T> ygl::subdivide_vert_bezier(const std::vector<T>& vert,
    const std::vector<int>& verts, const std::vector<vec4i>& segments, bool
    normalized=false);

Subdivide vertex properties for Bezier subdivision. This is implemente
for vecs and floats.

#### Function make_uvquads()

    std::tuple<std::vector<vec4i>, std::vector<vec2f>> ygl::make_uvquads(int
    usteps, int vsteps, bool uwrap=false, bool vwrap=false, bool vpole0=false,
    bool vpole1=false);

Generate a rectangular grid of usteps x vsteps uv values for
parametric surface generation. Values cam wrap and have poles.

#### Function make_uvlines()

    std::tuple<std::vector<vec2i>, std::vector<vec2f>> ygl::make_uvlines(int
    num, int usteps);

Generate parametric num lines of usteps segments.

#### Function make_uvpoints()

    std::tuple<std::vector<int>, std::vector<vec2f>> ygl::make_uvpoints(int
    num);

Generate a parametric point set. Mostly here for completeness.

#### Function merge_elems()

    std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>>
    ygl::merge_elems(int nverts, const std::vector<vec2i>& lines1, const
    std::vector<vec3i>& triangles1, const std::vector<vec4i>& quads1, const
    std::vector<vec2i>& lines2, const std::vector<vec3i>& triangles2, const
    std::vector<vec4i>& quads2);

Merge elements between shapes. The elements are merged by increasing
the array size of the second array by the number of vertices of the
first. Vertex data can then be concatenated successfully.

#### Function facet_elems()

    std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>,
    std::vector<int>> ygl::facet_elems(const std::vector<vec2i>& lines, const
    std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

Unshare shape data by duplicating all vertex data for each element,
giving a faceted look. Note that faceted tangents are not computed.
Returns the indices to copy from.

#### Function facet_vert()

    template<typename T>
    std::vector<T> ygl::facet_vert(const std::vector<T>& vert, const
    std::vector<int>& vmap);

Unshare vertices for faceting. This is implemented for vec and float
types.

## Shape sampling

#### Function sample_points()

    int ygl::sample_points(int npoints, float re);

Pick a point.

#### Function sample_points_cdf()

    std::vector<float> ygl::sample_points_cdf(int npoints);

Compute a distribution for sampling points uniformly.

#### Function sample_points()

    int ygl::sample_points(const std::vector<float>& cdf, float re);

Pick a point uniformly.

#### Function sample_lines_cdf()

    std::vector<float> ygl::sample_lines_cdf(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos);

Compute a distribution for sampling lines uniformly.

#### Function sample_lines()

    std::pair<int, float> ygl::sample_lines(const std::vector<float>& cdf,
    float re, float ru);

Pick a point on lines uniformly.

#### Function sample_triangles_cdf()

    std::vector<float> ygl::sample_triangles_cdf(const std::vector<vec3i>&
    triangles, const std::vector<vec3f>& pos);

Compute a distribution for sampling triangle meshes uniformly.

#### Function sample_triangles()

    std::pair<int, vec2f> ygl::sample_triangles(const std::vector<float>& cdf,
    float re, const vec2f& ruv);

Pick a point on a triangle mesh uniformly.

#### Function sample_quads_cdf()

    std::vector<float> ygl::sample_quads_cdf(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos);

Compute a distribution for sampling quad meshes uniformly.

#### Function sample_quads()

    std::pair<int, vec2f> ygl::sample_quads(const std::vector<float>& cdf,
    float re, const vec2f& ruv);

Pick a point on a quad mesh uniformly.

#### Function sample_triangles_points()

    std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
    ygl::sample_triangles_points(const std::vector<vec3i>& triangles, const
    std::vector<vec3f>& pos, const std::vector<vec3f>& norm, const
    std::vector<vec2f>& texcoord, int npoints, uint64_t seed=0);

Samples a set of points over a triangle mesh uniformly. Returns pos,
norm and tecoord of the sampled points.

## Example shapes

#### Function make_sphere()

    std::tuple<std::vector<vec3i>, std::vector<vec3f>> ygl::make_sphere(int
    tesselation);

Make a sphere. Returns quads, pos.

#### Function make_geodesicsphere()

    std::tuple<std::vector<vec3i>, std::vector<vec3f>>
    ygl::make_geodesicsphere(int tesselation);

Make a geodesic sphere. Returns quads, pos.

#### Function make_cube()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>> ygl::make_cube(int
    tesselation);

Make a cube with unique vertices. This is watertight but has no
texture coordinates or normals. Returns quads, pos.

#### Function make_uvsphere()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvsphere(int tesselation, bool
    flipped=false);

Make a sphere. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvhemisphere()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvhemisphere(int tesselation, bool
    flipped=false);

Make a sphere. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvquad()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvquad(int tesselation);

Make a quad. Returns quads, pos, norm, texcoord.

#### Function make_fvsphere()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec4i>,
    std::vector<vec3f>, std::vector<vec4i>, std::vector<vec2f>>
    ygl::make_fvsphere(int tesselation);

Make a facevarying sphere with unique vertices but different texture
coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).

#### Function make_fvcube()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec4i>,
    std::vector<vec3f>, std::vector<vec4i>, std::vector<vec2f>>
    ygl::make_fvcube(int tesselation);

Make a facevarying cube with unique vertices but different texture
coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).

#### Function make_suzanne()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>> ygl::make_suzanne(int
    tesselation);

Make a suzanne monkey model for testing. Note that some quads are
degenerate. Returns quads, pos.

#### Function make_uvcube()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvcube(int tesselation);

Make a cube. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvspherecube()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvspherecube(int tesselation);

Make a sphere from a cube. This is not watertight. Returns quads, pos,
norm, texcoord.

#### Function make_uvspherizedcube()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvspherizedcube(int tesselation, float
    radius);

Make a cube than stretch it towards a sphere. This is not watertight.
Returns quads, pos, norm, texcoord.

#### Function make_uvflipcapsphere()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvflipcapsphere(int tesselation, float z,
    bool flipped=false);

Make a sphere with caps flipped. This is not watertight. Returns
quads, pos, norm, texcoord.

#### Function make_uvcutsphere()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvcutsphere(int tesselation, float z, bool
    flipped=false);

Make a cutout sphere. This is not watertight. Returns quads, pos,
norm, texcoord.

#### Function make_uvseashell()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>> ygl::make_uvseashell(int tesselation, const
    make_seashell_params& params);

Make a seashell. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_bezier_circle()

    std::tuple<std::vector<vec4i>, std::vector<vec3f>>
    ygl::make_bezier_circle();

Make a bezier circle. Returns bezier, pos.

#### Function make_hair()

    std::tuple<std::vector<vec2i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>, std::vector<float>> ygl::make_hair(int num, int
    tesselation, const std::vector<vec3i>& striangles, const
    std::vector<vec4i>& squads, const std::vector<vec3f>& spos, const
    std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord, const
    make_hair_params& params);

Make a hair ball around a shape. Returns lines, pos, norm, texcoord,
radius.

## Example shape type support

#### Function visit()

    template<typename Visitor>
    void ygl::visit(make_hair_params& val, Visitor& &visitor);

Visit struct elements.

## Image containers

#### Typedef image4f

    using ygl::image4f = image<vec4f>;

HDR image.

#### Typedef image4b

    using ygl::image4b = image<vec4b>;

LDR image.

#### Function begin()

    template<typename T>
    T* ygl::begin(image<T>& a);

Pixel iteration.

#### Function begin()

    template<typename T>
    const T* ygl::begin(const image<T>& a);

Pixel iteration.

#### Function end()

    template<typename T>
    T* ygl::end(image<T>& a);

Pixel iteration.

#### Function end()

    template<typename T>
    const T* ygl::end(const image<T>& a);

Pixel iteration.

#### Function data()

    template<typename T>
    T* ygl::data(image<T>& a);

Pixel access.

#### Function data()

    template<typename T>
    const T* ygl::data(const image<T>& a);

Pixel access.

#### Function size()

    template<typename T>
    int ygl::size(image<T>& a);

Number of pixels.

#### Function empty()

    template<typename T>
    bool ygl::empty(image<T>& a);

Chech if an image is empty.

#### Function make_image()

    template<typename T>
    image<T> ygl::make_image(int width, int height, T* vals);

Create an image with values stored in an array in scanliine order.

#### Function make_image()

    template<typename T>
    image<vec<T, 4>> ygl::make_image(int w, int h, int nc, const T* vals,
    const vec<T, 4>& def);

Create a 4 channel image with the given number of channels.

## Image operations

#### Function srgb_to_linear()

    vec4f ygl::srgb_to_linear(const vec4b& srgb);

Approximate conversion from srgb.

#### Function linear_to_srgb()

    vec4b ygl::linear_to_srgb(const vec4f& lin);

Approximate conversion to srgb.

#### Function xyz_to_xyY()

    vec3f ygl::xyz_to_xyY(const vec3f& xyz);

Convert between CIE XYZ and xyY.

#### Function xyY_to_xyz()

    vec3f ygl::xyY_to_xyz(const vec3f& xyY);

Convert between CIE XYZ and xyY.

#### Function xyz_to_rgb()

    vec3f ygl::xyz_to_rgb(const vec3f& xyz);

Convert between CIE XYZ and RGB.

#### Function rgb_to_xyz()

    vec3f ygl::rgb_to_xyz(const vec3f& rgb);

Convert between CIE XYZ and RGB.

#### Function tonemap_image()

    image4b ygl::tonemap_image(const image4f& hdr, float exposure, float
    gamma, bool filmic=false);

Tone mapping HDR to LDR images.

#### Function image_over()

    void ygl::image_over(vec4f* img, int width, int height, int nlayers,
    vec4f* *layers);

Image over operator.

#### Function image_over()

    void ygl::image_over(vec4b* img, int width, int height, int nlayers,
    vec4b* *layers);

Image over operator.

#### Function hsv_to_rgb()

    vec4b ygl::hsv_to_rgb(const vec4b& hsv);

Converts HSV to RGB.

## Example images

#### Function make_grid_image()

    image4b ygl::make_grid_image(int width, int height, int tile=64, const
    vec4b& c0={90, 90, 90, 255}, const vec4b& c1={128, 128, 128, 255});

Make a grid image.

#### Function make_checker_image()

    image4b ygl::make_checker_image(int width, int height, int tile=64, const
    vec4b& c0={90, 90, 90, 255}, const vec4b& c1={128, 128, 128, 255});

Make a checkerboard image.

#### Function make_bumpdimple_image()

    image4b ygl::make_bumpdimple_image(int width, int height, int tile=64);

Make an image with bumps and dimples.

#### Function make_ramp_image()

    image4b ygl::make_ramp_image(int width, int height, const vec4b& c0, const
    vec4b& c1, bool srgb=false);

Make a uv colored grid.

#### Function make_gammaramp_image()

    image4b ygl::make_gammaramp_image(int width, int height);

Make a gamma ramp image.

#### Function make_gammaramp_imagef()

    image4f ygl::make_gammaramp_imagef(int width, int height);

Make a gamma ramp image.

#### Function make_uv_image()

    image4b ygl::make_uv_image(int width, int height);

Make an image color with red/green in the [0,1] range. Helpful to
visualize uv texture coordinate application.

#### Function make_uvgrid_image()

    image4b ygl::make_uvgrid_image(int width, int height, int tile=64, bool
    colored=true);

Make a uv colored grid.

#### Function make_recuvgrid_image()

    image4b ygl::make_recuvgrid_image(int width, int height, int tile=64, bool
    colored=true);

Make a uv recusive colored grid.

#### Function bump_to_normal_map()

    image4b ygl::bump_to_normal_map(const image4b& img, float scale=1);

Comvert a bump map to a normal map.

#### Function make_sunsky_image()

    image4f ygl::make_sunsky_image(int res, float thetaSun, float turbidity=3,
    bool has_sun=false, bool has_ground=true);

Make a sunsky HDR model with sun at theta elevation in [0,pi/2],
turbidity in [1.7,10] with or without sun.

#### Function make_noise_image()

    image4b ygl::make_noise_image(int resx, int resy, float scale=1, bool
    wrap=true);

Make a noise image. Wrap works only if both resx and resy are powers
of two.

#### Function make_fbm_image()

    image4b ygl::make_fbm_image(int resx, int resy, float scale=1, float
    lacunarity=2, float gain=0.5f, int octaves=6, bool wrap=true);

Make a noise image. Wrap works only if both resx and resy are powers
of two.

#### Function make_ridge_image()

    image4b ygl::make_ridge_image(int resx, int resy, float scale=1, float
    lacunarity=2, float gain=0.5f, float offset=1.0f, int octaves=6, bool
    wrap=true);

Make a noise image. Wrap works only if both resx and resy are powers
of two.

#### Function make_turbulence_image()

    image4b ygl::make_turbulence_image(int resx, int resy, float scale=1,
    float lacunarity=2, float gain=0.5f, int octaves=6, bool wrap=true);

Make a noise image. Wrap works only if both resx and resy are powers
of two.

## Image loading and saving

#### Enum resize_filter

    enum struct resize_filter {

        def,
        box,
        triangle,
        cubic_spline,
        catmull_rom,
        mitchell,
    }

Filter type for resizing.

Members:
    - def: default
    - box: box filter
    - triangle: triangle filter
    - cubic_spline: cubic spline
    - catmull_rom: Catmull-Rom interpolating sline.
    - mitchell: Mitchel-Netrevalli filter with B=1/3, C=1/3.

#### Enum resize_edge

    enum struct resize_edge {

        def,
        clamp,
        reflect,
        wrap,
        zero,
    }

Edge mode for resizing.

Members:
    - def: default
    - clamp: clamp
    - reflect: reflect
    - wrap: wrap
    - zero: zero

#### Function is_hdr_filename()

    bool ygl::is_hdr_filename(const std::string& filename);

Check if an image is HDR based on filename.

#### Function load_image4b()

    image4b ygl::load_image4b(const std::string& filename);

Loads a 4 channel ldr image.

#### Function load_image4f()

    image4f ygl::load_image4f(const std::string& filename);

Loads a 4 channel hdr image.

#### Function save_image4b()

    bool ygl::save_image4b(const std::string& filename, const image4b& img);

Saves a 4 channel ldr image.

#### Function save_image4f()

    bool ygl::save_image4f(const std::string& filename, const image4f& img);

Saves a 4 channel hdr image.

#### Function load_imagef()

    std::vector<float> ygl::load_imagef(const std::string& filename, int&
    width, int& height, int& ncomp);

Loads an image with variable number of channels.

#### Function load_image()

    std::vector<byte> ygl::load_image(const std::string& filename, int& width,
    int& height, int& ncomp);

Loads an image with variable number of channels.

#### Function load_imagef_from_memory()

    std::vector<float> ygl::load_imagef_from_memory(const std::string&
    filename, const byte* data, int length, int& width, int& height, int&
    ncomp);

Loads an image from memory with variable number of channels.

#### Function load_image_from_memory()

    std::vector<byte> ygl::load_image_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp);

Loads an image from memory with variable number of channels.

#### Function save_imagef()

    bool ygl::save_imagef(const std::string& filename, int width, int height,
    int ncomp, const float* hdr);

Saves an image with variable number of channels.

#### Function save_image()

    bool ygl::save_image(const std::string& filename, int width, int height,
    int ncomp, const byte* ldr);

Saves an image with variable number of channels.

#### Function save_image()

    bool ygl::save_image(const std::string& filename, const image4f& hdr,
    float exposure, float gamma, bool filmic=false);

Save a 4 channel HDR or LDR image with tonemapping based on filename.

#### Function resize_image()

    void ygl::resize_image(const image4f& img, image4f& res_img, resize_filter
    filter=resize_filter::def, resize_edge edge=resize_edge::def, bool
    premultiplied_alpha=false);

Resize an image.

#### Function resize_image()

    void ygl::resize_image(const image4b& img, image4b& res_img, resize_filter
    filter=resize_filter::def, resize_edge edge=resize_edge::def, bool
    premultiplied_alpha=false);

Resize an image.

## Ray-primitive intersection

#### Function intersect_point()

    bool ygl::intersect_point(const ray3f& ray, const vec3f& p, float r,
    float& ray_t);

Intersect a ray with a point (approximate). Based on
http://geomalgorithms.com/a02-lines.html.

#### Function intersect_line()

    bool ygl::intersect_line(const ray3f& ray, const vec3f& v0, const vec3f&
    v1, float r0, float r1, float& ray_t, vec2f& euv);

Intersect a ray with a line (approximate). Based on
http://geomalgorithms.com/a05-intersect-1.html and
http://geomalgorithms.com/a07-distance.html#
dist3D_Segment_to_Segment.

#### Function intersect_triangle()

    bool ygl::intersect_triangle(const ray3f& ray, const vec3f& v0, const
    vec3f& v1, const vec3f& v2, float& ray_t, vec2f& euv);

Intersect a ray with a triangle.

#### Function intersect_quad()

    bool ygl::intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f&
    v1, const vec3f& v2, const vec3f& v3, float& ray_t, vec2f& euv);

Intersect a ray with a quad represented as two triangles (0,1,3) and
(2,3,1), with the uv coordinates of the second triangle corrected by u
= 1-u' and v = 1-v' to produce a quad parametrization where u and v go
from 0 to 1. This is equivalent to Intel's Embree.

#### Function intersect_bbox()

    bool ygl::intersect_bbox(const ray3f& ray, const bbox3f& bbox);

Intersect a ray with a axis-aligned bounding box.

#### Function intersect_bbox()

    bool ygl::intersect_bbox(const ray3f& ray, const vec3f& ray_dinv, const
    vec3i& ray_dsign, const bbox3f& bbox);

Intersect a ray with a axis-aligned bounding box, implemented as
"Robust BVH Ray Traversal" by T. Ize published at
http://jcgt.org/published/0002/02/02/paper.pdf.

## Point-primitive overlap

#### Function overlap_point()

    bool ygl::overlap_point(const vec3f& pos, float dist_max, const vec3f& v0,
    float r0, float& dist);

Check if a point overlaps a position within a max distance.

#### Function closestuv_line()

    float ygl::closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f&
    v1);

Find closest line point to a position.

#### Function overlap_line()

    bool ygl::overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv);

Check if a line overlaps a position within a max distance.

#### Function closestuv_triangle()

    vec2f ygl::closestuv_triangle(const vec3f& pos, const vec3f& v0, const
    vec3f& v1, const vec3f& v2);

Find closest triangle point to a position.

#### Function overlap_triangle()

    bool ygl::overlap_triangle(const vec3f& pos, float dist_max, const vec3f&
    v0, const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float&
    dist, vec2f& euv);

Check if a triangle overlaps a position within a max distance.

#### Function overlap_quad()

    bool ygl::overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& euv);

Check if a quad overlaps a position within a max distance.

#### Function overlap_bbox()

    bool ygl::overlap_bbox(const vec3f& pos, float dist_max, const bbox3f&
    bbox);

Check if a bouning box overlaps a position within a max distance.

#### Function overlap_bbox()

    bool ygl::overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

Check if two bouning boxes overlap.

## Bounding volume hierarchy

#### Enum bvh_node_type

    enum struct bvh_node_type {

        internal,
        point,
        line,
        triangle,
        quad,
        vertex,
        instance,
    }

Type of BVH node.

Members:
    - internal: Internal.
    - point: Points.
    - line: Lines.
    - triangle: Triangles.
    - quad: Quads.
    - vertex: Vertices.
    - instance: Instances.

#### Function make_bvh()

    bvh_tree* ygl::make_bvh(const std::vector<int>& points, const
    std::vector<vec2i>& lines, const std::vector<vec3i>& triangles, const
    std::vector<vec4i>& quads, const std::vector<vec3f>& pos, const
    std::vector<float>& radius, float def_radius, bool equalsize);

Build a shape BVH from a set of primitives.

#### Function make_bvh()

    bvh_tree* ygl::make_bvh(const std::vector<bvh_instance>& instances, const
    std::vector<bvh_tree*>& shape_bvhs, bool own_shape_bvhs, bool equal_size);

Build a scene BVH from a set of shape instances.

#### Function get_shape_bvhs()

    const std::vector<bvh_tree*>& ygl::get_shape_bvhs(const bvh_tree* bvh);

Grab the shape BVHs.

#### Function refit_bvh()

    void ygl::refit_bvh(bvh_tree* bvh, const std::vector<vec3f>& pos, const
    std::vector<float>& radius, float def_radius);

Update the node bounds for a shape bvh.

#### Function refit_bvh()

    void ygl::refit_bvh(bvh_tree* bvh, const std::vector<frame3f>& frames,
    const std::vector<frame3f>& frames_inv);

Update the node bounds for a scene bvh.

#### Function intersect_bvh()

    bool ygl::intersect_bvh(const bvh_tree* bvh, const ray3f& ray, bool
    find_any, float& ray_t, int& iid, int& sid, int& eid, vec2f& euv);

Intersect ray with a bvh returning either the first or any
intersection depending on find_any. Returns the ray distance ray_t,
the instance id iid, the shape id sid, the shape element index eid and
the shape barycentric coordinates euv.

#### Function overlap_bvh()

    bool ygl::overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float
    max_dist, bool find_any, float& dist, int& iid, int& sid, int& eid, vec2f&
    euv);

Find a shape element that overlaps a point within a given distance
max_dist, returning either the closest or any overlap depending on
find_any. Returns the point distance dist, the instance id iid, the
shape id sid, the shape element index eid and the shape barycentric
coordinates euv.

#### Function intersect_bvh()

    intersection_point ygl::intersect_bvh(const bvh_tree* bvh, const ray3f&
    ray, bool early_exit);

Intersect a ray with a bvh (convenience wrapper).

#### Function overlap_bvh()

    intersection_point ygl::overlap_bvh(const bvh_tree* bvh, const vec3f& pos,
    float max_dist, bool early_exit);

Finds the closest element with a bvh (convenience wrapper).

## Simple scene

#### Enum material_type

    enum struct material_type {

        specular_roughness,
        metallic_roughness,
        specular_glossiness,
    }

Material type.

Members:
    - specular_roughness: Microfacet material type (OBJ).
    - metallic_roughness: Base and metallic material (metallic-roughness in glTF).
    - specular_glossiness: Diffuse and specular material (specular-glossness in glTF).

#### Enum keyframe_type

    enum struct keyframe_type {

        linear,
        step,
        catmull_rom,
        bezier,
    }

Keyframe type.

Members:
    - linear: Linear interpolation.
    - step: Step function.
    - catmull_rom: Catmull-Rom interpolation.
    - bezier: Cubic Bezier interpolation.

#### Function eval_pos()

    vec3f ygl::eval_pos(const shape* shp, int eid, const vec2f& euv);

Shape position interpolated using barycentric coordinates.

#### Function eval_norm()

    vec3f ygl::eval_norm(const shape* shp, int eid, const vec2f& euv);

Shape normal interpolated using barycentric coordinates.

#### Function eval_texcoord()

    vec2f ygl::eval_texcoord(const shape* shp, int eid, const vec2f& euv);

Shape texcoord interpolated using barycentric coordinates.

#### Function eval_color()

    vec4f ygl::eval_color(const shape* shp, int eid, const vec2f& euv);

Shape color interpolated using barycentric coordinates.

#### Function eval_radius()

    float ygl::eval_radius(const shape* shp, int eid, const vec2f& euv);

Shape radius interpolated using barycentric coordinates.

#### Function eval_tangsp()

    vec4f ygl::eval_tangsp(const shape* shp, int eid, const vec2f& euv);

Shape tangent space interpolated using barycentric coordinates.

#### Function eval_pos()

    vec3f ygl::eval_pos(const instance* ist, int sid, int eid, const vec2f&
    euv);

Instance position interpolated using barycentric coordinates.

#### Function eval_norm()

    vec3f ygl::eval_norm(const instance* ist, int sid, int eid, const vec2f&
    euv);

Instance normal interpolated using barycentric coordinates.

#### Function eval_texture()

    vec4f ygl::eval_texture(const texture* txt, const texture_info& info,
    const vec2f& texcoord, bool srgb=true, const vec4f& def={1, 1, 1, 1});

Evaluate a texture.

#### Function eval_texture()

    vec4f ygl::eval_texture(const texture* txt, const optional<texture_info>&
    info, const vec2f& texcoord, bool srgb=true, const vec4f& def={1, 1, 1,
    1});

Evaluate a texture.

#### Function eval_camera_ray()

    ray3f ygl::eval_camera_ray(const camera* cam, const vec2f& uv, const
    vec2f& luv);

Generates a ray from a camera for image plane coordinate uv and the
lens coordinates luv.

#### Function eval_camera_ray()

    ray3f ygl::eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv);

Generates a ray from a camera for pixel coordinates ij, the resolution
res, the sub-pixel coordinates puv and the lens coordinates luv and
the image resolution res.

#### Function sync_camera_aspect()

    void ygl::sync_camera_aspect(camera* cam, int& width, int& height);

Synchronizes a camera aspect with image width and height. Set image
values any one is 0 or less. Set camera aspect otherwise.

#### Function find_named_elem()

    template<typename T>
    T* ygl::find_named_elem(const std::vector<T*>& elems, const std::string&
    name);

Finds an element by name.

#### Function add_named_elem()

    template<typename T>
    T* ygl::add_named_elem(std::vector<T*>& elems, const std::string& name);



#### Function subdivide_shape_once()

    void ygl::subdivide_shape_once(shape* shp, bool subdiv=false);

Subdivides shape elements. Apply subdivision surface rules if
subdivide is true.

#### Function facet_shape()

    void ygl::facet_shape(shape* shp, bool recompute_normals=true);

Facet a shape. Supports only non-facevarying shapes.

#### Function tesselate_shape()

    void ygl::tesselate_shape(shape* shp, bool subdivide, bool
    facevarying_to_sharedvertex, bool quads_to_triangles, bool
    bezier_to_lines);

Tesselate a shape into basic primitives.

#### Function tesselate_shapes()

    void ygl::tesselate_shapes(scene* scn, bool subdivide, bool
    facevarying_to_sharedvertex, bool quads_to_triangles, bool
    bezier_to_lines);

Tesselate scene shapes.

#### Function update_transforms()

    void ygl::update_transforms(scene* scn, float time=0);

Update node transforms.

#### Function compute_animation_range()

    vec2f ygl::compute_animation_range(const scene* scn);

Compute animation range.

#### Function make_view_camera()

    camera* ygl::make_view_camera(const scene* scn, int camera_id);

Make a view camera either copying a given one or building a default
one.

#### Function compute_bounds()

    bbox3f ygl::compute_bounds(const shape* shp);

Computes a shape bounding box using a quick computation that ignores
radius.

#### Function compute_bounds()

    bbox3f ygl::compute_bounds(const scene* scn);

Compute a scene bounding box.

#### Function flatten_instances()

    void ygl::flatten_instances(scene* scn);

Flatten scene instances into separate shapes.

#### Function print_info()

    void ygl::print_info(const scene* scn);

Print scene information.

#### Function make_bvh()

    bvh_tree* ygl::make_bvh(const shape* shp, float def_radius=0.001f, bool
    equalsize=true);

Build a shape BVH.

#### Function make_bvh()

    bvh_tree* ygl::make_bvh(const scene* scn, float def_radius=0.001f, bool
    equalsize=true);

Build a scene BVH.

#### Function refit_bvh()

    void ygl::refit_bvh(bvh_tree* bvh, const shape* shp, float
    def_radius=0.001f);

Refits a scene BVH.

#### Function refit_bvh()

    void ygl::refit_bvh(bvh_tree* bvh, const scene* scn, bool do_shapes, float
    def_radius=0.001f);

Refits a scene BVH.

#### Function add_elements()

    void ygl::add_elements(scene* scn, const add_elements_options& opts={});

Add elements.

#### Function merge_into()

    void ygl::merge_into(scene* merge_into, scene* merge_from);

Merge scene into one another. Note that the objects are moved from
merge_from to merged_into, so merge_from will be empty after this
function.

#### Function load_scene()

    scene* ygl::load_scene(const std::string& filename, const load_options&
    opts={});

Loads a scene. For now OBJ or glTF are supported. Throws an exception
if an error occurs.

#### Function save_scene()

    void ygl::save_scene(const std::string& filename, const scene* scn, const
    save_options& opts);

Saves a scene. For now OBJ and glTF are supported.

## Scene type support

#### Function enum_names< material_type >()

    const std::vector<std::pair<std::string, material_type>>&
    ygl::enum_names<material_type>();

Names of enum values.

#### Function enum_names< keyframe_type >()

    const std::vector<std::pair<std::string, keyframe_type>>&
    ygl::enum_names<keyframe_type>();

Names of enum values.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(camera& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(texture& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(texture_info& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(material& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(shape& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(shape_group& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(instance& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(environment& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(node& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(animation& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(animation_group& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(scene& val, Visitor& &visitor);

Visit struct elements.

## Example scenes

#### Enum test_texture_type

    enum struct test_texture_type {

        none,
        grid,
        checker,
        colored,
        rcolored,
        bump,
        uv,
        gamma,
        noise,
        ridge,
        fbm,
        turbulence,
        gammaf,
        sky,
    }

Test texture type.

Members:
    - none: None (empty texture).
    - grid: Grid image.
    - checker: Checker image.
    - colored: Colored checker image.
    - rcolored: Detailed colored checker image.
    - bump: Bump and dimple imahe.
    - uv: Uv debug image.
    - gamma: Gamma ramp.
    - noise: Perlin noise image.
    - ridge: Perlin ridge image.
    - fbm: Perlin fbm image.
    - turbulence: Perlin turbulence image.
    - gammaf: Gamma ramp (HDR).
    - sky: Sky (HDR).

#### Enum test_material_type

    enum struct test_material_type {

        none,
        emission,
        matte,
        plastic,
        metal,
        transparent,
    }

Test material type.

Members:
    - none: None (empty material).
    - emission: Emission.
    - matte: Matte (diffuse).
    - plastic: Plastic.
    - metal: Matetal.
    - transparent: Transparent (diffuse with opacity).

#### Enum test_shape_type

    enum struct test_shape_type {

        floor,
        quad,
        cube,
        sphere,
        spherecube,
        spherizedcube,
        geosphere,
        flipcapsphere,
        suzanne,
        cubep,
        fvcube,
        fvsphere,
        matball,
        point,
        pointscube,
        hairball,
        beziercircle,
    }

Test shape type.

Members:
    - floor: Floor (shared vertex, 20x20 size).
    - quad: Quad (shared vertex).
    - cube: Cube (shared vertex, not watertight).
    - sphere: Sphere (shared vertex, not watertight).
    - spherecube: Sphere with cube uvs (shared vertex, not watertight).
    - spherizedcube: Spherized cube (shared vertex, not watertight).
    - geosphere: Geodesic sphere (shared vertex, watertight, no texcoord).
    - flipcapsphere: Sphere with flipped cap (shared vertex, not watertight).
    - suzanne: Suzanne (shared vertex, no texcoord).
    - cubep: Position-only cube (shared vertex).
    - fvcube: Face-varying cube (shared vertex).
    - fvsphere: Face-varying sphere (shared vertex).
    - matball: Matball (shared vertex, not watertight).
    - point: Single point.
    - pointscube: Random points in a cube.
    - hairball: Random lines on a sphere.
    - beziercircle: Bezier circle.

#### Function make_cornell_box_scene()

    scene* ygl::make_cornell_box_scene();

Makes the Cornell Box scene.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, camera* cam, const
    test_camera_params& tcam);

Updates a test camera.

#### Function test_camera_presets()

    std::unordered_map<std::string, test_camera_params>&
    ygl::test_camera_presets();

Test camera presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, texture* txt, const
    test_texture_params& ttxt);

Updates a test texture.

#### Function test_texture_presets()

    std::unordered_map<std::string, test_texture_params>&
    ygl::test_texture_presets();

Test texture presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, material* mat, const
    test_material_params& tmat);

Updates a test material.

#### Function test_material_presets()

    std::unordered_map<std::string, test_material_params>&
    ygl::test_material_presets();

Test material presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, shape* shp, const
    test_shape_params& tshp);

Updates a test shape, adding it to the scene if missing.

#### Function test_shape_presets()

    std::unordered_map<std::string, test_shape_params>&
    ygl::test_shape_presets();

Test shape presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, instance* ist, const
    test_instance_params& tist);

Updates a test instance.

#### Function test_instance_presets()

    std::unordered_map<std::string, test_instance_params>&
    ygl::test_instance_presets();

Test instance presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, environment* env, const
    test_environment_params& tenv);

Updates a test instance.

#### Function test_environment_presets()

    std::unordered_map<std::string, test_environment_params>&
    ygl::test_environment_presets();

Test environment presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, node* nde, const
    test_node_params& tndr);

Updates a test node.

#### Function test_node_presets()

    std::unordered_map<std::string, test_node_params>&
    ygl::test_node_presets();

Test nodes presets.

#### Function update_test_elem()

    void ygl::update_test_elem(const scene* scn, animation_group* anm, const
    test_animation_params& tndr);

Updates a test node.

#### Function update_test_scene()

    void ygl::update_test_scene(scene* scn, const test_scene_params& tscn,
    const std::unordered_set<void*>& refresh={});

Updates a test scene, adding missing objects. Objects are only added
and never removed.

#### Function make_test_scene()

    scene* ygl::make_test_scene(const test_scene_params& tscn);

Makes a test scene. Convenience wrapper around update_test_scene().

#### Function test_scene_presets()

    std::unordered_map<std::string, test_scene_params>&
    ygl::test_scene_presets();

Test scene presets.

#### Function remove_duplicates()

    void ygl::remove_duplicates(test_scene_params& tscn);

Remove duplicates based on name.

#### Function load_test_scene()

    test_scene_params ygl::load_test_scene(const std::string& filename);

Load test scene.

#### Function save_test_scene()

    void ygl::save_test_scene(const std::string& filename, const
    test_scene_params& scn);

Save test scene.

## Example scenes type support

#### Function enum_names< test_texture_type >()

    const std::vector<std::pair<std::string, test_texture_type>>&
    ygl::enum_names<test_texture_type>();

Names of enum values.

#### Function enum_names< test_material_type >()

    const std::vector<std::pair<std::string, test_material_type>>&
    ygl::enum_names<test_material_type>();

Names of enum values.

#### Function enum_names< test_shape_type >()

    const std::vector<std::pair<std::string, test_shape_type>>&
    ygl::enum_names<test_shape_type>();

Names of enum values.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_camera_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_texture_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_material_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_shape_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_instance_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_environment_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_node_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_animation_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(test_scene_params& val, Visitor& &visitor);

Visit struct elements.

## Path-tracing support

#### Function specular_exponent_to_roughness()

    float ygl::specular_exponent_to_roughness(float n);

Phong exponent to roughness.

#### Function specular_fresnel_from_ks()

    void ygl::specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f&
    esk);

Specular to fresnel eta.

#### Function fresnel_dielectric()

    vec3f ygl::fresnel_dielectric(float cosw, const vec3f& eta_);

Compute the fresnel term for dielectrics. Implementation from
https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-
equations/.

#### Function fresnel_metal()

    vec3f ygl::fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);

Compute the fresnel term for metals. Implementation from
https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-
equations/.

#### Function fresnel_schlick()

    vec3f ygl::fresnel_schlick(const vec3f& ks, float cosw);

Schlick approximation of Fresnel term.

#### Function fresnel_schlick()

    vec3f ygl::fresnel_schlick(const vec3f& ks, float cosw, float rs);

Schlick approximation of Fresnel term weighted by roughness. This is a
hack, but works better than not doing it.

#### Function eval_ggx()

    float ygl::eval_ggx(float rs, float ndh, float ndi, float ndo);

Evaluates the GGX distribution and geometric term.

#### Function sample_ggx()

    vec3f ygl::sample_ggx(float rs, const vec2f& rn);

Sample the GGX distribution.

#### Function sample_ggx_pdf()

    float ygl::sample_ggx_pdf(float rs, float ndh);

Evaluates the GGX pdf.

#### Function filter_triangle()

    float ygl::filter_triangle(float x);

Triangle filter. Ppublic domain from stb_image_resize.

#### Function filter_cubic()

    float ygl::filter_cubic(float x);

Cubic filter. Ppublic domain from stb_image_resize.

#### Function filter_catmullrom()

    float ygl::filter_catmullrom(float x);

Catmull-rom filter. Ppublic domain from stb_image_resize.

#### Function filter_mitchell()

    float ygl::filter_mitchell(float x);

Mitchell filter. Ppublic domain from stb_image_resize.

## Path tracing

#### Enum trace_shader_type

    enum struct trace_shader_type {

        pathtrace,
        eyelight,
        direct,
        pathtrace_nomis,
        debug_normal,
        debug_albedo,
        debug_texcoord,
    }

Type of rendering algorithm.

Members:
    - pathtrace: Pathtrace.
    - eyelight: Eye light for previews.
    - direct: Direct illumination.
    - pathtrace_nomis: Pathtrace without MIS.
    - debug_normal: Debug normal.
    - debug_albedo: Debug albedo.
    - debug_texcoord: Debug texcoord.

#### Enum trace_rng_type

    enum struct trace_rng_type {

        uniform,
        stratified,
    }

Random number generator type.

Members:
    - uniform: Uniform random numbers.
    - stratified: Stratified random numbers.

#### Enum trace_filter_type

    enum struct trace_filter_type {

        box,
        triangle,
        cubic,
        catmull_rom,
        mitchell,
    }

Filter type.

Members:
    - box: Box filter.
    - triangle: Hat filter.
    - cubic: Cubic spline.
    - catmull_rom: Catmull-Rom spline.
    - mitchell: Mitchell-Netrevalli.

#### Function make_trace_pixels()

    image<trace_pixel> ygl::make_trace_pixels(const image4f& img, const
    trace_params& params);

Initialize trace pixels.

#### Function make_trace_lights()

    trace_lights ygl::make_trace_lights(const scene* scn);

Initialize trace lights.

#### Function trace_samples()

    void ygl::trace_samples(const scene* scn, const camera* cam, const
    bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, int nsamples, const trace_params& params);

Trace the next nsamples samples.

#### Function trace_samples_filtered()

    void ygl::trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, int nsamples, const trace_params& params);

Trace the next nsamples samples with image filtering.

#### Function trace_image()

    image4f ygl::trace_image(const scene* scn, const camera* cam, const
    bvh_tree* bvh, const trace_params& params);

Trace the whole image.

#### Function trace_async_start()

    void ygl::trace_async_start(const scene* scn, const camera* cam, const
    bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, std::vector<std::thread>& threads, bool&
    stop_flag, const trace_params& params);

Starts an anyncrhounous renderer.

#### Function trace_async_stop()

    void ygl::trace_async_stop(std::vector<std::thread>& threads, bool&
    stop_flag);

Stop the asynchronous renderer.

## Path tracing type support

#### Function enum_names< trace_shader_type >()

    const std::vector<std::pair<std::string, trace_shader_type>>&
    ygl::enum_names<trace_shader_type>();

Names of enum values.

#### Function enum_names< trace_rng_type >()

    const std::vector<std::pair<std::string, trace_rng_type>>&
    ygl::enum_names<trace_rng_type>();

Names of enum values.

#### Function enum_names< trace_filter_type >()

    const std::vector<std::pair<std::string, trace_filter_type>>&
    ygl::enum_names<trace_filter_type>();

Names of enum values.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(trace_params& val, Visitor& &visitor);

Visit struct elements.

## Wavefront OBJ

#### Enum obj_element_type

    enum struct obj_element_type {

        point,
        line,
        face,
        bezier,
    }

Obj element type.

Members:
    - point: List of points.
    - line: Polyline.
    - face: Polygon face.
    - bezier: Bezier segments.

#### Function operator==()

    bool ygl::operator==(const obj_vertex& a, const obj_vertex& b);



#### Function operator==()

    bool ygl::operator==(const obj_texture_info& a, const obj_texture_info&
    b);



#### Function load_obj()

    obj_scene* ygl::load_obj(const std::string& filename, bool
    load_textures=false, bool skip_missing=false, bool flip_texcoord=true,
    bool flip_tr=true);

Load an OBJ from file filename. Load textures if load_textures is
true, and report errors only if skip_missing is false. Texture
coordinates and material Tr are flipped if flip_texcoord and flip_tp
are respectively true.

#### Function save_obj()

    void ygl::save_obj(const std::string& filename, const obj_scene* model,
    bool save_textures=false, bool skip_missing=false, bool
    flip_texcoord=true, bool flip_tr=true);

Save an OBJ to file filename. Save textures if save_textures is true,
and report errors only if skip_missing is false. Texture coordinates
and material Tr are flipped if flip_texcoord and flip_tp are
respectively true.

## Khronos glTF

#### Enum glTFAccessorSparseIndicesComponentType

    enum struct glTFAccessorSparseIndicesComponentType {

        NotSet,
        UnsignedByte,
        UnsignedShort,
        UnsignedInt,
    }

Values for glTFAccessorSparseIndices::componentType.

Members:
    - NotSet: Not set.
    - UnsignedByte: 
    - UnsignedShort: 
    - UnsignedInt: 

#### Enum glTFAccessorComponentType

    enum struct glTFAccessorComponentType {

        NotSet,
        Byte,
        UnsignedByte,
        Short,
        UnsignedShort,
        UnsignedInt,
        Float,
    }

Values for glTFAccessor::componentType.

Members:
    - NotSet: Not set.
    - Byte: 
    - UnsignedByte: 
    - Short: 
    - UnsignedShort: 
    - UnsignedInt: 
    - Float: 

#### Enum glTFAccessorType

    enum struct glTFAccessorType {

        NotSet,
        Scalar,
        Vec2,
        Vec3,
        Vec4,
        Mat2,
        Mat3,
        Mat4,
    }

Values for glTFAccessor::type.

Members:
    - NotSet: Not set.
    - Scalar: 
    - Vec2: 
    - Vec3: 
    - Vec4: 
    - Mat2: 
    - Mat3: 
    - Mat4: 

#### Enum glTFAnimationChannelTargetPath

    enum struct glTFAnimationChannelTargetPath {

        NotSet,
        Translation,
        Rotation,
        Scale,
        Weights,
    }

Values for glTFAnimationChannelTarget::path.

Members:
    - NotSet: Not set.
    - Translation: 
    - Rotation: 
    - Scale: 
    - Weights: 

#### Enum glTFAnimationSamplerInterpolation

    enum struct glTFAnimationSamplerInterpolation {

        NotSet,
        Linear,
        Step,
        CatmullRomSpline,
        CubicSpline,
    }

Values for glTFAnimationSampler::interpolation.

Members:
    - NotSet: Not set.
    - Linear: 
    - Step: 
    - CatmullRomSpline: 
    - CubicSpline: 

#### Enum glTFBufferViewTarget

    enum struct glTFBufferViewTarget {

        NotSet,
        ArrayBuffer,
        ElementArrayBuffer,
    }

Values for glTFBufferView::target.

Members:
    - NotSet: Not set.
    - ArrayBuffer: 
    - ElementArrayBuffer: 

#### Enum glTFCameraType

    enum struct glTFCameraType {

        NotSet,
        Perspective,
        Orthographic,
    }

Values for glTFCamera::type.

Members:
    - NotSet: Not set.
    - Perspective: 
    - Orthographic: 

#### Enum glTFImageMimeType

    enum struct glTFImageMimeType {

        NotSet,
        ImageJpeg,
        ImagePng,
    }

Values for glTFImage::mimeType.

Members:
    - NotSet: Not set.
    - ImageJpeg: 
    - ImagePng: 

#### Enum glTFMaterialAlphaMode

    enum struct glTFMaterialAlphaMode {

        NotSet,
        Opaque,
        Mask,
        Blend,
    }

Values for glTFMaterial::alphaMode.

Members:
    - NotSet: Not set.
    - Opaque: 
    - Mask: 
    - Blend: 

#### Enum glTFMeshPrimitiveMode

    enum struct glTFMeshPrimitiveMode {

        NotSet,
        Points,
        Lines,
        LineLoop,
        LineStrip,
        Triangles,
        TriangleStrip,
        TriangleFan,
    }

Values for glTFMeshPrimitive::mode.

Members:
    - NotSet: Not set.
    - Points: 
    - Lines: 
    - LineLoop: 
    - LineStrip: 
    - Triangles: 
    - TriangleStrip: 
    - TriangleFan: 

#### Enum glTFSamplerMagFilter

    enum struct glTFSamplerMagFilter {

        NotSet,
        Nearest,
        Linear,
    }

Values for glTFSampler::magFilter.

Members:
    - NotSet: Not set.
    - Nearest: 
    - Linear: 

#### Enum glTFSamplerMinFilter

    enum struct glTFSamplerMinFilter {

        NotSet,
        Nearest,
        Linear,
        NearestMipmapNearest,
        LinearMipmapNearest,
        NearestMipmapLinear,
        LinearMipmapLinear,
    }

Values for glTFSampler::minFilter.

Members:
    - NotSet: Not set.
    - Nearest: 
    - Linear: 
    - NearestMipmapNearest: 
    - LinearMipmapNearest: 
    - NearestMipmapLinear: 
    - LinearMipmapLinear: 

#### Enum glTFSamplerWrapS

    enum struct glTFSamplerWrapS {

        NotSet,
        ClampToEdge,
        MirroredRepeat,
        Repeat,
    }

glTFSampler::wrapS

Members:
    - NotSet: Not set.
    - ClampToEdge: 
    - MirroredRepeat: 
    - Repeat: 

#### Enum glTFSamplerWrapT

    enum struct glTFSamplerWrapT {

        NotSet,
        ClampToEdge,
        MirroredRepeat,
        Repeat,
    }

glTFSampler::wrapT

Members:
    - NotSet: Not set.
    - ClampToEdge: 
    - MirroredRepeat: 
    - Repeat: 

#### Typedef buffer_data

    using ygl::buffer_data = std::vector<unsigned char>;

Generic buffer data.

#### Function load_gltf()

    glTF* ygl::load_gltf(const std::string& filename, bool load_bin=true, bool
    load_img=false, bool skip_missing=false);

Load a gltf file filename from disk. Load binaries and images only if
load_bin and load_img are true, reporting errors only if skip_missing
is false.

#### Function load_binary_gltf()

    glTF* ygl::load_binary_gltf(const std::string& filename, bool
    load_bin=true, bool load_img=false, bool skip_missing=false);

Load a binary gltf file filename from disk. Load binaries and images
only if load_bin and load_img are true, reporting errors only if
skip_missing is false.

#### Function save_gltf()

    void ygl::save_gltf(const std::string& filename, const glTF* gltf, bool
    save_bin=true, bool save_img=false);

Save a gltf file filename to disk. Save binaries and images only if
save_bin and save_img are true.

#### Function save_binary_gltf()

    void ygl::save_binary_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin=true, bool save_img=false);

Save a gltf file filename to disk. Save binaries and images only if
save_bin and save_img are true.

#### Function node_transform()

    mat4f ygl::node_transform(const glTFNode* node);

Computes the local node transform and its inverse.

## Svg

#### Function load_svg()

    svg_scene* ygl::load_svg(const std::string& filename);

Load an SVG.

#### Function save_svg()

    void ygl::save_svg(const std::string& filename, const svg_scene* svg);

Save an SVG.

## String, path and file functions

#### Function startswith()

    bool ygl::startswith(const std::string& str, const std::string& substr);

Checks if a string starts with a prefix.

#### Function endswith()

    bool ygl::endswith(const std::string& str, const std::string& substr);

Checks if a string ends with a prefix.

#### Function contains()

    bool ygl::contains(const std::string& str, const std::string& substr);

Check is a string contains a substring.

#### Function splitlines()

    std::vector<std::string> ygl::splitlines(const std::string& str, bool
    keep_newline=false);

Splits a string into lines at the ' ' character. The line terminator
is kept if keep_newline. This function does not work on Window if
keep_newline is true.

#### Function partition()

    std::vector<std::string> ygl::partition(const std::string& str, const
    std::string& split);

Partition the string.

#### Function split()

    std::vector<std::string> ygl::split(const std::string& str);

Splits the string.

#### Function split()

    std::vector<std::string> ygl::split(const std::string& str, const
    std::string& substr);

Splits the string.

#### Function split()

    std::vector<std::string> ygl::split(const std::string& str, char substr);

Splits the string.

#### Function rstrip()

    std::string ygl::rstrip(const std::string& str);

Strip the string.

#### Function lstrip()

    std::string ygl::lstrip(const std::string& str);

Strip the string.

#### Function strip()

    std::string ygl::strip(const std::string& str);

Strip the string.

#### Function join()

    std::string ygl::join(const std::vector<std::string>& strs, const
    std::string& sep);

Joins a list of string with a string as separator.

#### Function lower()

    std::string ygl::lower(const std::string& str);

Converts an ASCII string to lowercase.

#### Function upper()

    std::string ygl::upper(const std::string& str);

Converts an ASCII string to uppercase.

#### Function isspace()

    bool ygl::isspace(const std::string& str);

Check if a string is space.

#### Function replace()

    std::string ygl::replace(const std::string& str, const std::string& s1,
    const std::string& s2);

Replace s1 with s2 in str.

#### Function path_dirname()

    std::string ygl::path_dirname(const std::string& filename);

Get directory name (including '/').

#### Function path_extension()

    std::string ygl::path_extension(const std::string& filename);

Get extension (including '.').

#### Function path_basename()

    std::string ygl::path_basename(const std::string& filename);

Get file basename.

#### Function path_filename()

    std::string ygl::path_filename(const std::string& filename);

Get filename without directory (equiv to get_basename() +
get_extension()).

#### Function replace_path_extension()

    std::string ygl::replace_path_extension(const std::string& filename, const
    std::string& ext);

Replace extension.

#### Function prepend_path_extension()

    std::string ygl::prepend_path_extension(const std::string& filename, const
    std::string& prep);

Prepend a string to the extension.

#### Function split_path()

    void ygl::split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext);

Splits a path calling the above functions.

#### Function path_convert_eparator()

    std::string ygl::path_convert_eparator(const std::string& path_);

Convert from Windows to Unix/OsX path separator.

#### Function format()

    std::string ygl::format(const std::string& fmt, const
    std::vector<std::string>& args);

Really-minimal Python like string format. The implementation is not
fast nor memory efficient. But it is good enough for some needs.

#### Function _format_one()

    void ygl::_format_one(std::vector<std::string>& vals);



#### Function _format_one()

    template<typename Arg, typename...>
    void ygl::_format_one(std::vector<std::string>& vals, const Arg& arg,
    const Args& ... args);



#### Function format()

    template<typename...>
    std::string ygl::format(const std::string& fmt, const Args& ... args);

Really-minimal Python like string format. Internally uses streams for
generality and supports for now only the '{}' operator. The
implementation is not fast nor memory efficient. But it is good enough
for some needs.

#### Function print()

    template<typename...>
    void ygl::print(const std::string& fmt, const Args& ... args);

Wrapper for the above function that prints to stdout.

#### Function println()

    template<typename...>
    void ygl::println(const std::string& fmt, const Args& ... args);

Wrapper for the above function that prints to stdout with endline.

## File loading and saving

#### Function load_binary()

    std::vector<unsigned char> ygl::load_binary(const std::string& filename);

Loads the contents of a binary file in an in-memory array.

#### Function load_text()

    std::string ygl::load_text(const std::string& filename);

Loads the contents of a text file into a string.

#### Function save_binary()

    void ygl::save_binary(const std::string& filename, const
    std::vector<unsigned char>& data);

Saves binary data to a file.

#### Function save_text()

    void ygl::save_text(const std::string& filename, const std::string& str);

Saves a string to a text file.

## Immediate-mode command line parser

#### Function should_exit()

    bool ygl::should_exit(cmdline_parser& parser);

Check unused arguments.

#### Function get_usage()

    std::string ygl::get_usage(const cmdline_parser& parser);

Returns the usage string.

#### Function parse_flag()

    bool ygl::parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def=false, bool
    req=false);

Pase a flag from the command line.

#### Function parse_opt()

    template<typename T>
    T ygl::parse_opt(cmdline_parser& parser, const std::string& name, const
    std::string& flag, const std::string& help, const T& def={}, bool
    req=false, const std::vector<T>& choices={});

Pase an option from the command line.

#### Function parse_opt()

    template<typename T>
    T ygl::parse_opt(cmdline_parser& parser, const std::string& name, const
    std::string& flag, const std::string& help, const
    std::vector<std::pair<std::string, T>>& key_values, const T& def, bool
    req=false, const std::vector<T>& choices={});

Parse an enum option from the command line.

#### Function parse_arg()

    template<typename T>
    T ygl::parse_arg(cmdline_parser& parser, const std::string& name, const
    std::string& help, const T& def={}, bool req=true, const std::vector<T>&
    choices={});

Parse positional argument from the command line.

#### Function parse_args()

    template<typename T>
    std::vector<T> ygl::parse_args(cmdline_parser& parser, const std::string&
    name, const std::string& help, const std::vector<T>& def={}, bool
    req=true, const std::vector<T>& choices={});

Parse all remaining positional argument from the command line.

#### Function parse_params()

    template<typename T>
    T ygl::parse_params(cmdline_parser& parser, const std::string& name, const
    T& def={}, bool req=false);

Parse options generated with a visit over the parameters.

#### Function make_parser()

    cmdline_parser ygl::make_parser(int argc, char* *argv, const std::string&
    prog, const std::string& help);

Initialize a command line parser.

## Simple logging

#### Function make_logger()

    logger* ygl::make_logger(const std::string& filename="", bool
    console=true, bool verbose=true, bool file_append=true);

Make a logger with an optional console stream, an optional file stram
and the specified verbosity level.

#### Function get_default_logger()

    logger* ygl::get_default_logger();

Get the default logger.

#### Function _log_msg()

    void ygl::_log_msg(logger* lgr, const std::string& msg, const char* type);



#### Function log_info()

    template<typename...>
    void ygl::log_info(logger* lgr, const std::string& msg, const Args& ...
    args);

Log an info message.

#### Function log_warning()

    template<typename...>
    void ygl::log_warning(logger* lgr, const std::string& msg, const Args& ...
    args);

Log an info message.

#### Function log_error()

    template<typename...>
    void ygl::log_error(logger* lgr, const std::string& msg, const Args& ...
    args);

Log an error message.

#### Function log_fatal()

    template<typename...>
    void ygl::log_fatal(logger* lgr, const std::string& msg, const Args& ...
    args);

Log a fatal message and exit.

#### Function log_info()

    template<typename...>
    void ygl::log_info(const std::string& msg, const Args& ... args);

Logs a message to the default loggers.

#### Function log_error()

    template<typename...>
    void ygl::log_error(const std::string& msg, const Args& ... args);

Logs a message to the default loggers.

#### Function log_fatal()

    template<typename...>
    void ygl::log_fatal(const std::string& msg, const Args& ... args);

Logs a message to the default loggers.

## Simple timer

## OpenGL utilities

#### Enum gl_elem_type

    enum struct gl_elem_type {

        point,
        line,
        triangle,
    }

OpenGL shape element types.

Members:
    - point: Points.
    - line: Lines.
    - triangle: Triangles.

#### Enum gl_light_type

    enum struct gl_light_type {

        point,
        directional,
    }

OpenGL light types.

Members:
    - point: Point lights.
    - directional: Directional lights.

#### Function gl_check_error()

    bool ygl::gl_check_error(bool print=true);

Checks for GL error and then prints.

#### Function gl_clear_buffers()

    void ygl::gl_clear_buffers(const vec4f& background={0, 0, 0, 0});

Clear window.

#### Function gl_enable_depth_test()

    void ygl::gl_enable_depth_test(bool enabled);

Enable/disable depth test.

#### Function gl_enable_culling()

    void ygl::gl_enable_culling(bool enabled, bool front=false, bool
    back=true);

Enable/disable culling.

#### Function gl_enable_wireframe()

    void ygl::gl_enable_wireframe(bool enabled);

Enable/disable wireframe.

#### Function gl_enable_blending()

    void ygl::gl_enable_blending(bool enabled);

Enable/disable blending.

#### Function gl_set_blend_over()

    void ygl::gl_set_blend_over();

Set blending to over operator.

#### Function gl_line_width()

    void ygl::gl_line_width(float w);

Line width.

#### Function gl_set_viewport()

    void ygl::gl_set_viewport(const vec4i& v);

Set viewport.

#### Function gl_set_viewport()

    void ygl::gl_set_viewport(const vec2i& v);

Set viewport.

#### Function gl_read_imagef()

    void ygl::gl_read_imagef(float* pixels, int w, int h, int nc);

Reads an image from the the framebuffer.

## OpenGL textures

#### Enum gl_texture_wrap

    enum struct gl_texture_wrap {

        not_set,
        repeat,
        clamp,
        mirror,
    }

Wrap values for OpenGL texture.

Members:
    - not_set: Not set.
    - repeat: Repeat.
    - clamp: Clamp to edge.
    - mirror: Mirror.

#### Enum gl_texture_filter

    enum struct gl_texture_filter {

        not_set,
        linear,
        nearest,
        linear_mipmap_linear,
        nearest_mipmap_nearest,
        linear_mipmap_nearest,
        nearest_mipmap_linear,
    }

Filter values for OpenGL texture.

Members:
    - not_set: Not set.
    - linear: Linear.
    - nearest: Nearest.
    - linear_mipmap_linear: Mip-mapping.
    - nearest_mipmap_nearest: Mip-mapping.
    - linear_mipmap_nearest: Mip-mapping.
    - nearest_mipmap_linear: Mip-mapping.

#### Function _update_texture()

    void ygl::_update_texture(gl_texture& txt, int w, int h, int nc, const
    void* pixels, bool floats, bool linear, bool mipmap, bool as_float, bool
    as_srgb);



#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, int w, int h, int nc, const
    float* pixels, bool linear, bool mipmap, bool as_float);

Updates a texture with pixels values of size w, h with nc number of
components (1-4). Internally use float if as_float and filtering if
filter.

#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, int w, int h, int nc, const
    unsigned char* pixels, bool linear, bool mipmap, bool as_srgb);

Updates a texture with pixels values of size w, h with nc number of
components (1-4). Internally use float if as_float and filtering if
filter.

#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, const image4f& img, bool linear,
    bool mipmap, bool as_float);

Updates a texture with pixels values from an image. Internally use
float if as_float and filtering if filter.

#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, const image4b& img, bool linear,
    bool mipmap, bool as_srgb);

Updates a texture with pixels values from an image. Internally use
float if as_float and filtering if filter.

#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, const image4f& img);

Updates a texture with pixels values from an image.

#### Function update_texture()

    void ygl::update_texture(gl_texture& txt, const image4b& img);

Updates a texture with pixels values from an image.

#### Function make_texture()

    gl_texture ygl::make_texture(const image4f& img, bool linear, bool mipmap,
    bool as_float);

Creates a texture from an image. Convenience wrapper to
update_texture().

#### Function make_texture()

    gl_texture ygl::make_texture(const image4b& img, bool linear, bool mipmap,
    bool as_srgb);

Creates a texture from an image. Convenience wrapper to
update_texture().

#### Function bind_texture()

    void ygl::bind_texture(const gl_texture& txt, uint unit);

Binds a texture to a texture unit.

#### Function unbind_texture()

    void ygl::unbind_texture(const gl_texture& txt, uint unit);

Unbinds a texture.

#### Function get_texture_id()

    uint ygl::get_texture_id(const gl_texture& txt);

Get texture id.

#### Function is_texture_valid()

    bool ygl::is_texture_valid(const gl_texture& txt);

Check if defined.

#### Function clear_texture()

    void ygl::clear_texture(gl_texture& txt);

Destroys the texture tid.

## OpenGL vertex array buffers

#### Function _update_vertex_buffer()

    void ygl::_update_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic);



#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, int num, int ncomp,
    const float* values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, int num, int ncomp,
    const int* values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<float>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec2f>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec3f>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec4f>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<int>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec2i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec3i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_vertex_buffer()

    void ygl::update_vertex_buffer(gl_vertex_buffer& buf, const
    std::vector<vec4i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function make_vertex_buffer()

    template<typename T>
    gl_vertex_buffer ygl::make_vertex_buffer(const std::vector<T>& values,
    bool dynamic=false);

Make a buffer with new data.

#### Function bind_vertex_buffer()

    void ygl::bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);

Bind the buffer at a particular attribute location.

#### Function unbind_vertex_buffer()

    void ygl::unbind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);

Unbind the buffer.

#### Function unbind_vertex_buffer()

    void ygl::unbind_vertex_buffer(uint vattr);

Unbind the buffer.

#### Function get_vertex_buffer_id()

    uint ygl::get_vertex_buffer_id(const gl_vertex_buffer& buf);

Get buffer id.

#### Function is_vertex_buffer_valid()

    bool ygl::is_vertex_buffer_valid(const gl_vertex_buffer& buf);

Check if defined.

#### Function clear_vertex_buffer()

    void ygl::clear_vertex_buffer(gl_vertex_buffer& buf);

Destroys the buffer.

## OpenGL element array buffers

#### Function _update_element_buffer()

    void ygl::_update_element_buffer(gl_element_buffer& buf, int n, int nc,
    const int* values, bool dynamic);



#### Function update_element_buffer()

    void ygl::update_element_buffer(gl_element_buffer& buf, int num, int
    ncomp, const int* values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_element_buffer()

    void ygl::update_element_buffer(gl_element_buffer& buf, const
    std::vector<int>& values, bool dynamic=false);

Updates the bufferwith new data.

#### Function update_element_buffer()

    void ygl::update_element_buffer(gl_element_buffer& buf, const
    std::vector<vec2i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_element_buffer()

    void ygl::update_element_buffer(gl_element_buffer& buf, const
    std::vector<vec3i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function update_element_buffer()

    void ygl::update_element_buffer(gl_element_buffer& buf, const
    std::vector<vec4i>& values, bool dynamic=false);

Updates the buffer with new data.

#### Function make_element_buffer()

    template<typename T>
    gl_element_buffer ygl::make_element_buffer(const std::vector<T>& values,
    bool dynamic=false);

Make a buffer with new data.

#### Function draw_elems()

    void ygl::draw_elems(const gl_element_buffer& buf);

Draws elements.

#### Function get_element_buffer_id()

    uint ygl::get_element_buffer_id(const gl_element_buffer& buf);

Get id.

#### Function is_element_buffer_valid()

    bool ygl::is_element_buffer_valid(const gl_element_buffer& buf);

Check if defined.

#### Function clear_element_buffer()

    void ygl::clear_element_buffer(gl_element_buffer& buf);

Destroys the buffer.

## OpenGL programs

#### Function make_program()

    gl_program ygl::make_program(const std::string& vertex, const std::string&
    fragment);

Creates an OpenGL program from vertex and fragment code.

#### Function clear_program()

    void ygl::clear_program(gl_program& prog);

Destroys the program.

#### Function get_program_uniform_location()

    int ygl::get_program_uniform_location(const gl_program& prog, const
    std::string& name);

Get uniform location.

#### Function get_program_attrib_location()

    int ygl::get_program_attrib_location(const gl_program& prog, const
    std::string& name);

Get attribute location.

#### Function get_program_uniforms_names()

    std::vector<std::pair<std::string, int>>
    ygl::get_program_uniforms_names(const gl_program& prog);

Get the names of all uniforms.

#### Function get_program_attributes_names()

    std::vector<std::pair<std::string, int>>
    ygl::get_program_attributes_names(const gl_program& prog);

Get the names of all attributes.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int pos, const int* val,
    int ncomp, int count);

Set uniform integer values val for program prog and variable pos. The
values have nc number of components (1-4) and count elements.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int pos, const float* val,
    int ncomp, int count);

Set uniform float values val for program prog and variable pos. The
values have nc number of components (1-4) and count elements.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int var, bool val);

Set uniform value.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int var, int val);

Set uniform value.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int var, float val);

Set uniform value.

#### Function set_program_uniform()

    template<typename T, int>
    bool ygl::set_program_uniform(gl_program& prog, int var, const vec<T, N>&
    val);

Set uniform value.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, int var, const mat<T, 4>&
    val);

Set uniform value.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, int var, const frame<T,
    3>& val);

Set uniform value.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int var, const int* val,
    int num);

Set uniform array.

#### Function set_program_uniform()

    bool ygl::set_program_uniform(gl_program& prog, int var, const float* val,
    int num);

Set uniform array.

#### Function set_program_uniform()

    template<typename T, int>
    bool ygl::set_program_uniform(gl_program& prog, int var, const vec<T, N>*
    val, int num);

Set uniform array.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, int var, const mat<T, 4>*
    val, int num);

Set uniform array.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, int var, const frame<T,
    3>* val, int num);

Set uniform array.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, const std::string& var,
    const T& val);

Set uniform value for names variable.

#### Function set_program_uniform()

    template<typename T>
    bool ygl::set_program_uniform(gl_program& prog, const std::string& var,
    const T* val, int num);

Set uniform array for names variable.

#### Function set_program_uniform_texture()

    bool ygl::set_program_uniform_texture(gl_program& prog, int pos, const
    gl_texture_info& tinfo, uint tunit);

Set uniform texture.

#### Function set_program_uniform_texture()

    bool ygl::set_program_uniform_texture(gl_program& prog, int var, int
    varon, const gl_texture_info& tinfo, uint tunit);

Set uniform texture with an additionasl texture enable flags.

#### Function set_program_uniform_texture()

    bool ygl::set_program_uniform_texture(gl_program& prog, const std::string&
    var, const gl_texture_info& tinfo, uint tunit);

Set uniform texture.

#### Function set_program_uniform_texture()

    bool ygl::set_program_uniform_texture(gl_program& prog, const std::string&
    var, const std::string& varon, const gl_texture_info& tinfo, uint tunit);

Set uniform texture with an additionasl texture enable flags.

#### Function set_program_vertattr()

    bool ygl::set_program_vertattr(gl_program& prog, int pos, const float*
    value, int nc);

Sets a constant value of nc components for the vertex attribute at pos
location.

#### Function set_program_vertattr()

    bool ygl::set_program_vertattr(gl_program& prog, int pos, const int*
    value, int nc);

Sets a constant value of nc components for the vertex attribute at pos
location.

#### Function set_program_vertattr()

    bool ygl::set_program_vertattr(gl_program& prog, const std::string& var,
    const gl_vertex_buffer& buf);

Binds a buffer to a vertex attribute.

#### Function set_program_vertattr()

    bool ygl::set_program_vertattr(gl_program& prog, int pos, const
    gl_vertex_buffer& buf, int nc, const float* def);

Binds a buffer to a vertex attribute, or a constant value if the
buffer is empty.

#### Function set_program_vertattr()

    template<typename T, int>
    bool ygl::set_program_vertattr(gl_program& prog, int var, const
    gl_vertex_buffer& buf, const vec<T, N>& def);

Binds a buffer or constant to a vertex attribute.

#### Function set_program_vertattr()

    template<typename T, int>
    bool ygl::set_program_vertattr(gl_program& prog, const std::string& var,
    const gl_vertex_buffer& buf, const vec<T, N>& def);

Binds a buffer or constant to a vertex attribute.

#### Function is_program_valid()

    bool ygl::is_program_valid(const gl_program& prog);

Check whether the program is valid.

#### Function bind_program()

    void ygl::bind_program(const gl_program& prog);

Binds a program.

#### Function unbind_program()

    void ygl::unbind_program(const gl_program& prog);

Unbind a program.

## OpenGL scene shader support

#### Function clear_shape()

    void ygl::clear_shape(gl_shape& shp);

Clear shape.

#### Function make_gl_lights()

    gl_lights ygl::make_gl_lights(const scene* scn);

Initialize gl lights.

#### Function clear_textures()

    void ygl::clear_textures(std::unordered_map<texture*, gl_texture>&
    textures);

Clear scene textures on the GPU.

#### Function clear_shapes()

    void ygl::clear_shapes(std::unordered_map<shape*, gl_shape>& shapes);

Clear scene shapes on the GPU.

#### Function update_textures()

    void ygl::update_textures(const scene* scn, std::unordered_map<texture*,
    gl_texture>& textures, const std::unordered_set<texture*>& refresh={},
    bool clear=false);

Update scene textures on the GPU.

#### Function update_shapes()

    void ygl::update_shapes(const scene* scn, std::unordered_map<shape*,
    gl_shape>& shapes, const std::unordered_set<shape*>& refresh={}, const
    std::unordered_set<shape_group*>& refreshg={}, bool clear=false);

Update scene shapes on the GPU.

## OpenGL image shader

#### Function make_stdimage_program()

    gl_stdimage_program ygl::make_stdimage_program();

Initialize a stdimage program.

#### Function draw_image()

    void ygl::draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma, bool filmic);

Draws an image texture the stdimage program.

#### Function draw_image()

    void ygl::draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom);

Draws an image texture the stdimage program.

#### Function draw_image()

    void ygl::draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const gl_stdimage_params& params, bool
    clear_background=true);

Draws an image texture the stdimage program.

## OpenGL surface shader

#### Function make_stdsurface_program()

    gl_stdsurface_program ygl::make_stdsurface_program();

Initialize a stdsurface shader.

#### Function is_program_valid()

    bool ygl::is_program_valid(const gl_stdsurface_program& prog);

Check if the program is valid.

#### Function begin_stdsurface_frame()

    void ygl::begin_stdsurface_frame(gl_stdsurface_program& prog, bool
    shade_eyelight, float tonemap_exposure, float tonemap_gamma, bool
    tonemap_filmic, const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj);

Starts a frame by setting exposure/gamma values, camera transforms and
projection. Sets also whether to use full shading or a quick eye light
preview.

#### Function end_stdsurface_frame()

    void ygl::end_stdsurface_frame(gl_stdsurface_program& prog);

Ends a frame.

#### Function set_stdsurface_lights()

    void ygl::set_stdsurface_lights(gl_stdsurface_program& prog, const vec3f&
    amb, const gl_lights& lights);

Set shading lights and ambient.

#### Function begin_stdsurface_shape()

    void ygl::begin_stdsurface_shape(gl_stdsurface_program& prog, const mat4f&
    xform, float normal_offset=0);

Begins drawing a shape with transform xform.

#### Function end_stdsurface_shape()

    void ygl::end_stdsurface_shape(gl_stdsurface_program& prog);

End shade drawing.

#### Function set_stdsurface_normaloffset()

    void ygl::set_stdsurface_normaloffset(gl_stdsurface_program& prog, float
    normal_offset);

Sets normal offset.

#### Function set_stdsurface_highlight()

    void ygl::set_stdsurface_highlight(gl_stdsurface_program& prog, const
    vec4f& highlight);

Set the object as highlighted.

#### Function set_stdsurface_material()

    void ygl::set_stdsurface_material(gl_stdsurface_program& prog,
    material_type type, gl_elem_type etype, const vec3f& ke, const vec3f& kd,
    const vec3f& ks, float rs, float op, const gl_texture_info& ke_txt, const
    gl_texture_info& kd_txt, const gl_texture_info& ks_txt, const
    gl_texture_info& rs_txt, const gl_texture_info& norm_txt, const
    gl_texture_info& occ_txt, bool use_phong, bool double_sided, bool
    alpha_cutout);

Set material values with emission ke, diffuse kd, specular ks and
specular roughness rs, opacity op. Indicates textures ids with the
correspoinding XXX_txt variables. Sets also normal and occlusion maps.
Works for points/lines/triangles indicated by etype, (diffuse for
points, Kajiya-Kay for lines, GGX/Phong for triangles). Material type
matches the scene material type.

#### Function set_stdsurface_constmaterial()

    void ygl::set_stdsurface_constmaterial(gl_stdsurface_program& prog, const
    vec3f& ke, float op);

Set constant material with emission ke and opacity op.

#### Function set_stdsurface_vert()

    void ygl::set_stdsurface_vert(gl_stdsurface_program& prog, const
    gl_vertex_buffer& pos, const gl_vertex_buffer& norm, const
    gl_vertex_buffer& texcoord, const gl_vertex_buffer& color, const
    gl_vertex_buffer& tangsp);

Set vertex data with buffers for position pos, normals norm, texture
coordinates texcoord, per-vertex color color and tangent space tangsp.

#### Function set_stdsurface_vert_skinning()

    void ygl::set_stdsurface_vert_skinning(gl_stdsurface_program& prog, const
    gl_vertex_buffer& weights, const gl_vertex_buffer& joints, int nxforms,
    const mat4f* xforms);

Set vertex data with buffers for skinning.

#### Function set_stdsurface_vert_gltf_skinning()

    void ygl::set_stdsurface_vert_gltf_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints, int
    nxforms, const mat4f* xforms);

Set vertex data with buffers for skinning.

#### Function set_stdsurface_vert_skinning_off()

    void ygl::set_stdsurface_vert_skinning_off(gl_stdsurface_program& prog);

Disables vertex skinning.

#### Function draw_stdsurface_scene()

    void ygl::draw_stdsurface_scene(const scene* scn, const camera* cam,
    gl_stdsurface_program& prog, std::unordered_map<shape*, gl_shape>& shapes,
    std::unordered_map<texture*, gl_texture>& textures, const gl_lights&
    lights, const vec2i& viewport_size, const gl_stdsurface_params& params);

Draw scene with stdsurface program.

## OpenGL standard shaders type support

#### Function visit()

    template<typename Visitor>
    void ygl::visit(gl_stdimage_params& val, Visitor& &visitor);

Visit struct elements.

#### Function visit()

    template<typename Visitor>
    void ygl::visit(gl_stdsurface_params& val, Visitor& &visitor);

Visit struct elements.

## OpenGL window

#### Typedef gl_text_callback

    void(* ygl::gl_text_callback) (gl_window*, unsigned int);

Text callback.

#### Typedef gl_mouse_callback

    void(* ygl::gl_mouse_callback) (gl_window*, int button, bool
    press, int mods);

Mouse callback.

#### Typedef gl_refresh_callback

    void(* ygl::gl_refresh_callback) (gl_window* );

Window refresh callback.

#### Function make_window()

    gl_window* ygl::make_window(int width, int height, const std::string&
    title, void* user_pointer=nullptr);

Initialize a window.

#### Function set_window_callbacks()

    void ygl::set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb);

Set window callbacks.

#### Function clear_window()

    void ygl::clear_window(gl_window* win);

Clear window.

#### Function get_user_pointer()

    void* ygl::get_user_pointer(gl_window* win);

Gets the user poiner.

#### Function set_window_title()

    void ygl::set_window_title(gl_window* win, const std::string& title);

Set window title.

#### Function wait_events()

    void ygl::wait_events(gl_window* win);

Wait events.

#### Function poll_events()

    void ygl::poll_events(gl_window* win);

Poll events.

#### Function swap_buffers()

    void ygl::swap_buffers(gl_window* win);

Swap buffers.

#### Function should_close()

    bool ygl::should_close(gl_window* win);

Should close.

#### Function get_window_size()

    vec2i ygl::get_window_size(gl_window* win);

Window size.

#### Function get_framebuffer_size()

    vec2i ygl::get_framebuffer_size(gl_window* win);

Framebuffer size.

#### Function get_mouse_button()

    int ygl::get_mouse_button(gl_window* win);

Mouse button.

#### Function get_mouse_pos()

    vec2i ygl::get_mouse_pos(gl_window* win);

Mouse position.

#### Function get_mouse_posf()

    vec2f ygl::get_mouse_posf(gl_window* win);

Mouse position.

#### Function get_key()

    bool ygl::get_key(gl_window* win, int key);

Check if a key is pressed (not all keys are supported)

#### Function get_screenshot()

    std::vector<vec4b> ygl::get_screenshot(gl_window* win, vec2i& wh, bool
    flipy=true, bool back=false);

Read pixels.

#### Function save_screenshot()

    void ygl::save_screenshot(gl_window* win, const std::string& imfilename);

Save a screenshot to disk.

#### Function handle_camera_navigation()

    bool ygl::handle_camera_navigation(gl_window* win, camera* cam, bool
    navigation_fps);

Handle camera navigation.

## OpenGL widgets

#### Function init_widgets()

    void ygl::init_widgets(gl_window* win, bool light_style=false, bool
    extra_font=true);

Initialize widgets.

#### Function begin_widgets()

    bool ygl::begin_widgets(gl_window* win, const std::string& title);

Begin draw widgets.

#### Function end_widgets()

    void ygl::end_widgets(gl_window* win);

End draw widgets.

#### Function get_widget_active()

    bool ygl::get_widget_active(gl_window* win);

Whether widgets are active.

#### Function draw_separator_widget()

    void ygl::draw_separator_widget(gl_window* win);

Horizontal separator.

#### Function draw_indent_widget_begin()

    void ygl::draw_indent_widget_begin(gl_window* win);

Indent widget.

#### Function draw_indent_widget_end()

    void ygl::draw_indent_widget_end(gl_window* win);

Indent widget.

#### Function draw_continue_widget()

    void ygl::draw_continue_widget(gl_window* win);

Continue line with next widget.

#### Function draw_label_widget()

    void ygl::draw_label_widget(gl_window* win, const std::string& lbl, const
    std::string& msg);

Label widget.

#### Function draw_label_widget()

    template<typename...>
    void ygl::draw_label_widget(gl_window* win, const std::string& lbl, const
    std::string& fmt, const Args& ... args);

Label widget.

#### Function draw_label_widget()

    template<typename T>
    void ygl::draw_label_widget(gl_window* win, const std::string& lbl, const
    T& val);

Label widget.

#### Function draw_checkbox_widget()

    bool ygl::draw_checkbox_widget(gl_window* win, const std::string& lbl,
    bool& val);

Checkbox widget.

#### Function draw_text_widget()

    bool ygl::draw_text_widget(gl_window* win, const std::string& lbl,
    std::string& str);

Text widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl, int&
    val, int min=0, int max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec2i& val, int min=0, int max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec3i& val, int min=0, int max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec4i& val, int min=0, int max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    float& val, float min=0, float max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec2f& val, float min=0, float max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec3f& val, float min=0, float max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    vec4f& val, float min=0, float max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    mat<float, 4>& val, float min=0, float max=1);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    frame<float, 3>& val, float min=-10, float max=10);

Slider widget.

#### Function draw_slider_widget()

    bool ygl::draw_slider_widget(gl_window* win, const std::string& lbl,
    quat4f& val, float min=-1, float max=1);

Slider widget.

#### Function draw_color_widget()

    bool ygl::draw_color_widget(gl_window* win, const std::string& lbl, vec4f&
    val);

Color widget.

#### Function draw_color_widget()

    bool ygl::draw_color_widget(gl_window* win, const std::string& lbl, vec4b&
    val);

Color widget.

#### Function draw_color_widget()

    bool ygl::draw_color_widget(gl_window* win, const std::string& lbl, vec3f&
    val);

Color widget.

#### Function draw_combo_widget_begin()

    bool ygl::draw_combo_widget_begin(gl_window* win, const std::string& lbl,
    const std::string& label);

Combo widget.

#### Function draw_combo_widget_item()

    bool ygl::draw_combo_widget_item(gl_window* win, const std::string& label,
    int idx, bool selected);

Combo widget.

#### Function draw_combo_widget_end()

    void ygl::draw_combo_widget_end(gl_window* win);

Combo widget.

#### Function draw_combo_widget_item()

    template<typename T>
    bool ygl::draw_combo_widget_item(gl_window* win, const std::string& label,
    int idx, T& val, const T& item);

Combo widget.

#### Function draw_combo_widget()

    template<typename T, typename T1>
    bool ygl::draw_combo_widget(gl_window* win, const std::string& lbl, T&
    val, const std::vector<T1>& vals, const std::function<T(const T1& )>&
    value_func, const std::function<std::string(const T1& )>& label_func);

Combo widget.

#### Function draw_combo_widget()

    bool ygl::draw_combo_widget(gl_window* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels);

Combo widget.

#### Function draw_combo_widget()

    template<typename T>
    bool ygl::draw_combo_widget(gl_window* win, const std::string& lbl, T&
    val, const std::vector<std::pair<std::string, T>>& labels);

Combo widget.

#### Function draw_combo_widget()

    template<typename T>
    bool ygl::draw_combo_widget(gl_window* win, const std::string& lbl, T*&
    val, const std::vector<T*>& vals, bool extra=true, T* extra_val=nullptr);

Combo widget.

#### Function draw_button_widget()

    bool ygl::draw_button_widget(gl_window* win, const std::string& lbl);

Button widget.

#### Function draw_header_widget()

    bool ygl::draw_header_widget(gl_window* win, const std::string& lbl);

Collapsible header widget.

#### Function draw_tree_widget_begin()

    bool ygl::draw_tree_widget_begin(gl_window* win, const std::string& lbl);

Start tree widget.

#### Function draw_tree_widget_end()

    void ygl::draw_tree_widget_end(gl_window* win);

End tree widget.

#### Function draw_tree_widget_begin()

    bool ygl::draw_tree_widget_begin(gl_window* win, const std::string& lbl,
    void*& selection, void* content);

Start selectable tree node widget.

#### Function draw_tree_widget_begin()

    bool ygl::draw_tree_widget_begin(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);

Start selectable tree node widget.

#### Function draw_tree_widget_end()

    void ygl::draw_tree_widget_end(gl_window* win, void* content);

End selectable tree node widget.

#### Function draw_tree_widget_leaf()

    void ygl::draw_tree_widget_leaf(gl_window* win, const std::string& lbl,
    void*& selection, void* content);

Selectable tree leaf nodewidget.

#### Function draw_tree_widget_leaf()

    void ygl::draw_tree_widget_leaf(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);

Selectable tree leaf node widget.

#### Function draw_tree_widget_color_begin()

    void ygl::draw_tree_widget_color_begin(gl_window* win, const vec4f&
    color);

Text color widget.

#### Function draw_tree_widget_color_end()

    void ygl::draw_tree_widget_color_end(gl_window* win);

Text color widget.

#### Function draw_image_widget()

    void ygl::draw_image_widget(gl_window* win, int tid, const vec2i& size,
    const vec2i& imsize);

Image widget.

#### Function draw_image_widget()

    void ygl::draw_image_widget(gl_window* win, gl_texture& txt, const vec2i&
    size);

Image widget.

#### Function draw_scroll_widget_begin()

    void ygl::draw_scroll_widget_begin(gl_window* win, const std::string& lbl,
    int height, bool border);

Scroll region widget.

#### Function draw_scroll_widget_end()

    void ygl::draw_scroll_widget_end(gl_window* win);

Scroll region widget.

#### Function draw_scroll_widget_here()

    void ygl::draw_scroll_widget_here(gl_window* win);

Scroll region widget.

#### Function draw_groupid_widget_begin()

    void ygl::draw_groupid_widget_begin(gl_window* win, int gid);

Group ids widget.

#### Function draw_groupid_widget_begin()

    void ygl::draw_groupid_widget_begin(gl_window* win, void* gid);

Group ids widget.

#### Function draw_groupid_widget_begin()

    void ygl::draw_groupid_widget_begin(gl_window* win, const char* gid);

Group ids widget.

#### Function draw_groupid_widget_end()

    void ygl::draw_groupid_widget_end(gl_window* win);

Group ids widget.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, bool&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Min, max and color are
ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl,
    std::string& val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Min, max and color are
ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, int&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    template<typename T, int>
    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl,
    vec<int, N>& val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, float&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, vec2f&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, vec3f&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, vec4f&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl,
    frame3f& val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max for frame
origin, or a deafult range when their are the same. Color is ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, mat4f&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl,
    quat4f& val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Uses min and max, or a deafult
range when their are the same. Color is ignored.

#### Function draw_value_widget()

    template<typename T, typename std::enable_if<std::is_enum<T>::value,
    int>::type>
    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl, T&
    val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Min, max and color are
ignored.

#### Function draw_value_widget()

    bool ygl::draw_value_widget(gl_window* win, const std::string& lbl,
    uint32_t& val, float min=0, float max=0, bool color=false);

Generic widget used for templated code. Internally convert to int
loosing precision. See the int version.

#### Function draw_imageinspect_widgets()

    void ygl::draw_imageinspect_widgets(gl_window* win, const std::string&
    lbl, const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos, const
    gl_stdimage_params& params);

Image inspection widgets.

#### Function draw_params_widgets()

    template<typename T>
    bool ygl::draw_params_widgets(gl_window* win, const std::string& lbl, T&
    params);

Draws a widget that sets params in non-recursive trivial structures.
Internally uses visit to implement the view.

#### Function draw_camera_selection_widget()

    bool ygl::draw_camera_selection_widget(gl_window* win, const std::string&
    lbl, camera*& cam, scene* scn, camera* view);

Draws a widget that can selected the camera.

#### Function draw_camera_widgets()

    bool ygl::draw_camera_widgets(gl_window* win, const std::string& lbl,
    camera* cam);

Draws widgets for a camera. Used for quickly making demos.

#### Function draw_scene_widgets()

    bool ygl::draw_scene_widgets(gl_window* win, const std::string& lbl,
    scene* scn, void*& selection, const std::unordered_map<texture*,
    gl_texture>& gl_txt, test_scene_params* test_scn=nullptr);

Draws widgets for a whole scene. Used for quickly making demos.

