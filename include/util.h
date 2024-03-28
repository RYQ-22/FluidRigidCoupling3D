#ifndef UTIL_H
#define UTIL_H

#include "field.h"

#include <algorithm>
#include <vector>
#include <cmath>
#include <iostream>

namespace appconstants {
    const double PI = 3.1415926535897932384626433832795;
    const double G = 9.8;
}


#ifdef WIN32
#undef min
#undef max
#endif

using std::min;
using std::max;
using std::swap;


inline int sign(const float& x)
{
    return std::signbit(x) ? -1 : 1;
}

template<class T>
inline T sqr(const T& x)
{
    return x * x;
}

template<class T>
inline T cube(const T& x)
{
    return x * x * x;
}

template<class T>
inline T min(T a1, T a2, T a3)
{
    return min(a1, min(a2, a3));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4)
{
    return min(min(a1, a2), min(a3, a4));
}

template<class T>
inline T second(T a1, T a2, T a3) {
    return a1 + a2 + a3 - max(a1, max(a2, a3)) - min(a1, min(a2, a3));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5)
{
    return min(min(a1, a2), min(a3, a4), a5);
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return min(min(a1, a2), min(a3, a4), min(a5, a6));
}

template<class T>
inline T min(T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
{
    return min(min(a1, a2), min(a3, a4), min(a5, a6, a7, a8));
}

template<class T>
inline T max(T a1, T a2, T a3)
{
    return max(a1, max(a2, a3));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4)
{
    return max(max(a1, a2), max(a3, a4));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5)
{
    return max(max(a1, a2), max(a3, a4), a5);
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5, T a6)
{
    return max(max(a1, a2), max(a3, a4), max(a5, a6));
}

template<class T>
inline T max(T a1, T a2, T a3, T a4, T a5, T a6, T a7, T a8)
{
    return max(max(a1, a2), max(a3, a4), max(a5, a6, a7, a8));
}

template<class T>
inline void minmax(T a1, T a2, T& amin, T& amax)
{
    if (a1 < a2) {
        amin = a1;
        amax = a2;
    }
    else {
        amin = a2;
        amax = a1;
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T& amin, T& amax)
{
    if (a1 < a2) {
        if (a1 < a3) {
            amin = a1;
            if (a2 < a3) amax = a3;
            else amax = a2;
        }
        else {
            amin = a3;
            if (a1 < a2) amax = a2;
            else amax = a1;
        }
    }
    else {
        if (a2 < a3) {
            amin = a2;
            if (a1 < a3) amax = a3;
            else amax = a1;
        }
        else {
            amin = a3;
            amax = a1;
        }
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T& amin, T& amax)
{
    if (a1 < a2) {
        if (a3 < a4) {
            amin = min(a1, a3);
            amax = max(a2, a4);
        }
        else {
            amin = min(a1, a4);
            amax = max(a2, a3);
        }
    }
    else {
        if (a3 < a4) {
            amin = min(a2, a3);
            amax = max(a1, a4);
        }
        else {
            amin = min(a2, a4);
            amax = max(a1, a3);
        }
    }
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T& amin, T& amax)
{
    //@@@ the logic could be shortcircuited a lot!
    amin = min(a1, a2, a3, a4, a5);
    amax = max(a1, a2, a3, a4, a5);
}

template<class T>
inline void minmax(T a1, T a2, T a3, T a4, T a5, T a6, T& amin, T& amax)
{
    //@@@ the logic could be shortcircuited a lot!
    amin = min(a1, a2, a3, a4, a5, a6);
    amax = max(a1, a2, a3, a4, a5, a6);
}

template<class T>
inline void update_minmax(T a1, T& amin, T& amax)
{
    if (a1 < amin) amin = a1;
    else if (a1 > amax) amax = a1;
}

template<class T>
inline void sort(T& a, T& b, T& c)
{
    T temp;
    if (a < b) {
        if (a < c) {
            if (c < b) { // a<c<b
                temp = c; c = b; b = temp;
            } // else: a<b<c
        }
        else { // c<a<b
            temp = c; c = b; b = a; a = temp;
        }
    }
    else {
        if (b < c) {
            if (a < c) { //b<a<c
                temp = b; b = a; a = temp;
            }
            else { // b<c<a
                temp = b; b = c; c = a; a = temp;
            }
        }
        else { // c<b<a
            temp = c; c = a; a = temp;
        }
    }
}

template<class T>
inline T clamp(T a, T lower, T upper)
{
    if (a < lower) return lower;
    else if (a > upper) return upper;
    else return a;
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r)
{
    if (r < 0) return 0;
    else if (r > 1) return 1;
    return r * r * r * (10 + r * (-15 + r * 6));
}

// only makes sense with T=float or double
template<class T>
inline T smooth_step(T r, T r_lower, T r_upper, T value_lower, T value_upper)
{
    return value_lower + smooth_step((r - r_lower) / (r_upper - r_lower)) * (value_upper - value_lower);
}

// only makes sense with T=float or double
template<class T>
inline T ramp(T r)
{
    return smooth_step((r + 1) / 2) * 2 - 1;
}

#ifdef _MSC_VER

inline double remainder(double x, double y)
{
    return x - std::floor(x / y + 0.5) * y;
}
#endif

inline unsigned int round_up_to_power_of_two(unsigned int n)
{
    int exponent = 0;
    --n;
    while (n) {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

inline unsigned int round_down_to_power_of_two(unsigned int n)
{
    int exponent = 0;
    while (n > 1) {
        ++exponent;
        n >>= 1;
    }
    return 1 << exponent;
}

// Transforms even the sequence 0,1,2,3,... into reasonably good random numbers 
// Challenge: improve on this in speed and "randomness"!
// This seems to pass several statistical tests, and is a bijective map (of 32-bit unsigned ints)
inline unsigned int randhash(unsigned int seed)
{
    unsigned int i = (seed ^ 0xA3C59AC3u) * 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    i ^= (i >> 16);
    i *= 2654435769u;
    return i;
}

// the inverse of randhash
inline unsigned int unhash(unsigned int h)
{
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= (h >> 16);
    h *= 340573321u;
    h ^= 0xA3C59AC3u;
    return h;
}

// returns repeatable stateless pseudo-random number in [0,1]
inline double randhashd(unsigned int seed)
{
    return randhash(seed) / (double)UINT_MAX;
}
inline float randhashf(unsigned int seed)
{
    return randhash(seed) / (float)UINT_MAX;
}

// returns repeatable stateless pseudo-random number in [a,b]
inline double randhashd(unsigned int seed, double a, double b)
{
    return (b - a) * randhash(seed) / (double)UINT_MAX + a;
}
inline float randhashf(unsigned int seed, float a, float b)
{
    return ((b - a) * randhash(seed) / (float)UINT_MAX + a);
}

inline int intlog2(int x)
{
    int exp = -1;
    while (x) {
        x >>= 1;
        ++exp;
    }
    return exp;
}

template<class T>
inline void get_barycentric(T x, int& i, T& f, int i_low, int i_high)
{
    T s = std::floor(x);
    i = (int)s;
    if (i < i_low) {
        i = i_low;
        f = 0;
    }
    else if (i > i_high - 2) {
        i = i_high - 2;
        f = 1;
    }
    else
        f = (T)(x - s);
}

template<class S, class T>
inline S lerp(const S& value0, const S& value1, T f)
{
    return (1 - f) * value0 + f * value1;
}

template<class S, class T>
inline S bilerp(const S& v00, const S& v10,
    const S& v01, const S& v11,
    T fx, T fy)
{
    return lerp(lerp(v00, v10, fx),
        lerp(v01, v11, fx),
        fy);
}

template<class S, class T>
inline S trilerp(const S& v000, const S& v100,
    const S& v010, const S& v110,
    const S& v001, const S& v101,
    const S& v011, const S& v111,
    T fx, T fy, T fz)
{
    return lerp(bilerp(v000, v100, v010, v110, fx, fy),
        bilerp(v001, v101, v011, v111, fx, fy),
        fz);
}

template<class S, class T>
inline S quadlerp(const S& v0000, const S& v1000,
    const S& v0100, const S& v1100,
    const S& v0010, const S& v1010,
    const S& v0110, const S& v1110,
    const S& v0001, const S& v1001,
    const S& v0101, const S& v1101,
    const S& v0011, const S& v1011,
    const S& v0111, const S& v1111,
    T fx, T fy, T fz, T ft)
{
    return lerp(trilerp(v0000, v1000, v0100, v1100, v0010, v1010, v0110, v1110, fx, fy, fz),
        trilerp(v0001, v1001, v0101, v1101, v0011, v1011, v0111, v1111, fx, fy, fz),
        ft);
}

// f should be between 0 and 1, with f=0.5 corresponding to balanced weighting between w0 and w2
template<class T>
inline void quadratic_bspline_weights(T f, T& w0, T& w1, T& w2)
{
    w0 = T(0.5) * sqr(f - 1);
    w1 = T(0.75) - sqr(f - T(0.5));;
    w2 = T(0.5) * sqr(f);
}

// f should be between 0 and 1
template<class T>
inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
{
    T f2(f * f), f3(f2 * f);
    wneg1 = -T(1. / 3) * f + T(1. / 2) * f2 - T(1. / 6) * f3;
    w0 = 1 - f2 + T(1. / 2) * (f3 - f);
    w1 = f + T(1. / 2) * (f2 - f3);
    w2 = T(1. / 6) * (f3 - f);
}

template<class S, class T>
inline S cubic_interp(const S& value_neg1, const S& value0, const S& value1, const S& value2, T f)
{
    T wneg1, w0, w1, w2;
    cubic_interp_weights(f, wneg1, w0, w1, w2);
    return wneg1 * value_neg1 + w0 * value0 + w1 * value1 + w2 * value2;
}

template<class T>
void zero(std::vector<T>& v)
{
    for (int i = (int)v.size() - 1; i >= 0; --i) v[i] = 0;
}

template<class T>
T abs_max(const std::vector<T>& v)
{
    T m = 0;
    for (int i = (int)v.size() - 1; i >= 0; --i) {
        if (std::fabs(v[i]) > m)
            m = std::fabs(v[i]);
    }
    return m;
}

template<class T>
bool contains(const std::vector<T>& a, T e)
{
    for (unsigned int i = 0; i < a.size(); ++i)
        if (a[i] == e) return true;
    return false;
}

template<class T>
void add_unique(std::vector<T>& a, T e)
{
    for (unsigned int i = 0; i < a.size(); ++i)
        if (a[i] == e) return;
    a.push_back(e);
}

template<class T>
void insert(std::vector<T>& a, unsigned int index, T e)
{
    a.push_back(a.back());
    for (unsigned int i = (unsigned int)a.size() - 1; i > index; --i)
        a[i] = a[i - 1];
    a[index] = e;
}

template<class T>
void erase(std::vector<T>& a, unsigned int index)
{
    for (unsigned int i = index; i < a.size() - 1; ++i)
        a[i] = a[i + 1];
    a.pop_back();
}

template<class T>
void erase_swap(std::vector<T>& a, unsigned int index)
{
    for (unsigned int i = index; i < a.size() - 1; ++i)
        swap(a[i], a[i + 1]);
    a.pop_back();
}

template<class T>
void erase_unordered(std::vector<T>& a, unsigned int index)
{
    a[index] = a.back();
    a.pop_back();
}

template<class T>
void erase_unordered_swap(std::vector<T>& a, unsigned int index)
{
    swap(a[index], a.back());
    a.pop_back();
}

template<class T>
void find_and_erase_unordered(std::vector<T>& a, const T& doomed_element)
{
    for (unsigned int i = 0; i < a.size(); ++i)
        if (a[i] == doomed_element) {
            erase_unordered(a, i);
            return;
        }
}

template<class T>
void replace_once(std::vector<T>& a, const T& old_element, const T& new_element)
{
    for (unsigned int i = 0; i < a.size(); ++i)
        if (a[i] == old_element) {
            a[i] = new_element;
            return;
        }
}

template<class T>
void write_matlab(std::ostream& output, const std::vector<T>& a, const char* variable_name, bool column_vector = true, int significant_digits = 18)
{
    output << variable_name << "=[";
    std::streamsize old_precision = output.precision();
    output.precision(significant_digits);
    for (unsigned int i = 0; i < a.size(); ++i) {
        output << a[i] << " ";
    }
    output << "]";
    if (column_vector)
        output << "'";
    output << ";" << std::endl;
    output.precision(old_precision);
}

template<class T>
T sample_min(const Eigen::Vector3d& pos, const Field3<T>& v) {
    assert(pos.size() == 3);
    int i0 = pos(0) > 0 ? (int)pos(0) : (int)pos(0) - 1;
    int j0 = pos(1) > 0 ? (int)pos(1) : (int)pos(1) - 1;
    int k0 = pos(2) > 0 ? (int)pos(2) : (int)pos(2) - 1;
    return min(
        v(i0, j0, k0),
        v(i0 + 1, j0, k0),
        v(i0, j0 + 1, k0),
        v(i0, j0, k0 + 1),
        v(i0 + 1, j0 + 1, k0),
        v(i0, j0 + 1, k0 + 1),
        v(i0 + 1, j0, k0 + 1),
        v(i0 + 1, j0 + 1, k0 + 1));
}

template<class T>
T sample_max(const Eigen::Vector3d& pos, const Field3<T>& v) {
    int i0 = pos(0) > 0 ? (int)pos(0) : (int)pos(0) - 1;
    int j0 = pos(1) > 0 ? (int)pos(1) : (int)pos(1) - 1;
    int k0 = pos(2) > 0 ? (int)pos(2) : (int)pos(2) - 1;
    return max(
        v(i0, j0, k0),
        v(i0 + 1, j0, k0),
        v(i0, j0 + 1, k0),
        v(i0, j0, k0 + 1),
        v(i0 + 1, j0 + 1, k0),
        v(i0, j0 + 1, k0 + 1),
        v(i0 + 1, j0, k0 + 1),
        v(i0 + 1, j0 + 1, k0 + 1));
}

template<class T>
T interpolate_value(double i, double j, double k, const Field3<T>& v) {
    int i0 = i > 0 ? (int)i : (int)i - 1;
    int j0 = j > 0 ? (int)j : (int)j - 1;
    int k0 = k > 0 ? (int)k : (int)k - 1;
    float a1 = i0 + 1.0f - i, a2 = i - i0;
    float b1 = j0 + 1.0f - j, b2 = j - j0;
    float c1 = k0 + 1.0f - k, c2 = k - k0;
    return v(i0, j0, k0) * a1 * b1 * c1
        + v(i0 + 1, j0 + 1, k0 + 1) * a2 * b2 * c2
        + v(i0 + 1, j0, k0) * a2 * b1 * c1
        + v(i0, j0 + 1, k0) * a1 * b2 * c1
        + v(i0, j0, k0 + 1) * a1 * b1 * c2
        + v(i0 + 1, j0 + 1, k0) * a2 * b2 * c1
        + v(i0 + 1, j0, k0 + 1) * a2 * b1 * c2
        + v(i0, j0 + 1, k0 + 1) * a1 * b2 * c2;
} 

template<class T>
T interpolate_value(const Eigen::Vector3d& pos, const Field3<T>& v) {
    assert(pos.size() == 3);
    return interpolate_value(pos(0), pos(1), pos(2), v);
}

template<class T>
void interpolate_gradient(Eigen::Matrix<T, 3, 1>& grad, double i, double j, double k, const Field3<T>& v) {
    int i0 = i > 0 ? (int)i : (int)i - 1;
    int j0 = j > 0 ? (int)j : (int)j - 1;
    int k0 = k > 0 ? (int)k : (int)k - 1;
    float a1 = i0 + 1.0f - i, a2 = i - i0;
    float b1 = j0 + 1.0f - j, b2 = j - j0;
    float c1 = k0 + 1.0f - k, c2 = k - k0;
    grad.resize(3);  
    // x
    T grad_x0 = v(i0 + 1, j0, k0) - v(i0, j0, k0);
    T grad_x1 = v(i0 + 1, j0 + 1, k0) - v(i0, j0 + 1, k0);
    T grad_x2 = v(i0 + 1, j0, k0 + 1) - v(i0, j0, k0 + 1);
    T grad_x3 = v(i0 + 1, j0 + 1, k0 + 1) - v(i0, j0 + 1, k0 + 1);
    grad(0) = grad_x0 * b1 * c1 + grad_x3 * b2 * c2 + grad_x2 * b1 * c2 + grad_x1 * b2 * c1;
    // y
    T grad_y0 = v(i0, j0 + 1, k0) - v(i0, j0, k0);
    T grad_y1 = v(i0 + 1, j0 + 1, k0) - v(i0 + 1, j0, k0);
    T grad_y2 = v(i0, j0 + 1, k0 + 1) - v(i0, j0, k0 + 1);
    T grad_y3 = v(i0 + 1, j0 + 1, k0 + 1) - v(i0 + 1, j0, k0 + 1);
    grad(1) = grad_y0 * a1 * c1 + grad_y3 * a2 * c2 + grad_y2 * a1 * c2 + grad_y1 * a2 * c1;
    // z
    T grad_z0 = v(i0, j0, k0 + 1) - v(i0, j0, k0);
    T grad_z1 = v(i0 + 1, j0, k0 + 1) - v(i0 + 1, j0, k0);
    T grad_z2 = v(i0, j0 + 1, k0 + 1) - v(i0, j0 + 1, k0);
    T grad_z3 = v(i0 + 1, j0 + 1, k0 + 1) - v(i0 + 1, j0 + 1, k0);
    grad(2) = grad_z0 * a1 * b1 + grad_z3 * a2 * b2 + grad_z2 * a1 * b2 + grad_z1 * a2 * b1;
    return;
}

template<class T>
void interpolate_gradient(Eigen::Matrix<T, 3, 1>& grad, const Eigen::Vector3d& pos, const Field3<T>& v) {
    assert(pos.size() == 3);
    interpolate_gradient(grad, pos(0), pos(1), pos(2), v);
}

inline double computeTriangleArea(double phi0, double phi1, double phi2) {
    sort(phi0, phi1, phi2);
    if(phi0 > 0){
        return 0;
    }
    if(phi2 < 0){
        return 1;
    }
    if(phi1 < 0){
        return phi0/(phi0-phi2) * phi1/(phi1-phi2);
    }
    else{
        return 1 - phi1/(phi1-phi0) * phi2/(phi2-phi0);
    }
}

template<class T>
inline std::vector<T> union_phi(const std::vector<T>& phi1, const std::vector<T>& phi2) {
    if (phi1.size() != phi2.size()) {
        abort();
    }
    std::vector<T> ans = phi1;
    for (int i = 0; i < ans.size(); i++) {
        if (abs(phi1[i]) > abs(phi2[i])) {
            ans[i] = phi2[i];
        }
    }
    return ans;
}


#endif