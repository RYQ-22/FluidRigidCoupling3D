#ifndef FIELD_H
#define FIELD_H

#include "log.h"

namespace backend {

template<class T>
class Field3 : public Eigen::Matrix<T, Eigen::Dynamic, 1> {
private:
    int N1_, N2_, N3_;    
    
public:
    Field3() : N1_(0), N2_(0), N3_(0) {}

    Field3(const int& N1, const int& N2, const int& N3) : N1_(N1), N2_(N2), N3_(N3) {  
        assert(N1>0 && N2>0 && N3>0);
        this->resize(N1 * N2 * N3);
    }

    T& operator()(int i, int j, int k) {
        assert(i>=0 && i<N1_ && j>=0 && j<N2_ && k>=0 && k<N3_);
        return this->coeffRef(i * N2_ * N3_ + j * N3_ + k);
    }

    const T& operator()(int i, int j, int k) const {
        assert(i>=0 && i<N1_ && j>=0 && j<N2_ && k>=0 && k<N3_);
        return this->coeffRef(i * N2_ * N3_ + j * N3_ + k);
    }

    Field3<T> operator+(const Field3<T>& other) const {
        assert(this->N1_ == other.N1_ && this->N2_ == other.N2_ && this->N3_ == other.N3_);
        Field3<T> result = *this;
        result += other;        
        result.N1_ = this->N1_;
        result.N2_ = this->N2_;
        result.N3_ = this->N3_;
        return result;
    }

    Field3<T> operator-(const Field3<T>& other) const {
        assert(this->N1_ == other.N1_ && this->N2_ == other.N2_ && this->N3_ == other.N3_);
        Field3<T> result = *this;
        result -= other;        
        result.N1_ = this->N1_;
        result.N2_ = this->N2_;
        result.N3_ = this->N3_;
        return result;
    }
    
    Field3<T> operator-() const {        
        Field3<T> result(this->N1_, this->N2_, this->N3_);
        result.setZero(); 
        result -= *this;        
        return result;
    }

    Field3<T> operator*(const T& scalar) const {        
        Field3<T> result = *this;
        result *= scalar;
        return result;
    }   

    void resize3D(const int& N1, const int& N2, const int& N3) {
        N1_ = N1;
        N2_ = N2;
        N3_ = N3;        
        this->resize(N1 * N2 * N3);
    }

    int getN1() const {return N1_;}
    int getN2() const {return N2_;}
    int getN3() const {return N3_;}    
};

template <class T>
Field3<T> operator*(const T& scalar, const Field3<T>& other) {
    Field3<T> result = other;
    result *= scalar;
    return result;
}

// Alias
using Field3d = Field3<double>;
using Field3i = Field3<int>;

}

#endif