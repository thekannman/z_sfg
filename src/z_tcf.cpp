//Copyright (c) 2015 Zachary Kann
//
//Permission is hereby granted, free of charge, to any person obtaining a copy
//of this software and associated documentation files (the "Software"), to deal
//in the Software without restriction, including without limitation the rights
//to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//copies of the Software, and to permit persons to whom the Software is
//furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all
//copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.

// ---
// Author: Zachary Kann

// Implementation of TCF class. See include/z_tcf.hpp for
// more details about the class.

#include "z_tcf.hpp"

template <class T>
void TCF<T>::Correlate(const arma::Mat<T>& mat1, const arma::Mat<T>& mat2) {
  int num_corr = CalculateNumCorr(mat1,mat2);
  for (int i1=0; i1<num_corr; i1++) {
    int i1corr = i1*interval_;
    correlation_function_(0) +=
        2.0*arma::dot(mat1.row(i1corr),mat2.row(i1corr))/
        static_cast<T>(num_corr);
    for (int i2 = 1; i2 < length_; i2++) {
      int i1corri2 = i1corr + i2;
      correlation_function_(i2) +=
          arma::dot(mat1.row(i1corr),mat2.row(i1corri2))/
          static_cast<T>(num_corr);
      correlation_function_(i2) +=
          arma::dot(mat1.row(i1corri2),mat2.row(i1corr))/
          static_cast<T>(num_corr);
    }
  }
}

template <class T>
void TCF<T>::Correlate(const arma::Mat<T>& mat1, const arma::Mat<T>& mat2,
                       const int mod) {
  int corrmod = mod + length_;
  for (int i2=mod; i2<corrmod; i2++) {
    int i2min = i2-mod;
    int i2mod = i2%length_;
    correlation_function_(i2min) += arma::dot(mat1.row(mod),mat2.row(i2mod));
    correlation_function_(i2min) += arma::dot(mat1.row(i2mod),mat2.row(mod));
  }
}

template <class T>
void TCF<T>::CorrelateOneDirection(const arma::Mat<T>& mat1,
                                   const arma::Mat<T>& mat2) {
  int num_corr = CalculateNumCorr(mat1,mat2);
  for (int i1=0; i1<num_corr; i1++) {
    int i1corr = i1*interval_;
    for (int i2=0; i2<length_; i2++) {
      int i1corri2 = i1corr + i2;
      correlation_function_(i2) +=
          arma::dot(mat1.row(i1corr),mat2.row(i1corri2))/
          static_cast<T>(num_corr);
    }
  }
}

template <class T>
void TCF<T>::CorrelateOneDirection(const arma::Mat<T>& mat1,
                                   const arma::Mat<T>& mat2, const int mod) {
  int corrmod = mod + length_;
  for (int i2=mod; i2<corrmod; i2++) {
    int i2min = i2-mod;
    int i2mod = i2%length_;
    correlation_function_(i2min) += arma::dot(mat1.row(mod),mat2.row(i2mod));
  }
}

template <class T>
void TCF<T>::Correlate(const arma::Row<T>& vec1, const arma::Row<T>& vec2) {
  int num_corr = CalculateNumCorr(vec1,vec2);
  for (int i1=0; i1<num_corr; i1++) {
    int i1corr = i1*interval_;
    correlation_function_(0) +=
        2.0*vec1[i1corr]*vec2[i1corr]/static_cast<T>(num_corr);;
    for (int i2=1; i2<length_; i2++) {
      int i1corri2 = i1corr + i2;
      correlation_function_(i2) +=
          vec1[i1corr]*vec2[i1corri2]/static_cast<T>(num_corr);
      correlation_function_(i2) +=
          vec1[i1corri2]*vec2[i1corr]/static_cast<T>(num_corr);
    }
  }
}

template <class T>
void TCF<T>::Correlate(const arma::Row<T>& vec1, const arma::Row<T>& vec2,
  const int mod) {
  int corrmod = mod + length_;
  for (int i2=mod; i2<corrmod; i2++) {
    int i2min = i2-mod;
    int i2mod = i2%length_;
    correlation_function_(i2min) += vec1[mod]*vec2[i2mod];
    correlation_function_(i2min) += vec1[i2mod]*vec2[mod];
  }
}

template <class T>
void TCF<T>::CorrelateOneDirection(const arma::Row<T>& vec1,
                                   const arma::Row<T>& vec2) {
  int num_corr = CalculateNumCorr(vec1,vec2);
  for (int i1=0; i1<num_corr; i1++) {
    int i1corr = i1*interval_;
    for (int i2=0; i2<length_; i2++) {
      int i1corri2 = i1corr + i2;
      correlation_function_(i2) +=
          vec1[i1corr]*vec2[i1corri2]/static_cast<T>(num_corr);
    }
  }
}

template <class T>
void TCF<T>::CorrelateOneDirection(const arma::Row<T>& vec1,
                                   const arma::Row<T>& vec2, const int mod) {
  int corrmod = mod + length_;
  for (int i2=mod; i2<corrmod; i2++) {
    int i2min = i2-mod;
    int i2mod = i2%length_;
    correlation_function_(i2min) += vec1[mod]*vec2[i2mod];
  }
}

// Added to avoid linker error.
template class TCF<double>;
template class TCF<arma::cx_double>;
template class TCF<int>;
