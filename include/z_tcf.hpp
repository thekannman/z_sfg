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

// This class is for dealing with time-correlation functions (TCFs) between
// scalars or vectors.

#include <armadillo>

#ifndef _Z_TCF_HPP_
#define _Z_TCF_HPP_

template <class T>
class TCF {
 public:
  TCF(int length, int x_spacing, int interval = 1, int zeros = 0)
      : length_(length), x_spacing_(x_spacing), interval_(interval),
        zeros_(zeros) {
    // TODO(Zak): I think this should be -1 not -2... test that.
    correlation_function_.zeros(length_+zeros_);
  }

  // Calculations correlation between two vector time series assuming
  // time-symmetry around t=0.
  void Correlate(const arma::Mat<T>& mat1, const arma::Mat<T>& mat2);

  // Calculations correlation between two vector time series assuming
  // time-symmetry around t=0. Allows time series to start at arbitrary index
  // and wrap around end of vector.
  void Correlate(const arma::Mat<T>& mat1, const arma::Mat<T>& mat2, const int mod);

  // Calculations correlation between two vector time series without assuming
  // time-symmetry around t=0.
  void CorrelateOneDirection(const arma::Mat<T>& mat1,
                             const arma::Mat<T>& mat2);

  // Calculations correlation between two vector time series without assuming
  // time-symmetry around t=0. Allows time series to start at arbitrary index
  // and wrap around end of vector.
  void CorrelateOneDirection(const arma::Mat<T>& mat1, const arma::Mat<T>& mat2,
                               const int mod);

  // Wrapper for Correlate to avoid ambiguity
  inline void Correlate(const arma::subview_cube<T>& vec1,
                        const arma::subview_cube<T>& vec2) {
    Correlate(arma::Row<T>(vec1), arma::Row<T>(vec2));
  }

  // Wrapper for Correlate to avoid ambiguity
  inline void Correlate(const arma::subview_cube<T>& vec1,
                        const arma::subview_cube<T>& vec2, const int mod) {
    Correlate(arma::Row<T>(vec1), arma::Row<T>(vec2), mod);
  }

  // Wrapper for CorrelateOneDirection to avoid ambiguity
  inline void CorrelateOneDirection(const arma::subview_cube<T>& vec1,
                                    const arma::subview_cube<T>& vec2) {
    CorrelateOneDirection(arma::Row<T>(vec1), arma::Row<T>(vec2));
  }

  // Wrapper for CorrelateOneDirection to avoid ambiguity
  inline void CorrelateOneDirection(const arma::subview_cube<T>& vec1,
                                    const arma::subview_cube<T>& vec2,
                                    const int mod) {
    CorrelateOneDirection(arma::Row<T>(vec1), arma::Row<T>(vec2), mod);
  }

  // Calculations correlation between two scalar time series assuming
  // time-symmetry around t=0.
  void Correlate(const arma::Row<T>& vec1, const arma::Row<T>& vec2);

  // Calculations correlation between two scalar time series assuming
  // time-symmetry around t=0. Allows time series to start at arbitrary index
  // and wrap around end of vector.
  void Correlate(const arma::Row<T>& vec1, const arma::Row<T>& vec2,
                    const int mod);

  // Calculations correlation between two scalar time series without assuming
  // time-symmetry around t=0.
  void CorrelateOneDirection(const arma::Row<T>& vec1,
                             const arma::Row<T>& vec2);

  // Calculations correlation between two scalar time series without assuming
  // time-symmetry around t=0. Allows time series to start at arbitrary index
  // and wrap around end of vector.
  void CorrelateOneDirection(const arma::Row<T>& vec1,
                                  const arma::Row<T>& vec2, const int mod);

  // Wrapper for CorrelateOneDirection that allows vector self-correlation.
  inline void CorrelateOneDirection(const arma::Mat<T>& mat) {
    CorrelateOneDirection(mat, mat);
  }

  // Wrapper for CorrelateOneDirection that allows vector self-correlation and
  // allows time series to start at arbitrary index and wrap around end of
  // vector.
  inline void CorrelateOneDirection(const arma::Mat<T>& vec, const int mod) {
    CorrelateOneDirection(vec, vec, mod);
  }

  // Wrapper for CorrelateOneDirection that allows scalar self-correlation.
  inline void CorrelateOneDirection(const arma::Row<T>& vec) {
    CorrelateOneDirection(vec, vec);
  }

  // Wrapper for CorrelateOneDirection that allows scalar self-correlation and
  // allows time series to start at arbitrary index and wrap around end of
  // vector.
  inline void CorrelateOneDirection(const arma::Row<T>& vec,
                                         const int mod) {
    CorrelateOneDirection(vec, vec, mod);
  }

  inline int CalculateNumCorr(const arma::Row<T>& vec) {
    return (vec.n_elem - length_)/interval_ + 1;
  }

  inline int CalculateNumCorr(const int steps) {
    return (steps - length_)/interval_ + 1;
  }

  inline void MultiplyTCF(const double multiplier) {
    correlation_function_ *= multiplier;
  }

  inline void MultiplyTCF(const int i, const T multiplier) {
    correlation_function_(i) *= multiplier;
  }

  inline void AddTCF(const int i, const T addend) {
    correlation_function_(i) += addend;
  }

  // Mutators
  inline void set_tcf(const int i, const T data) {
    correlation_function_(i) = data;
  }

  // Accessors
  int length() const { return length_; }
  int x_spacing() const { return x_spacing_; }
  int zeros() const { return zeros_; }
  int total_length() const { return length_+zeros_; }
  int interval() const { return interval_; }
  inline T tcf(int i) const { return correlation_function_(i); }
  inline arma::Row<T> tcf() const { return correlation_function_; }

 protected:
  arma::Row<T> correlation_function_;

 private:
  int length_, x_spacing_, interval_, zeros_, number_;

  inline int CalculateNumCorr(const arma::Row<T>& vec1,
                              const arma::Row<T>& vec2) {
    int n1 = vec1.n_elem;
    int n2 = vec2.n_elem;
    int shorter = (n1 < n2) ? n1 : n2;
    return (shorter - length_)/interval_ + 1;
  }

  inline int CalculateNumCorr(const arma::Mat<T>& mat) {
    return (mat.n_rows - length_)/interval_ + 1;
  }

  inline int CalculateNumCorr(const arma::Mat<T>& mat1,
                              const arma::Mat<T>& mat2) {
    int n1 = mat1.n_rows;
    int n2 = mat2.n_rows;
    int shorter = (n1 < n2) ? n1 : n2;
    return (shorter - length_)/interval_ + 1;
  }
};
#endif
