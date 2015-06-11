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

// This class is for dealing with complex time-correlation functions (TCFs)
// and their fourier transforms.

#include <armadillo>
#include <cassert>
#include "z_tcf.hpp"
#include "z_constants.hpp"

#ifndef _Z_CX_TCF_HPP_
#define _Z_CX_TCF_HPP_

class CxTCF:
  public TCF<arma::cx_double> {
 public:
  inline CxTCF(int length, int x_spacing, int interval = 1, int zeros = 0)
      : TCF(length, x_spacing, interval, zeros),
        spec_length_(2*total_length()-2) {
    // TODO(Zak): I think this should be -1 not -2... test that.
    correlation_function_.zeros(spec_length_);
  }

  // The standard call to prep and run a fourier transform. Assumes symmetry
  // around t = 0.
  inline void FourierPlus() {
    for (int i = 0; i<length(); ++i) {
      correlation_function_(i) -= correlation_function_(length()-1);
    }
    for (int i = length()+2*zeros()-1; i<spec_length_; i++) {
      correlation_function_(i) =
        std::conj(correlation_function_(spec_length_-i));
    }
    fourier_transform_ = fft(correlation_function_);
    freq_spacing_nu_ = 1.0/2.0/C_SPEED/(x_spacing()*total_length());
    freq_spacing_omega_ = freq_spacing_nu_*2.0*M_PI*C_SPEED;
  }

  // Similar to FourierPlus, but without the assumption of symmetry
  // around t = 0.
  inline void FourierNoSymmPlus() {
    for (int i = 0; i<length(); ++i) {
      correlation_function_(i) -= correlation_function_(length()-1);
    }
    fourier_transform_ = fft(correlation_function_);
    freq_spacing_nu_ = 1.0/2.0/C_SPEED/(x_spacing()*total_length());
    freq_spacing_omega_ = freq_spacing_nu_*2.0*M_PI*C_SPEED;
  }

  // Similar to FourierPlus, but doesn't shift the time correlation function
  // to be zero after corr steps. Formerly called run_fftw.
  inline void Fourier() {
    for (int i = length()+2*zeros()-1; i<spec_length_; i++) {
      correlation_function_(i) =
        std::conj(correlation_function_(spec_length_-i));
    }
    fourier_transform_ = fft(correlation_function_);
    freq_spacing_nu_ = 1.0/2.0/C_SPEED/(x_spacing()*total_length());
    freq_spacing_omega_ = freq_spacing_nu_*2.0*M_PI*C_SPEED;
  }

  inline void MultiplyFT(const double multiplier) {
    fourier_transform_ *= multiplier;
  }

  inline void FreqWeighting(const double power) {
    for (int i = 0; i < spec_length_; ++i) {
      fourier_transform_(i) *= std::pow(i*freq_spacing_omega_, power);
    }
  }

  inline void PrintFTReal(const std::string& filename, double x_shift = 0.0) {
    std::ofstream output_file(filename.c_str());
    assert(output_file.is_open());
    for (int i = total_length(); i < spec_length_; ++i) {
      output_file << (i-spec_length_)*freq_spacing_nu_ + x_shift << " " <<
                     real(fourier_transform_(i)) << std::endl;
    }
    for (int i = 0; i < total_length(); ++i) {
      output_file << i*freq_spacing_nu_ + x_shift << " " <<
                     real(fourier_transform_(i)) << std::endl;
    }
    output_file.close();
  }

  // Mutators
  inline void set_ft(const int i, const arma::cx_double data) {
    fourier_transform_(i) = data;
  }

  // Accessors
  inline int spec_length() const { return spec_length_; }
  inline double freq_spacing_nu() const { return freq_spacing_nu_; }
  inline double freq_spacing_omega() const { return freq_spacing_omega_; }
  inline arma::cx_double ft(int i) const { return fourier_transform_(i); }

 private:
  arma::cx_rowvec fourier_transform_;
  int spec_length_;
  double freq_spacing_nu_, freq_spacing_omega_;
};
#endif
