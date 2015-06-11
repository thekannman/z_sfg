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

// Class for histogramming data

#include <fstream>
#include <string>
#include <armadillo>

#ifndef _Z_HISTOGRAM_HPP_
#define _Z_HISTOGRAM_HPP_

class Histogram {
 public:
  inline Histogram(const int number_of_bins, const double bin_width,
                   const bool centered = false)
      : number_of_bins_(number_of_bins), bin_width_(bin_width),
        centered_(centered) {
    array_ = arma::zeros<arma::rowvec>(number_of_bins_);
    counts_ = arma::zeros<arma::irowvec>(number_of_bins_);
    if (centered)
      min_ = -static_cast<double>(number_of_bins/2)*bin_width_;
    else
      min_ = 0.0;
  }

  inline Histogram(const int number_of_bins, const double min_value,
                   const double max_value)
      : number_of_bins_(number_of_bins), min_(min_value),
        centered_(false) {
    array_ = arma::zeros<arma::rowvec>(number_of_bins_);
    counts_ = arma::zeros<arma::irowvec>(number_of_bins_);
    bin_width_ = (max_value - min_value)/number_of_bins;
  }

  inline Histogram(const double bin_width, const double min_value,
                   const double max_value)
      : bin_width_(bin_width), min_(min_value),
        centered_(false) {
    number_of_bins_ = (max_value - min_value) / bin_width_;
    array_ = arma::zeros<arma::rowvec>(number_of_bins_);
    counts_ = arma::zeros<arma::irowvec>(number_of_bins_);
    min_ = min_value;
  }

  inline Histogram();

  // A wrapper for Print that first sets the filename.
  inline void Print(const std::string& filename, const bool normalize = false) {
    filename_ = filename;
    Print(normalize);
  }

  // Adds a a (optionally weighted) value
  inline void Add(const double x_value, const double weight = 1.0,
                  const bool increment = true) {
    int bin_number = FindBinNumber(x_value);
    if (bin_number < 0 || bin_number >= number_of_bins_) return;
    array_(bin_number) += weight;
    if (increment)
      counts_(bin_number)++;
  }

  inline void Multiply(const double rhs) {
    array_ *= rhs;
  }

  // Weights histogram by volume of spherical shell
  inline void WeightByShellVolume() {
    double sphere_factor = 4.0*M_PI/3.0*bin_width_*bin_width_*bin_width_;
    for (int i = 0; i < number_of_bins_; ++i) {
      array_(i) *= 1.0/sphere_factor/(((i+1)*(i+1)*(i+1)-i*i*i));
    }
  }

  inline void PrintFraction(const std::string& filename) {
    filename_ = filename;
    PrintFraction();
  }

 private:
  inline int FindBinNumber (const double x_value) {
    int bin_number;
    if (centered_) {
      bin_number =
          static_cast<int>(floor(x_value/bin_width_)) + number_of_bins_/2;
    } else {
      bin_number = static_cast<int>((x_value-min_)*bin_width_);
    }
    return bin_number;
  }

  // Prints the histogram to the preset file.
  void Print(const bool normalize);

  void PrintFraction();

  int number_of_bins_;
  int total_;
  double bin_width_;
  double min_;
  arma::rowvec array_;
  arma::irowvec counts_;
  std::string filename_;
  std::ofstream file_;
  bool centered_;


};

#endif
