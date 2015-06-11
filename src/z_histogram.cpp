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

// Implementation of Hist class. See include/z_hist.hpp for
// more details about the class.

#include "z_histogram.hpp"

void Histogram::Print(const bool normalize) {
  file_.open(filename_.c_str());
  for (int i=0; i<number_of_bins_; i++) {
    double normed;
    if (normalize && (counts_(i) > 0))
      normed = array_(i)/counts_(i);
    else
      normed = array_(i);
    file_ << std::fixed <<
             min_ + bin_width_*(static_cast<double>(i)+0.5) << " ";
    file_ << std::scientific << normed << std::endl;
  }
  file_.close();
}

void Histogram::PrintFraction() {
  file_.open(filename_.c_str());
  double sum = arma::sum(array_);
  for (int i=0; i<number_of_bins_; i++) {
    file_ << std::fixed <<
             min_ + bin_width_*(static_cast<double>(i)+0.5) << " ";
    file_ << std::scientific << array_(i)/sum << std::endl;
  }
  file_.close();
}
