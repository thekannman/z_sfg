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

#include "z_frequency_map.hpp"
#include "z_vec.hpp"

void PrintSpectra(const arma::cx_rowvec& spectra, const std::string& filename,
                  const int corr, const int nzeros, const double deltaT,
                  const double avgFreq) {
  std::ofstream file;
  file.open(filename.c_str());
  double factor = 1.0/2.0/C_SPEED/(deltaT*static_cast<double>(corr+nzeros));

  int imax = 2*(corr+nzeros)-2;
  for (int i=corr+nzeros; i<imax; i++) {
    file << std::fixed <<
            static_cast<double>(i+2-2*(corr+nzeros))*factor+avgFreq << " ";
    file << std::scientific << real(spectra(i)) << std::endl;
  }
  imax = corr+nzeros;
  for (int i=0; i<imax; i++) {
    file << std::fixed << static_cast<double>(i)*factor + avgFreq << " ";
    file << std::scientific << real(spectra(i)) << std::endl;
  }
  file.close();
}
