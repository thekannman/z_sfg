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

#include "boost/multi_array.hpp"
#include <armadillo>
#include "z_constants.hpp"

#ifndef _Z_VEC_HPP_
#define _Z_VEC_HPP_

typedef boost::multi_array<double, 4> hypercube;

extern void SetupDx(arma::cube& xoo, arma::cube& xio, arma::cube& xii,
                    arma::mat& roo2, arma::mat& rio2, arma::mat& rii2,
                    const arma::mat& x, const arma::mat& xion,
                    arma::icube& shiftoo, arma::icube& shiftio,
                    arma::icube& shiftii, const int numMols,
                    const int numIons, const arma::rowvec& box);

extern void SetupDxIO(arma::cube& xio, const arma::mat& x,
                       const arma::mat& xion, arma::icube& shiftio,
                       const int numMols, const int  numIons,
                       const arma::rowvec& box);

extern void SetupShiftOO(const arma::mat& x, arma::icube& shiftoo,
                          const int numMols, const arma::rowvec& box);

extern void FindDxNoShift(arma::rowvec& result, const arma::rowvec& vec1,
                           const arma::rowvec& vec2, const arma::rowvec& box);

extern void FindDx(arma::rowvec& result, const arma::rowvec& vec1,
                   const arma::rowvec& vec2, const arma::rowvec& box,
                   arma::irowvec& shift);

template <typename type1>
arma::rowvec VeclikeToRow(const type1& vec1, const int  length)
{
    arma::rowvec result(length);
    for (int i=0; i<length; i++)
    {
        result[i] = vec1[i];
    }
    return result;
}

inline arma::rowvec UseDx(const arma::rowvec& vec1, const arma::rowvec& vec2,
                          const arma::rowvec& box, const arma::irowvec& shift) {
  return vec1 - vec2 + shift%box; }

#endif
