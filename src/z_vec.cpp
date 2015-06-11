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

#include "z_vec.hpp"

void SetupDx(arma::cube& xoo, arma::cube& xio, arma::cube& xii, arma::mat& roo2, arma::mat& rio2, arma::mat& rii2, const arma::mat& x, const arma::mat& xion, arma::icube& shiftoo,
                    arma::icube& shiftio, arma::icube& shiftii, const int numMols, const int numIons, const arma::rowvec& box)
{
    arma::rowvec x_row(DIMS);
    arma::irowvec shift_row(DIMS);

    for (int i1=0; i1<numMols; i1++)
        for (int i2=i1+1; i2<numMols; i2++)
        {
            FindDx(x_row, x.row(i1), x.row(i2), box, shift_row);
            xoo.tube(i1,i2) = x_row;
            shiftoo.tube(i1,i2) = shift_row;
            xoo.tube(i2,i1) = -1.0*x_row;
            shiftoo.tube(i2,i1) = -1*shift_row;
            roo2(i2,i1) = roo2(i1,i2) = dot(x_row,x_row);
        }
    for (int i1=0; i1<numIons; i1++)
    {
        for (int i2=0; i2<numMols; i2++)
        {
            FindDx(x_row, xion.row(i1), x.row(i2), box, shift_row);
            xio.tube(i1,i2) = x_row;
            shiftio.tube(i1,i2) = shift_row;
            rio2(i2,i1) = rio2(i1,i2) = dot(x_row,x_row);
        }
        for (int i2=i1+1; i2<numIons; i2++)
        {
            FindDx(x_row, xion.row(i1), xion.row(i2), box, shift_row);
            xii.tube(i1,i2) = x_row;
            shiftii.tube(i1,i2) = shift_row;
            xii.tube(i2,i1) = -1.0*x_row;
            shiftii.tube(i2,i1) = -1*shift_row;
            rii2(i2,i1) = rii2(i1,i2) = dot(x_row,x_row);
        }
    }
}

void SetupDxIO(arma::cube& xio, const arma::mat& x, const arma::mat& xion, arma::icube& shiftio, const int numMols, const int numIons, const arma::rowvec& box)
{
    arma::rowvec x_row(DIMS);
    arma::irowvec shift_row(DIMS);

    for (int i1=0; i1<numIons; i1++)
        for (int i2=0; i2<numMols; i2++)
        {
            FindDx(x_row, xion.row(i1), x.row(i2), box, shift_row);
            xio.tube(i1,i2) = x_row;
            shiftio.tube(i1,i2) = shift_row;
        }
}

void SetupShiftOO(const arma::mat& x, arma::icube& shiftoo, const int numMols, const arma::rowvec& box)
{
    arma::rowvec x_row(DIMS);
    arma::irowvec shift_row(DIMS);

    for (int i1=0; i1<numMols; i1++)
        for (int i2=i1+1; i2<numMols; i2++)
        {
            FindDx(x_row, x.row(i1), x.row(i2), box, shift_row);
            shiftoo.tube(i1,i2) = shift_row;
            shiftoo.tube(i2,i1) = -1*shift_row;
        }
}

void FindDxNoShift(arma::rowvec& result, const arma::rowvec& vec1, const arma::rowvec& vec2, const arma::rowvec& box)
{
    result = vec1 - vec2;
    for (int i1=0; i1<DIMS; i1++)
    {
        if (result(i1)>box(i1)/2.0)
            result(i1) -= box(i1);
        else if (result(i1)<-box(i1)/2.0)
            result(i1) += box(i1);
    }
}

void FindDx(arma::rowvec& result, const arma::rowvec& vec1, const arma::rowvec& vec2, const arma::rowvec& box, arma::irowvec& shift)
{
    result = vec1 - vec2;
    for (int i1=0; i1<DIMS; i1++)
    {
        if (result(i1)>box(i1)/2.0)
        {
            result(i1) -= box(i1);
            shift(i1) = -1;
        }
        else if (result(i1)<-box(i1)/2.0)
        {
            result(i1) += box(i1);
            shift(i1) = 1;
        }
        else shift(i1) = 0;
    }
}
