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

// Some useful unit conversions.

#ifndef _Z_CONVERSIONS_HPP_
#define _Z_CONVERSIONS_HPP_

const double NM_TO_AU = 18.89726125; // NM is nanometer, AU is atomic units
const double AU_TO_NM = 0.052917721092; // NM is nanometer, AU is atomic units
const double NM_TO_ANG = 10.0; // ANG is angstrom
const double ANG_TO_NM = 0.1;
const double NM_TO_M = 1.0e-9;
const double M_TO_NM = 1.0e9;
const double PS_TO_S = 1.0e-12;
const double S_TO_PS = 1.0e12;
  const double ENM_TO_D = 48.03204544; //Converts from e (charge), nm, to Debye
const double D_TO_CM = 3.33564e-30; //Converts from Debye to Coulomb, m
const double AMU_TO_KG = 1.66053886e-27;
// TODO(Zak) Replace THZ_CONVERT with ENM_TO_D in rest of files.
// TODO(Zak) Replace FIELD_CONVERT with (AU_TO_NM)^2 in rest of files.
//const double FIELD_CONVERT = 0.002800314724;
//const double THZ_CONVERT = 0.048032045;

#endif
