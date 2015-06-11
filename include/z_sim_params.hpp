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

// Provides convenient access to simulation parameters as well as the ability
// to read those parameters from an intuitive parameter file.

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include "xdrfile_xtc.h"
#include "z_constants.hpp"

#ifndef _Z_SIMPARAMS_HPP_
#define _Z_SIMPARAMS_HPP_

enum CombRule {kC6C12, kLB, kGeometric};

class SimParams {
 public:
  // Wrapper for ReadParams
  inline void ReadParams(const std::string& filename) {
    set_filename(filename);
    ReadParams();
  }

  // Extracts box and timestep data from xtc file
  void ExtractTrajMetadata(char *traj, rvec **x_in, arma::rowvec& box);


  // Basically a mutator
  inline void set_temperature(const double temperature) {
    temperature_ = temperature;
    kT_ = KB*temperature;
    beta_ = 1.0/kT_;
  }

  // Basically a mutator
  inline void set_max_time(const double max_time) {
    if (max_time < 1.0e-5)
      max_steps_ = std::numeric_limits<int>::max();
    else
      max_steps_ = max_time / dt_;
  }

  // Mutators
  inline void set_filename(const std::string& filename) {
    filename_ = filename;
  }

  inline void set_comb_rule(const CombRule& comb_rule) {
    comb_rule_ = comb_rule;
  }

  inline void set_box(arma::rowvec box) { box_ = box; }

  inline void set_box(matrix box_mat) {
    for(int i=0; i<DIMS; i++)
      box_(i) = box_mat[i][i];
  }

  // Basically an accessor
  inline int max_time() const { return (max_steps_*dt_); }
  inline double volume() const { return box_(0)*box_(1)*box(2); }

  // Accessors
  inline arma::rowvec box() const { return box_; }
  inline double box(int i) const { return box_(i); }
  inline double temperature() const { return temperature_;}
  inline double kT() const { return kT_;}
  inline double beta() const { return beta_;}
  inline double dt() const { return dt_; }
  inline int steps() const { return steps_; }
  inline int num_atoms() const { return num_atoms_; }
  inline int max_steps() const { return max_steps_; }

 private:
  // Uses parameter file to initialize class members.
  void ReadParams();

  CombRule comb_rule_;
  arma::rowvec box_;
  std::string filename_;
  std::ifstream file_;
  double temperature_;
  double kT_;
  double beta_;
  double dt_;
  int steps_;
  int num_atoms_;
  int max_steps_;
  int numMols_, numChromos_, wAtoms_;
  int numCations_, numAnions_;
  int numIons_, numOthers_;
  double rminOO_, rminOO2_;
  double rminC_, rminA_, rminAH_, rminC2_, rminA2_, rminAH2_;
  double rmin2C_, rmin2A_, rmin2C2_, rmin2A2_;
  double rmin3C_, rmin3A_, rmin3C2_, rmin3A2_;
  double gamm_, gammC_, gammA_;
};
#endif
