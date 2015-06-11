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

// A set of functions for calculation of electric field / frequency maps
// as outlined in the papers of Skinner et al. See, for example, S.M. Gruenbaum,
// et al. J. Chem. Theory Comput. 2013 9 (7), 3109. These maps are necessary
// for calculation of vibrational spectra using the methodology  outline in the
// aforementioned papers.

// Abbreviations:
//   IR = infrared
//   SFG = sum-frequency generation

#include <armadillo>
#include "z_frequency_map.hpp"
#include "z_histogram.hpp"
#include "z_molecule_group.hpp"
#include "z_cx_tcf.hpp"

#ifndef _Z_SFG_MAP_HPP_
#define _Z_SFG_MAP_HPP_

class SfgMap:
  public FrequencyMap {
 public:
  inline SfgMap(const MoleculeGroup& group, const Chromophore chromophore,
                const double timestep, const int correlation_length = 200,
                const double slab_center_z = 0.0)
      : FrequencyMap(group, chromophore, timestep, correlation_length),
        slab_center_z_(slab_center_z) {
    mu_01_ = arma::zeros<arma::mat>(num_chromophores(), steps_guess_);
    alpha_01_ = arma::zeros<arma::mat>(num_chromophores(), steps_guess_);
  }

 private:
  arma::mat mu_01_;
  arma::mat alpha_01_;
  double slab_center_z_;
  static const double r_switch_ = 0.4;
  static const double r_switch_squared_ = 0.16;
  static const double r_switch_cubed_ = 0.064;

  inline void ResizeArrays() {
    omega_01_.resize(num_chromophores(), steps_guess_);
    mu_01_.resize(num_chromophores(), steps_guess_);
    alpha_01_.resize(num_chromophores(), steps_guess_);
  }

  inline double MuOrAlphaProduct(const int chromophore, const int corr_start,
                                 const int i_corr) {
    return mu_01_(chromophore,corr_start)*mu_01_(chromophore,corr_start+i_corr);
  }

  // Calculates switching function needed for SFG calculations in slab geometry.
  // Allows smooth switching between contribution to upper and lower surface as
  // one moves through the slab.
  double SwitchingFunction(const arma::rowvec& box);

  // Calculates omega and mu/alpha using the frequency maps.
  void UseMapping(const int chromophore,
                  const arma::rowvec& box);
};
#endif
