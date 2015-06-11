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

#include <complex>
#include <cmath>
#include <armadillo>
#include "z_histogram.hpp"
#include "z_atom_group.hpp"
#include "z_molecule_group.hpp"
#include "z_cx_tcf.hpp"
#include "z_conversions.hpp"
#ifndef _Z_FREQUENCY_MAP_HPP_
#define _Z_FREQUENCY_MAP_HPP_

enum Chromophore {kOH, kOD};
enum ModelType {kThreeSite, kFourSite};

class FrequencyMap {
 public:
  inline FrequencyMap(const MoleculeGroup& water_group,
                      const Chromophore chromophore, const double timestep,
                      const int correlation_length = 200)
      : step_(0), correlation_length_(correlation_length),
        num_chromophores_(2*water_group.num_molecules()),
        chromophore_(chromophore), timestep_(timestep), frequency_bins_(4000),
        min_frequency_(chromophore == kOH ? 1650.0 : 850.0),
        max_frequency_(chromophore == kOH ? 4950.0 : 4150.0),
        frequency_distribution_(frequency_bins_, min_frequency_,
                                max_frequency_),
        spectral_density_(frequency_bins_, min_frequency_, max_frequency_),
        spectra_tcf_(correlation_length,timestep_, 1,
                    (4096+2)/2-correlation_length) {
    steps_guess_ = 1000;
    omega_01_ = arma::zeros<arma::mat>(num_chromophores_, steps_guess_);

    switch(chromophore) {
      case kOH:
        SetupOH();
        break;
      case kOD:
        SetupOD();
        break;
      default:
        assert(false && "Unrecognized chromophore option.");
        break;
    }
    switch(water_group.size()/water_group.num_molecules()) {
      case 3:
        model_type_ = kThreeSite;
        break;
      case 4:
        model_type_ = kFourSite;
        break;
      default:
        assert(false && "Number of atoms in water model is not one 3 or 4.");
        break;
    }
  }

  inline void CalculateFrequencies(const MoleculeGroup& hydrogen_group,
                                   const AtomGroup& oxygen_group,
                                   const arma::rowvec& box) {
    if(step_ == steps_guess_) {
      steps_guess_ += kStepsGuessIncrement_;
      ResizeArrays();
    }
    for (int i_chrom = 0; i_chrom < hydrogen_group.size(); ++i_chrom) {
      int i_oxygen = i_chrom/2;
      oh_vector_ =
          hydrogen_group.position(i_chrom) - oxygen_group.position(i_oxygen);
      oh_unit_vector_ = arma::normalise(oh_vector_);
      field_projection_ =
          arma::dot(oh_vector_, hydrogen_group.electric_field(i_chrom))*
          AU_TO_NM*AU_TO_NM;
      transition_dipole_position_ = oxygen_group.position(i_oxygen) +
                                    transition_dipole_shift_*oh_vector_;
      UseMapping(i_chrom, box);
      frequency_distribution_.Add(omega_01_(i_chrom));
      double spectral_density_weight = SpectralDensityWeight(i_chrom);
      spectral_density_.Add(omega_01_(i_chrom), spectral_density_weight);
    }
    step_++;
  }

  inline void CalculateSpectra() {
    // Two mean functions needed because of armadillo convention.
    avg_frequency_ = arma::mean(arma::mean(omega_01_));
    omega_01_ -= avg_frequency_;
    int num_corr = spectra_tcf_.CalculateNumCorr(step_);
    std::complex<double> exp_argument, exp_omega;
    for (int i_chrom = 0; i_chrom < num_chromophores_; ++i_chrom) {
      for (int i_num = 0; i_num < num_corr; ++i_num) {
        int corr_start = i_num*spectra_tcf_.interval();
        double omega_integral = 0.0;
        exp_argument = -1.0i*omega_integral;
        exp_omega = std::exp(exp_argument);
        double mu_or_alpha_product =
            MuOrAlphaProduct(i_chrom, corr_start, 0);
        spectra_tcf_.AddTCF(0,mu_or_alpha_product*exp_omega);
        for (int i_corr = 1; i_corr < correlation_length_; ++i_corr) {
          omega_integral += timestep_/2.0*TWO_PI_C *
                            (omega_01_(i_chrom,corr_start+i_corr)
                             +omega_01_(i_chrom,corr_start));
          exp_argument = -1.0i*omega_integral;
          exp_omega = std::exp(exp_argument);
          mu_or_alpha_product =
              MuOrAlphaProduct(i_chrom, corr_start, i_corr);
          spectra_tcf_.AddTCF(i_corr,mu_or_alpha_product*exp_omega);
        }
      }
    }
    for (int i_corr = 0; i_corr < correlation_length_; ++i_corr)
      spectra_tcf_.MultiplyTCF(std::exp(-i_corr*timestep_/2.0/lifetime_));
    spectra_tcf_.FourierPlus();
  }

  inline void PrintSpectra(const std::string& filename) {
    spectra_tcf_.PrintFTReal(filename, avg_frequency_);
  }

  //accessors
  inline int num_chromophores() const { return num_chromophores_; }

 protected:
  int step_;
  static const double mu_prime_0_ = 0.1646;
  static const double mu_prime_1_ = 11.39;
  static const double mu_prime_2_ = 63.41;
  int steps_guess_;
  static const int kStepsGuessIncrement_ = 1000;
  arma::mat omega_01_;
  arma::cube mu_01_;
  arma::rowvec transition_dipole_position_;
  arma::rowvec oh_unit_vector_;
  double field_projection_;
  double omega_01_0_;
  double omega_01_1_;
  double omega_01_2_;
  double chi_01_0_;
  double chi_01_1_;

  inline double SpectralDensityWeight(const int chromophore) {
    return MuOrAlphaProduct(chromophore, step_, 0);
  }

 private:
  double avg_frequency_;
  double lifetime_;
  int correlation_length_;
  double num_chromophores_;
  static const double transition_dipole_shift_ = 0.067; // Formerly dipPos
  Chromophore chromophore_;
  double timestep_;
  int frequency_bins_;
  double min_frequency_;
  double max_frequency_;
  ModelType model_type_;
  arma::mat alpha_;
  arma::mat switcher_;
  //arma::cube u; // needed for coupled spectra
  arma::rowvec oh_vector_;
  Histogram frequency_distribution_;
  Histogram spectral_density_;
  CxTCF spectra_tcf_;

  virtual void ResizeArrays() = 0;

  virtual double MuOrAlphaProduct(const int chromophore, const int corr_start,
                                  const int i_corr) = 0;

  inline void SetupOH() {
    lifetime_ = 0.700; //TODO(Zak): change this for coupled spectra
    omega_01_0_ = 3760.2;
    omega_01_1_ = -3541.7;
    omega_01_2_ = -152677.0;
    chi_01_0_ = 0.19285;
    chi_01_1_ = -1.7261E-5;
  }

  inline void SetupOD() {
    lifetime_ = 1.800; //TODO(Zak): change this for coupled spectra
    omega_01_0_ = 2767.8;
    omega_01_1_ = -2630.3;
    omega_01_2_ = -102601.0;
    chi_01_0_ = 0.16593;
    chi_01_1_ = -2.0632E-5;
  }

  // Calculates omega and mu/alpha using the frequency maps.
  virtual void UseMapping(const int chromophore,
                  const arma::rowvec& box) = 0;

  // TODO(Zak): Move this to a coupled spectra class
  // Calculates the intramolecular coupling between OH groups for inclusion
  // in the system Hamiltonian needed for coupled spectra calculations
  //extern double IntraCouple();

  // TODO(Zak): Move this to a coupled spectra class
  // Calculates the intermolecular coupling between OH groups for inclusion
  // in the system Hamiltonian needed for coupled spectra calculations
  //extern void InterCouple(arma::rowvec& omega, const arma::mat& xdipole,
  //                        const arma::rowvec& chi01, const arma::rowvec& muprime,
  //                        const arma::mat& u, const arma::rowvec& box,
  //                        const arma::icube& shift, const int numMols,
  //                        const int numChromos);
};
#endif
