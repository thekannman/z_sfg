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

#include "z_sfg_map.hpp"
#include "z_vec.hpp"

double SfgMap::SwitchingFunction(const arma::rowvec& box) {
  double z_shift = transition_dipole_position_(2) - slab_center_z_;
  if (z_shift > box(2)/2.0)
    z_shift -= box(2)/2.0;
  if (z_shift < -box(2)/2.0)
    z_shift += box(2)/2.0;

  if (z_shift > r_switch_)
    return 1.0;
  else if (z_shift < -r_switch_)
    return -1.0;
  else
    return 2.0*((3.0*r_switch_squared_*z_shift)-(z_shift*z_shift*z_shift))/
           (4.0*r_switch_cubed_);
}

void SfgMap::UseMapping(const int chromophore,
                        const arma::rowvec& box) {
  double chi_01, mu_prime;
  // The maps were fit to pure water simulations. When tested with ions, it
  // was found that the quadratic fit becomes problematic since chromophores
  // can experience fields more negative than the maximum of the quadratic.
  // this breaks the expected monotonocity of the frequency/field relation.
  // To resolve this issue, we use a linear extrapolation of the fit for
  // negative field values.
  if (field_projection_ >= 0) {
    omega_01_(chromophore,step_) =
        omega_01_0_ + omega_01_1_*field_projection_ +
        omega_01_2_*field_projection_*field_projection_;
    mu_prime = mu_prime_0_ + mu_prime_1_*field_projection_ +
               mu_prime_2_*field_projection_*field_projection_;
  } else {
    omega_01_(chromophore,step_) = omega_01_0_ + omega_01_1_*field_projection_;
    mu_prime = mu_prime_0_ + mu_prime_1_*field_projection_;
  }
  chi_01 = chi_01_0_ + chi_01_1_*omega_01_(chromophore,step_);
  mu_01_(chromophore,step_) =
      mu_prime*chi_01*oh_unit_vector_(2)*SwitchingFunction(box);
  alpha_01_ = chi_01*(4.6*oh_unit_vector_(0)*oh_unit_vector_(0)+1.0) +
              chi_01*(4.6*oh_unit_vector_(1)*oh_unit_vector_(1)+1.0);
}
