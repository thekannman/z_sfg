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

// Implementation of ParticleGroup class. See include/z_atom_group.hpp for
// more details about the class.

#include "z_particle_group.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "z_string.hpp"
#include "z_vec.hpp"
#include <algorithm>

void ParticleGroup::CalculateElectricField(const ParticleGroup& other_group,
                              const arma::rowvec& box,
                              const arma::imat& nearby_molecules,
                              arma::rowvec& dx) {
  for (int i_other = 0; i_other < other_group.size(); ++i_other) {
    if (included_in_field_[other_group.indices(i_other)])
      continue;
    else
      included_in_field_[other_group.indices(i_other)] = true;
    for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
      if (nearby_molecules(i_atom, other_group.index_to_molecule(i_other)) == 0)
        continue;
      //Skip neutral atoms
      if (std::abs(other_group.charge(i_other) < 0.0001)) continue;
      FindDxNoShift(dx, position(i_atom), other_group.position(i_other), box);
      double r2 = arma::dot(dx,dx);
      electric_fields_.row(i_atom) +=
          other_group.charge(i_other)*dx/r2/std::sqrt(r2);
    }
  }
}

void ParticleGroup::CalculateElectricField(const ParticleGroup& other_group,
                              const arma::rowvec& box,
                              const double cutoff_squared,
                              arma::rowvec& dx) {
  for (int i_other = 0; i_other < other_group.size(); ++i_other) {
    if (included_in_field_[other_group.indices(i_other)])
      continue;
    else
      included_in_field_[other_group.indices(i_other)] = true;
    for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
      // Ignore contribution from same molecule
      if (other_group.index_to_molecule(i_other) == index_to_molecule(i_other))
        continue;
      //Skip neutral atoms
      if (std::abs(other_group.charge(i_other) < 0.0001)) continue;
      FindDxNoShift(dx, position(i_atom), other_group.position(i_other), box);
      double r2 = arma::dot(dx,dx);
      if (r2 > cutoff_squared) continue;
      electric_fields_.row(i_atom) +=
          other_group.charge(i_other)*dx/r2/std::sqrt(r2);
    }
  }
}

void ParticleGroup::MarkNearbyAtoms(const ParticleGroup& other_group,
                                const arma::rowvec& box,
                                const double cutoff_squared,
                                arma::imat& nearby_atoms) const {
  nearby_atoms = arma::zeros<arma::imat>(group_size_,
                                         other_group.num_molecules());
  arma::rowvec dx;
  for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
    const int molecule_number = index_to_molecule_[i_atom];
    for (std::vector<int>::const_iterator i_other =
             other_group.molecule_position_markers_.begin();
         i_other != other_group.molecule_position_markers_.end(); ++i_other) {
      if (molecule_number == other_group.index_to_molecule(*i_other)) continue;
      FindDxNoShift(dx, position(i_atom), other_group.position(*i_other), box);
      double r2 = arma::dot(dx,dx);
      if (r2 < cutoff_squared)
        nearby_atoms(i_atom,*i_other) = 1;
    }
  }
}

void ParticleGroup::FindClusters(const double cutoff_distance_squared,
                             const arma::rowvec& box,
                             Histogram& clusters) const {
  arma::irowvec in_current_cluster;
  arma::irowvec in_any_cluster = arma::zeros<arma::irowvec>(group_size_);
  arma::rowvec dx;
  for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
    if (in_any_cluster(i_atom)) continue;
    in_current_cluster.zeros();
    in_current_cluster(i_atom) = 1;
    in_any_cluster(i_atom) = 1;
    for (int i_other = 0; i_other < group_size_; ++i_other) {
      if (in_any_cluster(i_other)) continue;
      FindDxNoShift(dx, position(i_atom), position(i_other), box);
      double r2 = arma::dot(dx,dx);
      if (r2 > cutoff_distance_squared) continue;
      in_current_cluster(i_other) = 1;
      in_any_cluster(i_other) = 1;
      FindCluster(cutoff_distance_squared, box, dx, in_current_cluster,
                  in_any_cluster, i_other);
    }
    clusters.Add(arma::sum(in_current_cluster));
  }
}

void ParticleGroup::FindNearestk(const arma::rowvec& point, const arma::rowvec& box,
                             const int exclude_molecule, const int k,
                             arma::rowvec& dx, arma::rowvec& r2,
                             arma::irowvec& nearest) const {
  for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
    if (index_to_molecule_[i_atom] == exclude_molecule) {
      r2(i_atom) = std::numeric_limits<double>::max();
    } else {
      FindDxNoShift(dx, point, positions_.row(i_atom), box);
      r2(i_atom) = arma::dot(dx,dx);
    }
  }
  arma::uword min_index;
  for (int i = 0; i < k; ++i) {
    r2.min(min_index);
    r2(min_index) = std::numeric_limits<double>::max();
    nearest(i) = min_index;
  }
}

void ParticleGroup::FindCluster(const double cutoff_distance_squared,
                             const arma::rowvec& box, arma::rowvec& dx,
                             arma::irowvec& in_current_cluster,
                             arma::irowvec& in_any_cluster,
                             const int start_atom) const {
  for (int i_other = start_atom + 1; i_other < group_size_; ++i_other) {
    if (in_any_cluster(i_other)) continue;
    FindDxNoShift(dx, position(start_atom), position(i_other), box);
    double r2 = arma::dot(dx,dx);
    if (r2 > cutoff_distance_squared) continue;
    in_current_cluster(i_other) = 1;
    in_any_cluster(i_other) = 1;
    FindCluster(cutoff_distance_squared, box, dx, in_current_cluster,
                  in_any_cluster, i_other);
  }
}
