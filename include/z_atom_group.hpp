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

// Contains information about groups of atoms as defined in
// supplied .ndx file. This allows the properties of a subset
// of the atoms to be treated with the same ease as the entire
// set of atoms. This class includes postion, velocity, and mass
// information for individual atoms as well as for molecules
// defined in the .top file.

#include "z_subsystem_group.hpp"
#ifndef _Z_ATOM_GROUP_HPP_
#define _Z_ATOM_GROUP_HPP_

//See description above
class AtomGroup:
  public SubsystemGroup {
 public:
  // Used for creating subsets of the all-atoms group.
  AtomGroup(std::string name, std::vector<int> indices,
            const SystemGroup& all_atoms) : SubsystemGroup(name, indices) {
    positions_ = arma::zeros(group_size_, DIMS);
    velocities_ = arma::zeros(group_size_, DIMS);
    Init(name, all_atoms);
    copy_positions(all_atoms);
    copy_velocities(all_atoms);
  }

  // Creates .gro file
  void WriteGro(const std::string& gro, const arma::rowvec& box,
                const std::string description);

  // Deletes a molecule from the group. Sends references to the
  // center of mass position and velocity in case the user wishes
  // to replace it with another molecule.
  void RemoveMolecule(const int molecule_id, arma::rowvec& old_com_position,
                      arma::rowvec& old_com_velocity);

  // Adds a molecule to the group at the given position.
  // Currently only supports monatomic molecules.
  void AddMolecule(const int molecule_id, const Molecule& mol_to_add,
                   const arma::rowvec& position, const arma::rowvec& velocity);

  // Removes old molecule and adds new one at same position.
  // TODO(Zak) Separate into versions for all_atoms group and smaller groups.
  inline void ReplaceMolecule(const Molecule& new_molecule,
                                      const int mol_to_remove) {
    arma::rowvec position (DIMS);
    arma::rowvec velocity(DIMS);
    RemoveMolecule(mol_to_remove, position, velocity);
    // TODO(Zak): allow adjustement of velocity based
    // on Maxwell-Boltzmann distribution instead of
    // simple replacement.
    AddMolecule(group_size_, new_molecule, position, velocity);
  }

  // Just here to allow AtomGroup to mimic MoleculeGroup.
  inline void ZeroCom() { return; }

  // Just here to allow AtomGroup to mimic MoleculeGroup.
  inline void UpdateCom() { return; };

  // Mutators
  inline void set_ion_gammas(const double cation_gamma,
                             const double anion_gamma) {
    double tol = 1.0e-3;
    index_to_gamma_.resize(group_size_);
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      if (charge(i_atom) > tol)
        index_to_gamma_[i_atom] = cation_gamma;
      else if (charge(i_atom) < -tol)
        index_to_gamma_[i_atom] = anion_gamma;
    }
  }

  // Accessors
  inline int num_molecules() const { return group_size_; }
  inline arma::mat com_positions() const { return positions(); }
  inline arma::mat com_velocities() const { return velocities(); }

  inline arma::rowvec com_position(int molecule) const {
    return position(molecule);
  }

  inline arma::rowvec com_velocity(int molecule) const {
    return velocity(molecule);
  }

  inline double com_position(int molecule, int dim) const {
    return position(molecule, dim);
  }

  inline double com_velocity(int molecule, int dim) const {
    return velocity(molecule, dim);
  }

  inline arma::rowvec com_velocity_xy(int molecule) const {
    return velocity_xy(molecule);
  }

  inline int molecular_index(int atom) const { return atom; }

  // Used to be molecule_mass
  inline double molecular_mass(int molecule) const {
    return mass(molecule);
  }

  inline double molecular_charge(int molecule) const {
    return charge(molecule);
  }

  inline int molecular_index_to_molecule(int molecule) const {
    return index_to_molecule_[molecule];
  }

  inline void SetElectricField(const ParticleGroup& other_group,
                               const arma::rowvec& box) {
    ZeroElectricField();
    if (field_check_) {
      ReadElectricField(field_file_);
    } else {
      std::fill(included_in_field_.begin(), included_in_field_.end(), false);
      if (other_group.has_molecules()) {
        MarkNearbyAtoms(other_group, box, kFieldCutoffSquared(),
                        nearby_molecules);
        CalculateElectricField(other_group, box, nearby_molecules, dx);
      } else {
        CalculateElectricField(other_group, box, kFieldCutoffSquared(), dx);
      }
    }
  }

  inline void UpdateElectricField(const ParticleGroup& other_group,
                                  const arma::rowvec& box) {
    if (field_check_)
      return;
    if (other_group.has_molecules()) {
      MarkNearbyAtoms(other_group, box, kFieldCutoffSquared(),
                      nearby_molecules);
      CalculateElectricField(other_group, box, nearby_molecules, dx);
    } else {
      CalculateElectricField(other_group, box, kFieldCutoffSquared(), dx);
    }
  }

  inline void InducedDipole(const int atom, arma::rowvec& dipole,
                            const bool zero_first = false) const {
    if(zero_first)
      dipole = arma::zeros<arma::rowvec>(DIMS);
    dipole += electric_field(atom)*gamma(atom);
  }

 private:
  // Called by all constructors to create needed initialize
  // vectors and matrices.
  void Init(const std::string& name, const SystemGroup& all_atoms);

  // Takes care of the atom-level removal for RemoveMolecule
  void RemoveAtom(const int index);

  // Takes care of the atom-level addition for AddMolecule
  void AddAtom(const Atom atom_to_add, const int molecule_id,
                       const arma::rowvec& position,
                       const arma::rowvec& velocity);

 // Just here to allow AtomGroup to mimic MoleculeGroup.
  inline void CalculateMolecularCharges (
      const std::vector<double>& molecular_index_to_charge) {
    return;
  }

  inline void PermanentDipole(const int molecule, const arma::rowvec& box,
                              arma::rowvec& dipole,
                              const bool zero_first = false) {};

};

#endif
