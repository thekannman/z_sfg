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
#ifndef _Z_MOLECULE_GROUP_HPP_
#define _Z_MOLECULE_GROUP_HPP_

//See description above
class MoleculeGroup:
  public SubsystemGroup {
 public:
  // Used for creating subsets of the all-atoms group.
  MoleculeGroup(std::string name, std::vector<int> indices,
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
    AddMolecule(num_molecules_, new_molecule, position, velocity);
  }

  // Resets center of mass position and velocity of all molecules.
  inline void ZeroCom() {
    com_positions_ = arma::zeros(num_molecules_, DIMS);
    com_velocities_ = arma::zeros(num_molecules_, DIMS);
  }

  // Sets center of mass position and velocity of all molecules.
  void UpdateCom();

  // Mutators
  inline void set_ion_gammas(const double cation_gamma,
                             const double anion_gamma) {
    double tol = 1.0e-3;
    index_to_gamma_.resize(group_size_);
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      if (index_to_molecular_charge_[i_atom] > tol)
        index_to_gamma_[i_atom] = cation_gamma;
      else if (index_to_molecular_charge_[i_atom] < -tol)
        index_to_gamma_[i_atom] = anion_gamma;
    }
  }

  // Accessors
  inline int first_index_in_molecule(const int molecule) const {
    return first_index_in_molecule_[molecule];
  }
  inline int last_index_in_molecule(const int molecule) const {
    return last_index_in_molecule_[molecule];
  }
  inline int num_molecules() const { return num_molecules_; }
  inline arma::mat com_positions() const { return com_positions_; }
  inline arma::mat com_velocities() const { return com_velocities_; }

  inline arma::rowvec com_position(int molecule) const {
    return com_positions_.row(molecule);
  }

  inline arma::rowvec com_velocity(int molecule) const {
    return com_velocities_.row(molecule);
  }

  inline double com_position(int molecule, int dim) const {
    return com_positions_(molecule, dim);
  }

  inline double com_velocity(int molecule, int dim) const {
    return com_velocities_(molecule, dim);
  }

  inline arma::rowvec com_velocity_xy(int molecule) const {
    return com_velocities_(molecule, arma::span(0,1));
  }

  inline int molecular_index(int atom) const {
    return index_to_molecular_index_[atom];
  }

  // Used to be molecule_mass
  inline double molecular_mass(int molecule) const {
    return molecular_index_to_mass_[molecule];
  }

  inline double molecular_charge(int molecule) const {
    return molecular_index_to_charge_[molecule];
  }

  inline int molecular_index_to_molecule(int molecule) const {
    return molecular_index_to_molecule_[molecule];
  }

  inline void PermanentDipole(const int molecule, const arma::rowvec& box,
                              arma::rowvec& dipole,
                              const bool zero_first = false) {
    if(zero_first)
      dipole = arma::zeros<arma::rowvec>(DIMS);
    int first_index = first_index_in_molecule_[molecule];
    int last_index = last_index_in_molecule_[molecule];
    for (int i_atom = first_index+1; i_atom <= last_index; ++i_atom) {
      FindDxNoShift(dx, position(first_index), position(i_atom), box);
      dipole += dx*index_to_charge_[i_atom];
    }
  }

  inline void InducedDipole(const int molecule, arma::rowvec& dipole,
                            const bool zero_first = false) const {
    if(zero_first)
      dipole = arma::zeros<arma::rowvec>(DIMS);
    int first_index = first_index_in_molecule_[molecule];
    int last_index = last_index_in_molecule_[molecule];
    for (int i_atom = first_index+1; i_atom <= last_index; ++i_atom) {
      dipole += electric_field(i_atom)*gamma(i_atom);
    }
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

 protected:

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

  inline void CalculateMolecularCharges (
      const std::vector<double>& molecular_index_to_charge) {
    index_to_molecular_charge_.resize(group_size_);
    for (int i_mol = 0; i_mol < num_molecules_; ++i_mol) {
      molecular_index_to_charge_[i_mol] =
          molecular_index_to_charge[molecular_index_to_molecule(i_mol)];
    }
    for (int i_atom = 0; i_atom < group_size_; ++i_atom) {
      index_to_molecular_charge_[i_atom] =
          molecular_index_to_charge[index_to_molecular_index_[i_atom]];
    }
  }

  std::vector<double> index_to_molecular_charge_;
  std::vector<int> index_to_molecular_index_;
  std::vector<double> molecular_index_to_mass_;
  std::vector<double> molecular_index_to_charge_;
  std::vector<int>molecular_index_to_molecule_;
  int num_molecules_;
  arma::mat com_positions_;
  arma::mat com_velocities_;
  std::vector<int> first_index_in_molecule_;
  std::vector<int> last_index_in_molecule_;
};

#endif
