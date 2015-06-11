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

#include <cassert>
#include <armadillo>
#include "z_constants.hpp"
#include "z_molecule.hpp"
#include "z_gromacs.hpp"
#include "z_string.hpp"
#include "z_histogram.hpp"
#include "z_vec.hpp"
#ifndef _Z_PARTICLE_GROUP_HPP_
#define _Z_PARTICLE_GROUP_HPP_

//See description above
class ParticleGroup {
 public:
  ParticleGroup(std::string name, std::vector<int> indices)
      : indices_(indices), group_size_(indices.size()),
        field_filename_(name+"_electric_field.dat"), name_(name) {};

  ParticleGroup() {};

  // Creates .gro file
  virtual void WriteGro(const std::string& gro, const arma::rowvec& box,
                        const std::string description) = 0;

  // Deletes a molecule from the group. Sends references to the
  // center of mass position and velocity in case the user wishes
  // to replace it with another molecule.
  virtual void RemoveMolecule(const int molecule_id,
                              arma::rowvec& old_com_position,
                              arma::rowvec& old_com_velocity) = 0;

  // Wrapper for the above function which throws away references
  // to position/velocity vectors. Used if replacement is not planned.
  inline void RemoveMolecule(const int molecule_id) {
    arma::rowvec position(DIMS);
    arma::rowvec velocity(DIMS);
    RemoveMolecule(molecule_id, position, velocity);
  }

  // Adds a molecule to the group at the given position.
  // Currently only supports monatomic molecules.
  virtual void AddMolecule(const int molecule_id, const Molecule& mol_to_add,
                           const arma::rowvec& position,
                           const arma::rowvec& velocity) = 0;

  // Removes old molecule and adds new one at same position.
  // TODO(Zak) Separate into versions for all_atoms group and smaller groups.
  virtual inline void ReplaceMolecule(const Molecule& new_molecule,
                                      const int mol_to_remove) = 0;

  // Resets center of mass position and velocity of all molecules.
  virtual inline void ZeroCom() = 0;

  // Sets center of mass position and velocity of all molecules.
  virtual void UpdateCom() = 0;

  inline void ZeroElectricField() {
    electric_fields_ = arma::zeros(group_size_, DIMS);
  }

  void CalculateElectricField(const ParticleGroup& other_group,
                              const arma::rowvec& box,
                              const arma::imat& nearby_molecules,
                              arma::rowvec& dx);

  void CalculateElectricField(const ParticleGroup& other_group,
                              const arma::rowvec& box,
                              const double cutoff_squared,
                              arma::rowvec& dx);

  void MarkNearbyAtoms(const ParticleGroup& other_group,
                              const arma::rowvec& box,
                              const double cutoff_squared,
                              arma::imat& nearby_atoms) const;

  // Mutators
  inline void set_has_molecules(const bool has_mol) {
    has_molecules_ = has_mol;
  }
  inline void set_name(const std::string& name) {
    name_ = name;
  }

  inline void set_position(const int i, const rvec& position) {
    positions_.row(i) = RvecToRow(position);
  }

  inline void set_velocity(const int i, const rvec& velocity) {
    velocities_.row(i) = RvecToRow(velocity);
  }

  inline void set_positions(const rvec x_in[]) {
    int i = 0;
    for (std::vector<int>::iterator i_atom = indices_.begin();
         i_atom != indices_.begin(); ++i_atom, ++i) {
      set_position(i, x_in[*i_atom]);
    }
  }

  inline void set_velocities(const rvec v_in[]) {
    int i = 0;
    for (std::vector<int>::iterator i_atom = indices_.begin();
         i_atom != indices_.begin(); ++i_atom, ++i) {
      set_velocity(i, v_in[*i_atom]);
    }
  }

  inline void set_phase_space(const rvec x_in[], const rvec v_in[]) {
    int i = 0;
    for (std::vector<int>::iterator i_atom = indices_.begin();
         i_atom != indices_.begin(); ++i_atom, ++i) {
      set_position(i, x_in[*i_atom]);
      set_velocity(i, v_in[*i_atom]);
    }
  }

  inline void set_gammas(const double gamma) {
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      index_to_gamma_.push_back(gamma);
    }
  }

  virtual inline void set_ion_gammas(const double cation_gamma,
                                     const double anion_gamma) = 0;

  // Accessors
  inline double kFieldCutoffSquared() const { return kFieldCutoffSquared_; }
  inline bool has_molecules() const { return has_molecules_; }
  inline bool field_check() const { return field_check_; }
  inline std::string name() const { return name_; }
  inline std::vector<int>::const_iterator begin() const {
    return indices_.begin();
  }
  inline std::vector<int>::const_iterator end() const { return indices_.end(); }
  inline int size() const { return group_size_; }
  inline std::vector<int> indices() const { return indices_; }
  inline int indices(int i) const { return indices_[i]; }
  virtual int num_molecules() const = 0;
  inline arma::mat positions() const { return positions_; }
  inline arma::mat velocities() const { return velocities_; }
  virtual inline arma::mat com_positions() const = 0;
  virtual inline arma::mat com_velocities() const = 0;
  inline arma::rowvec position(int atom) const { return positions_.row(atom); }
  inline arma::rowvec velocity(int atom) const { return velocities_.row(atom); }
  inline arma::rowvec electric_field(int atom) const {
    return electric_fields_.row(atom);
  }
  inline double mass(int atom) const { return index_to_mass_[atom]; }
  inline double sigma(int atom) const { return index_to_sigma_[atom]; }
  inline double epsilon(int atom) const { return index_to_epsilon_[atom]; }
  inline double charge(int atom) const { return index_to_charge_[atom]; }
  inline double gamma(int atom) const { return index_to_gamma_[atom]; }

  virtual inline arma::rowvec com_position(int molecule) const = 0;

  virtual inline arma::rowvec com_velocity(int molecule) const = 0;

  inline double position(int atom, int dim) const {
    return positions_(atom, dim);
  }

  inline double velocity(int atom, int dim) const {
    return velocities_(atom, dim);
  }

  inline arma::rowvec velocity_xy(int atom) const {
    return velocities_(atom, arma::span(0,1));
  }

  virtual inline double com_position(int molecule, int dim) const = 0;

  virtual inline double com_velocity(int molecule, int dim) const = 0;

  virtual inline arma::rowvec com_velocity_xy(int molecule) const = 0;

  virtual inline int molecular_index(int atom) const = 0;

  // Used to be molecule_mass
  virtual inline double molecular_mass(int molecule) const = 0;

  virtual inline double molecular_charge(int molecule) const = 0;

  virtual inline int molecular_index_to_molecule(int molecule) const = 0;

  inline int index_to_molecule(int index) const {
    return index_to_molecule_[index];
  }

  inline std::string atom_name(int index) const {
    return atom_names_[index];
  }

  inline std::string molecule_name(int index) const {
    return molecule_names_[index];
  }

  inline void WriteElectricField() {
    assert(field_file_.is_open());
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      for (int i_dim = 0; i_dim < DIMS; i_dim++)
        field_file_ << std::scientific << electric_fields_(i_atom,i_dim) << " ";
    }
    field_file_ << std::endl;
  }

  inline void ReadElectricField(std::fstream& field_file) {
    std::string line;
    assert(field_file.is_open());
    getline(field_file, line);
    const std::vector<std::string> split_line = Split(line, ' ');
    std::vector<std::string>::const_iterator i_split = split_line.begin();
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      for (int i_dim = 0; i_dim < DIMS; i_dim++) {
        electric_fields_(i_atom,i_dim) = atof((*i_split).c_str());
        i_split++;
      }
    }
  }

  void FindClusters(const double cutoff_distance_squared,
                    const arma::rowvec& box, Histogram& clusters) const;

  // Finds the k atoms that are nearest to some point. The atoms are returned
  // in distance order.
  void FindNearestk(const arma::rowvec& point, const arma::rowvec& box,
                    const int exclude_molecule, const int k, arma::rowvec& dx,
                    arma::rowvec& r2, arma::irowvec& nearest) const;

  inline void OpenFieldFile() {
    field_file_.open(field_filename_.c_str(), std::fstream::in);
    if (!field_file_.is_open()) {
      field_file_.open(field_filename_.c_str(), std::fstream::out);
      field_check_ = false;
    } else {
      std::cout << "Using previous electric field data from file: " <<
                 field_filename_ << "." << std::endl;
      field_check_ = true;
    }
  }

  virtual inline void PermanentDipole(const int molecule,
                                      const arma::rowvec& box,
                                      arma::rowvec& dipole,
                                      const bool zero_first = false) = 0;

  virtual inline void InducedDipole(const int molecule, arma::rowvec& dipole,
                                    const bool zero_first = false) const = 0;

  virtual inline void SetElectricField(const ParticleGroup& other_group,
                                       const arma::rowvec& box) = 0;

  virtual inline void UpdateElectricField(const ParticleGroup& other_group,
                                          const arma::rowvec& box) = 0;

 protected:
  std::vector<double> index_to_gamma_;
  std::vector<int> index_to_molecule_;
  std::vector<std::string> molecule_names_;
  arma::mat positions_;
  arma::mat velocities_;
  arma::mat electric_fields_;
  std::vector<int> indices_;
  std::vector<std::string> atom_names_;
  std::vector<double> index_to_mass_;
  std::vector<double> index_to_sigma_;
  std::vector<double> index_to_epsilon_;
  std::vector<double> index_to_charge_;
  int group_size_;
  bool has_molecules_;
  arma::imat nearby_molecules;
  arma::rowvec dx;
  bool field_check_;
  std::vector<bool> included_in_field_;
  std::fstream field_file_;
  std::string field_filename_;

  // Indices of atoms for use in finding whether molecule is nearby.
  // Primarily for electric field calculations.
  std::vector<int> molecule_position_markers_;

 private:
  static const double kFieldCutoff_ = 0.7831;
  static const double kFieldCutoffSquared_ = 0.61324561;
  std::string name_;

  // Takes care of the atom-level removal for RemoveMolecule
  virtual void RemoveAtom(const int index) = 0;

  // Takes care of the atom-level addition for AddMolecule
  virtual void AddAtom(const Atom atom_to_add, const int molecule_id,
                       const arma::rowvec& position,
                       const arma::rowvec& velocity) = 0;

  void FindCluster(const double cutoff_distance_squared,
                   const arma::rowvec& box, arma::rowvec& dx,
                   arma::irowvec& in_current_cluster,
                   arma::irowvec& in_any_cluster,
                   const int start_atom) const;

  virtual inline void CalculateMolecularCharges (
      const std::vector<double>& molecular_index_to_charge) = 0;

};

#endif
