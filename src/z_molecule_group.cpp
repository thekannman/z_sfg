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

// Implementation of MoleculeGroup class. See include/z_atom_group.hpp for
// more details about the class.

#include "z_molecule_group.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "z_string.hpp"
#include "z_vec.hpp"
#include <algorithm>

void MoleculeGroup::WriteGro(const std::string& gro, const arma::rowvec& box,
                             const std::string description) {
    std::ofstream gro_file;
    gro_file.open(gro.c_str());
    assert(gro_file.is_open());
    gro_file << description << std::endl;
    gro_file << group_size_ << std::endl;
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      gro_file << boost::format("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f") %
          (index_to_molecular_index_[i_atom]+1) %
          molecule_name(index_to_molecular_index_[i_atom]) %
          atom_name(i_atom) % (indices_[i_atom]+1) % positions_(i_atom,0) %
          positions_(i_atom,1) % positions_(i_atom,2) %
          velocities_(i_atom,0) % velocities_(i_atom,1) %
          velocities_(i_atom,2) << std::endl;
    }
    gro_file << boost::format("%10.5f%10.5f%10.5f") % box(0) % box(1) % box(2);
    gro_file.close();
}

void MoleculeGroup::RemoveMolecule(const int molecule_id,
                               arma::rowvec& old_com_position,
                               arma::rowvec& old_com_velocity) {
  std::vector<int> mark_for_rem;
  const int marked_molecule = molecular_index_to_molecule_[molecule_id];
  for (std::vector<int>::iterator i_mol = index_to_molecule_.begin();
       i_mol != index_to_molecule_.end(); i_mol++) {
    if (*i_mol == marked_molecule)
      mark_for_rem.push_back(std::distance(index_to_molecule_.begin(),i_mol));
    else if (*i_mol > marked_molecule)
      (*i_mol) -= 1;
  }
  for (std::vector<int>::iterator i_mol = index_to_molecular_index_.begin(),
           i_atom = indices_.begin();
       i_mol != index_to_molecular_index_.end(); ++i_mol, ++i_atom) {
    // Decrements molecular and atomic indices to eliminate gap
    if ((*i_mol) > molecule_id) {
      (*i_mol) -= 1;
      (*i_atom) -= mark_for_rem.size();
    }
  }
  assert(!mark_for_rem.empty());
  for (std::vector<int>::iterator i_mol =
           molecular_index_to_molecule_.begin()+molecule_id+1;
       i_mol != molecular_index_to_molecule_.end(); ++i_mol) {
    (*i_mol) -= 1;
  }
  molecular_index_to_mass_.erase(
      molecular_index_to_mass_.begin() + molecule_id);
  molecular_index_to_molecule_.erase(
      molecular_index_to_molecule_.begin() + molecule_id);
  molecule_names_.erase(
      molecule_names_.begin() + molecule_id);
  // Save references to position/velocity for later replacement.
  old_com_position = com_positions_.row(molecule_id);
  old_com_velocity = com_velocities_.row(molecule_id);
  com_positions_.shed_row(molecule_id);
  com_velocities_.shed_row(molecule_id);
  for (std::vector<int>::reverse_iterator i_atom = mark_for_rem.rbegin();
       i_atom != mark_for_rem.rend(); i_atom++) {
    RemoveAtom(*i_atom);
  }
  num_molecules_--;
}


void MoleculeGroup::AddMolecule(const int molecule_id, const Molecule& molecule,
                            const arma::rowvec& position,
                            const arma::rowvec& velocity)
{
    // TODO(Zak): Allow for polyatomic molecules by randomly
    // orienting around center of mass position.
    // This assert is a placeholder until then.
    assert (molecule.num_atoms() == 1);
    molecule_names_.push_back(molecule.name());
    molecular_index_to_mass_.push_back(molecule.mass());
    molecular_index_to_molecule_.push_back(molecule_id);
    for (std::vector<Atom>::const_iterator i_atom = molecule.begin();
         i_atom != molecule.end(); ++i_atom) {
      AddAtom(*i_atom, molecule_id, position, velocity);
    }
    num_molecules_++;
    com_positions_.insert_rows(com_positions_.n_rows, position);
    com_velocities_.insert_rows(com_velocities_.n_rows, velocity);
}

void MoleculeGroup::UpdateCom() {
  ZeroCom();
  for (int i_atom = 0; i_atom < group_size_; i_atom++) {
    com_positions_.row(index_to_molecular_index_[i_atom]) +=
        index_to_mass_[i_atom]*positions_.row(i_atom);
    com_velocities_.row(index_to_molecular_index_[i_atom]) +=
        index_to_mass_[i_atom]*velocities_.row(i_atom);
  }
  for (int i_mol = 0; i_mol < num_molecules_; ++i_mol) {
    com_positions_.row(i_mol) /= molecular_index_to_mass_[i_mol];
    com_velocities_.row(i_mol) /= molecular_index_to_mass_[i_mol];
  }
}

void MoleculeGroup::RemoveAtom(const int index) {
  indices_.erase(indices_.begin() + index);
  atom_names_.erase(atom_names_.begin() + index);
  index_to_mass_.erase(index_to_mass_.begin() + index);
  index_to_molecule_.erase(index_to_molecule_.begin() + index);
  index_to_sigma_.erase(index_to_sigma_.begin() + index);
  index_to_epsilon_.erase(index_to_epsilon_.begin() + index);
  index_to_molecular_index_.erase(
      index_to_molecular_index_.begin() + index);
  positions_.shed_row(index);
  velocities_.shed_row(index);
  group_size_--;
}

void MoleculeGroup::AddAtom(const Atom atom_to_add, const int molecule_id,
                        const arma::rowvec& position,
                        const arma::rowvec& velocity) {
  atom_names_.push_back(atom_to_add.name());
  index_to_mass_.push_back(atom_to_add.mass());
  index_to_molecule_.push_back(molecule_id);
  index_to_sigma_.push_back(atom_to_add.sigma());
  index_to_epsilon_.push_back(atom_to_add.epsilon());
  index_to_molecular_index_.push_back(num_molecules_);
  // TODO(Zak): Fix to work outside of all_atoms group.
  indices_.push_back(group_size_);
  group_size_++;
  positions_.insert_rows(positions_.n_rows, position);
  velocities_.insert_rows(velocities_.n_rows, velocity);
}

void MoleculeGroup::Init(const std::string& name,
                         const SystemGroup& all_atoms) {
  last_index_in_molecule_.resize(group_size_);
  for (std::vector<int>::iterator i_index = indices_.begin();
       i_index != indices_.end(); ++i_index) {
    index_to_mass_.push_back(all_atoms.mass(*i_index));
    index_to_molecule_.push_back(all_atoms.index_to_molecule(*i_index));
    index_to_sigma_.push_back(all_atoms.sigma(*i_index));
    index_to_epsilon_.push_back(all_atoms.epsilon(*i_index));
    index_to_charge_.push_back(all_atoms.charge(*i_index));
    atom_names_.push_back(all_atoms.atom_name(*i_index));
  }
  std::vector<double>::iterator i_mass = index_to_mass_.begin();
  for (std::vector<int>::iterator i_mol = index_to_molecule_.begin();
       i_mol != index_to_molecule_.end(); ++i_mol, ++i_mass) {
    std::vector<int>::iterator mol_index;
    mol_index = std::find(molecular_index_to_molecule_.begin(),
                          molecular_index_to_molecule_.end(), *i_mol);
    if (mol_index == molecular_index_to_molecule_.end()) {
      molecular_index_to_mass_.push_back(0.0);
      first_index_in_molecule_.push_back(index_to_molecular_index_.size());
      index_to_molecular_index_.push_back(molecular_index_to_molecule_.size());
      molecular_index_to_molecule_.push_back(*i_mol);
      molecule_names_.push_back(molecule_name(*i_mol));
      // TODO(Zak): Find better way to determine identity of marker atoms.
      //            For now, just use first atom as this works for water.
      molecule_position_markers_.push_back(
          std::distance(index_to_molecule_.begin(), i_mol));
    } else {
      double mol_number = std::distance(molecular_index_to_molecule_.begin(),
                                        mol_index);
      last_index_in_molecule_[mol_number] = (index_to_molecular_index_.size());
      index_to_molecular_index_.push_back(mol_number);
      }
    molecular_index_to_mass_.back() += *i_mass;
  }
  num_molecules_ = molecular_index_to_molecule_.size();
  com_positions_.set_size(num_molecules_, DIMS);
  com_velocities_.set_size(num_molecules_, DIMS);

  CalculateMolecularCharges(all_atoms.molecular_charges());
  set_has_molecules(true);
}
