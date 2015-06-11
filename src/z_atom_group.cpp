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

// Implementation of AtomGroup class. See include/z_atom_group.hpp for
// more details about the class.

#include "z_atom_group.hpp"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>
#include "z_string.hpp"
#include "z_vec.hpp"
#include <algorithm>

void AtomGroup::WriteGro(const std::string& gro, const arma::rowvec& box,
                         const std::string description) {
    std::ofstream gro_file;
    gro_file.open(gro.c_str());
    assert(gro_file.is_open());
    gro_file << description << std::endl;
    gro_file << group_size_ << std::endl;
    for (int i_atom = 0; i_atom < group_size_; i_atom++) {
      gro_file << boost::format("%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f") %
          (i_atom+1) %
          molecule_name(i_atom) %
          atom_name(i_atom) % (indices_[i_atom]+1) % positions_(i_atom,0) %
          positions_(i_atom,1) % positions_(i_atom,2) %
          velocities_(i_atom,0) % velocities_(i_atom,1) %
          velocities_(i_atom,2) << std::endl;
    }
    gro_file << boost::format("%10.5f%10.5f%10.5f") % box(0) % box(1) % box(2);
    gro_file.close();
}

void AtomGroup::RemoveMolecule(const int molecule_id,
                               arma::rowvec& old_com_position,
                               arma::rowvec& old_com_velocity) {
  const int marked_molecule = index_to_molecule(molecule_id);
  for (std::vector<int>::iterator i_mol = index_to_molecule_.begin() ;
       i_mol != index_to_molecule_.end(); i_mol++) {
    if (*i_mol > marked_molecule)
      (*i_mol) -= 1;
  }
  for (std::vector<int>::iterator i_atom = indices_.begin();
       i_atom != indices_.end(); ++i_atom) {
    // Decrements molecular and atomic indices to eliminate gap
    if ((*i_atom) > molecule_id) {
      (*i_atom) -= 1;
    }
  }
  molecule_names_.erase(molecule_names_.begin() + molecule_id);
  // Save references to position/velocity for later replacement.
  old_com_position = position(molecule_id);
  old_com_velocity = velocity(molecule_id);
  RemoveAtom(molecule_id);
}


void AtomGroup::AddMolecule(const int molecule_id, const Molecule& molecule,
                            const arma::rowvec& position,
                            const arma::rowvec& velocity) {
  assert (molecule.num_atoms() == 1);
  molecule_names_.push_back(molecule.name());
  AddAtom(molecule.front(), molecule_id, position, velocity);
}

void AtomGroup::RemoveAtom(const int index) {
  indices_.erase(indices_.begin() + index);
  atom_names_.erase(atom_names_.begin() + index);
  index_to_mass_.erase(index_to_mass_.begin() + index);
  index_to_molecule_.erase(index_to_molecule_.begin() + index);
  index_to_sigma_.erase(index_to_sigma_.begin() + index);
  index_to_epsilon_.erase(index_to_epsilon_.begin() + index);
  positions_.shed_row(index);
  velocities_.shed_row(index);
  group_size_--;
}

void AtomGroup::AddAtom(const Atom atom_to_add, const int molecule_id,
                        const arma::rowvec& position,
                        const arma::rowvec& velocity) {
  atom_names_.push_back(atom_to_add.name());
  index_to_mass_.push_back(atom_to_add.mass());
  index_to_molecule_.push_back(molecule_id);
  index_to_sigma_.push_back(atom_to_add.sigma());
  index_to_epsilon_.push_back(atom_to_add.epsilon());
  // TODO(Zak): Fix to work outside of all_atoms group.
  indices_.push_back(group_size_);
  group_size_++;
  positions_.insert_rows(positions_.n_rows, position);
  velocities_.insert_rows(velocities_.n_rows, velocity);
}

void AtomGroup::Init(const std::string& name, const SystemGroup& all_atoms) {
  for (std::vector<int>::iterator i_index = indices_.begin();
       i_index != indices_.end(); ++i_index) {
    index_to_mass_.push_back(all_atoms.mass(*i_index));
    index_to_molecule_.push_back(all_atoms.index_to_molecule(*i_index));
    index_to_sigma_.push_back(all_atoms.sigma(*i_index));
    index_to_epsilon_.push_back(all_atoms.epsilon(*i_index));
    index_to_charge_.push_back(all_atoms.charge(*i_index));
    atom_names_.push_back(all_atoms.atom_name(*i_index));
    molecule_names_.push_back(all_atoms.molecule_name(*i_index));
  }
  set_has_molecules(false);
  included_in_field_.resize(all_atoms.size());
}
