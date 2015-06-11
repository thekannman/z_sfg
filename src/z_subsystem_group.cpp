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

#include "z_subsystem_group.hpp"
#include "z_molecule_group.hpp"
#include "z_atom_group.hpp"
#include "z_string.hpp"
#include "z_vec.hpp"
#include <algorithm>

SubsystemGroup* SubsystemGroup::MakeSubsystemGroup(
    std::string name, std::vector<int> indices, const SystemGroup& all_atoms) {
  std::vector<int> atom_count(all_atoms.num_molecules());
  for (std::vector<int>::const_iterator i_index = indices.begin();
       i_index != indices.end(); ++i_index) {
    if (atom_count[all_atoms.index_to_molecule(*i_index)] > 0)
      return new MoleculeGroup(name, indices, all_atoms);
    atom_count[all_atoms.index_to_molecule(*i_index)]++;
  }
  return new AtomGroup(name, indices, all_atoms);
}

