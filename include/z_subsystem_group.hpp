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

#include "z_particle_group.hpp"
#include "z_system_group.hpp"

#ifndef _Z_SUBSYSTEM_GROUP_HPP_
#define _Z_SUBSYSTEM_GROUP_HPP_

//See description above
class SubsystemGroup:
  public ParticleGroup {
 public:
  SubsystemGroup(std::string name, std::vector<int> indices)
    : ParticleGroup(name, indices) {};


  static SubsystemGroup* MakeSubsystemGroup(std::string name,
                                            std::vector<int> indices,
                                            const SystemGroup& all_atoms);

  SubsystemGroup() {};

  virtual ~SubsystemGroup() {};

 protected:
  inline void copy_positions(const SystemGroup& all_atoms) {
    int i = 0;
    for (std::vector<int>::iterator i_index = indices_.begin();
         i_index != indices_.end(); ++i_index, i++) {
      positions_.row(i) = all_atoms.position(*i_index);
    }
  }

  inline void copy_velocities(const SystemGroup& all_atoms) {
    int i = 0;
    for (std::vector<int>::iterator i_index = indices_.begin();
         i_index != indices_.end(); ++i_index, ++i) {
      velocities_.row(i) = all_atoms.velocity(*i_index);
    }
  }

 private:
  // Called by all constructors to create needed initialize
  // vectors and matrices.
  virtual void Init(const std::string& name, const SystemGroup& all_atoms) = 0;
};

#endif
