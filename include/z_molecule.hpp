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

// Contains information about a type of molecule as read from .top file.
// This is necessary for creation of atom groups from .gro file.

#include <vector>
#include <fstream>
#include "z_atom.hpp"

#ifndef _Z_MOLECULE_HPP_
#define _Z_MOLECULE_HPP_

class Molecule {
 public:
  Molecule() : num_atoms_(0), mass_(0.0) {}

  void AddAtom(Atom atom);

  // Better than overloaded << operator for printing to std::cout
  void Print() const;

  //Mutators
  inline void set_name(const std::string& name) { name_ = name; }

  //Accessors
  inline int  num_atoms() const { return num_atoms_; }
  inline std::vector<Atom> atoms() const { return atoms_; }
  inline std::vector<Atom>::const_iterator begin() const {
    return atoms_.begin(); }
  inline const Atom front() const {
    return atoms_.front(); }
  inline std::vector<Atom>::const_iterator end() const { return atoms_.end(); }

  inline std::string name() const { return name_; }
  inline double mass() const { return mass_; }

 private:
  int num_atoms_;
  std::vector<Atom> atoms_;
  std::string name_;
  double mass_;

  friend std::ostream &operator<<( std::ostream &output, Molecule &molecule );
};

#endif
