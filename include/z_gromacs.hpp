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

// A set of i/o functions for use with file formats associated
// with the GROMACS molecular dynamics simulation package.

#include <map>
#include <string>
#include <armadillo>
#include "xdrfile_xtc.h"
#include "z_sim_params.hpp"

#ifndef _Z_GROMACS_HPP_
#define _Z_GROMACS_HPP_

// Determines types of atoms present in topology file
extern void ReadAtomTypes(std::ifstream& top_file,
                            std::vector<Atom>& atom_types);

// Makes vector of molecule types present in topology file
extern std::vector<Molecule> GenMolecules(const std::string& top,
                                           SimParams& params);

// Reads atom-index/group mapping from index file.
extern std::map<std::string, std::vector<int> > ReadNdx(
    const std::string& ndx);

// Extracts single group from group map.
extern std::vector<int> SelectGroup(
    std::map<std::string, std::vector<int> >& groups,
    const std::string& group_name);

// Allows conversion of rvec to arma::rowvec for compatibility
// between xtc library and armadillo library.
extern arma::rowvec RvecToRow(const rvec& vec1);

#endif
