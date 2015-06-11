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

// Implementation of SimParams class. See include/z_sim_params.hpp for
// more details about the class.

#include "z_sim_params.hpp"
#include <vector>
#include <cassert>
#include <boost/lexical_cast.hpp>
#include "z_string.hpp"
#include "z_conversions.hpp"

//TODO(Zak): replace this with a boost/program_options call
void SimParams::ReadParams() {
  std::string line;
  std::string subline;
  file_.open(filename_.c_str());
  while (getline(file_, line)) {
    if(line[0] == '#' || line[0] == '@') continue;
    std::vector<std::string> split_line = Split(line, '=');
    std::string variable = Trim(split_line[0]);
    std::string value = Trim(split_line[1]);
    assert(!variable.empty() && !value.empty());
    if (variable == "numMols")
      numMols_ = boost::lexical_cast<int>(value.c_str());
    else if (variable == "numCations")
      numCations_ = boost::lexical_cast<int>(value.c_str());
    else if (variable == "numAnions")
      numAnions_ = boost::lexical_cast<int>(value.c_str());
    else if (variable == "numOthers")
      numOthers_ = boost::lexical_cast<int>(value.c_str());
    else if (variable == "numSteps")
      max_steps_ = boost::lexical_cast<int>(value.c_str());
    else if (variable == "rminOO")
      rminOO_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rminC")
      rminC_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rminA")
      rminA_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rminAH")
      rminAH_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rmin2C")
      rmin2C_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rmin2A")
      rmin2A_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rmin3C")
      rmin3C_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "rmin3A")
      rmin3A_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "dt")
      dt_ = boost::lexical_cast<double>(value.c_str());
    else if (variable == "temp")
      set_temperature(boost::lexical_cast<double>(value.c_str()));
    else if (variable == "gamm")
      gamm_ =
          boost::lexical_cast<double>(value.c_str())*pow(ANG_TO_NM,3)*ENM_TO_D;
    else if (variable == "gammC")
      gammC_ =
          boost::lexical_cast<double>(value.c_str())*pow(ANG_TO_NM,3)*ENM_TO_D;
    else if (variable == "gammA")
      gammA_ =
          boost::lexical_cast<double>(value.c_str())*pow(ANG_TO_NM,3)*ENM_TO_D;
    else
      assert(false && "Param option not recognized");
  }
  numChromos_ = 2.0*numMols_;
  numIons_ = numCations_+numAnions_;
  rminOO2_ = rminOO_*rminOO_;
  rminC2_ = rminC_*rminC_;
  rminA2_ = rminA_*rminA_;
  rminAH2_ = rminAH_*rminAH_;
  rmin2C2_ = rmin2C_*rmin2C_;
  rmin2A2_ = rmin2A_*rmin2A_;
  rmin3C2_ = rmin3C_*rmin3C_;
  rmin3A2_ = rmin3A_*rmin3A_;
  file_.close();
}

void SimParams::ExtractTrajMetadata(char *traj, rvec **x_in,
                                       arma::rowvec& box) {
  int st;
  matrix box_mat;
  float t1, t2, prec;
  XDRFILE *traj_file;
  read_xtc_natoms(traj, &num_atoms_);
  traj_file = xdrfile_open(traj, "r");
  *x_in = new rvec [num_atoms_];
  read_xtc(traj_file, num_atoms_, &st, &t1, box_mat, *x_in, &prec);
  read_xtc(traj_file, num_atoms_, &st, &t2, box_mat, *x_in, &prec);
  dt_ = t2 - t1;
  for(int i=0; i<DIMS; i++)
    box(i) = box_mat[i][i];
  xdrfile_close(traj_file);
}
