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

#include <complex>
#include "xdrfile_trr.h"
#include "boost/program_options.hpp"
#include "z_sim_params.hpp"
#include "z_vec.hpp"
#include "z_cx_tcf.hpp"
#include "z_constants.hpp"
#include "z_conversions.hpp"
#include "z_molecule.hpp"
#include "z_subsystem_group.hpp"
#include "z_gromacs.hpp"
#include "z_sfg_map.hpp"

namespace po = boost::program_options;
// Units are nm, ps.

int main (int argc, char *argv[]) {
  SimParams params;

  po::options_description desc("Options");
  desc.add_options()
    ("help,h",  "Print help messages")
    ("oxygen,O", po::value<std::string>()->default_value("OW"),
     "Group consisting of water oxygens for distance calculations")
    ("water,W", po::value<std::string>()->default_value("SOL"),
     "Group consisting of water molecules for electric field and terahertz"
     "calculations")
    ("solute,S", po::value<std::string>()->default_value("Ion"),
     "Group consisting of charged solutes for electric field calculations")
    ("index,n", po::value<std::string>()->default_value("index.ndx"),
     ".ndx file containing atomic indices for groups")
    ("gro", po::value<std::string>()->default_value("conf.gro"),
     ".gro file containing list of atoms/molecules")
    ("top", po::value<std::string>()->default_value("topol.top"),
     ".top file containing atomic/molecular properties")
    ("max_time,t", po::value<double>()->default_value(0.0),
     "Maximum simulation time to use in calculations")
    ("chromophore,c", po::value<std::string>()->required(),
     "Identity of chromophore (OH or OD)");

  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);
  if (vm.count("help")) {
    std::cout << desc << "\n";
    exit(EXIT_SUCCESS);
  }

  Chromophore chromophore;
  if (vm["chromophore"].as<std::string>() == "OH")
    chromophore = kOH;
  else if (vm["chromophore"].as<std::string>() == "OD")
    chromophore = kOD;
  else
    assert(false && "Chromophore identity unknown.");

  std::map<std::string, std::vector<int> > groups;
  groups = ReadNdx(vm["index"].as<std::string>());
  std::vector<Molecule> molecules = GenMolecules(vm["top"].as<std::string>(),
                                                 params);

  SystemGroup all_atoms(vm["gro"].as<std::string>(), molecules);

  MoleculeGroup hydrogen_group(vm["hydrogen"].as<std::string>(),
                               SelectGroup(groups,
                                           vm["hydrogen"].as<std::string>()),
                               all_atoms);
  AtomGroup oxygen_group(vm["oxygen"].as<std::string>(),
                         SelectGroup(groups, vm["oxygen"].as<std::string>()),
                         all_atoms);
  MoleculeGroup water_group(vm["water"].as<std::string>(),
                            SelectGroup(groups, vm["water"].as<std::string>()),
                            all_atoms);

  bool just_water = true;
  std::vector<int> solute_indices;
  std::string solute_group_name;
  if (vm.count("solute")) {
    solute_indices = SelectGroup(groups, vm["solute"].as<std::string>());
    just_water = false;
  }
  solute_group_name = vm.count("solute") ? vm["solute"].as<std::string>() : "";

  SubsystemGroup *solute_group_pointer =
      SubsystemGroup::MakeSubsystemGroup(solute_group_name, solute_indices,
                                         all_atoms);
  SubsystemGroup &solute_group = *solute_group_pointer;

  //TODO(Zak): update to allow modified correlation_length or slab_center_z.
  SfgMap sfg_map(water_group, kOH, params.dt());

  rvec *x_in = NULL;
  matrix box_mat;
  arma::rowvec box = arma::zeros<arma::rowvec>(DIMS);
  std::string xtc_filename = "prod.xtc";
  XDRFILE *xtc_file;
  params.ExtractTrajMetadata(strdup(xtc_filename.c_str()), (&x_in), box);
  xtc_file = xdrfile_open(strdup(xtc_filename.c_str()), "r");
  params.set_box(box);
  params.set_max_time(vm["max_time"].as<double>());

  hydrogen_group.OpenFieldFile();

  arma::rowvec dx;
  int st;
  float time, prec;
  int step = 0;
  for (step = 0; step < params.max_steps(); ++step) {
    if(read_xtc(xtc_file, params.num_atoms(), &st, &time, box_mat, x_in, &prec))
      break;

    params.set_box(box_mat);

    oxygen_group.set_positions(x_in);
    hydrogen_group.set_positions(x_in);
    if (!just_water)
      solute_group.set_positions(x_in);

    hydrogen_group.SetElectricField(water_group, params.box());
    if (!just_water)
      hydrogen_group.UpdateElectricField(solute_group, params.box());
    if (!hydrogen_group.field_check())
      hydrogen_group.WriteElectricField();

    sfg_map.CalculateFrequencies(hydrogen_group, oxygen_group, params.box());
  }
  xdrfile_close(xtc_file);

  sfg_map.CalculateSpectra();
  sfg_map.PrintSpectra("sfg_spectra.txt");

}
 // main
