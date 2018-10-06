#ifndef VOLUME_CALC_H
#define VOLUME_CALC_H

#include <string>
#include <vector>
#include <map>

#ifdef _OPENMP
  #include <omp.h>
#endif

#include "openmc/constants.h" // FILTER_CELL, FILTER_MATERIAL, FILTER_UNIVERSE
#include "openmc/xml_interface.h"  // ???
#include "openmc/error.h"  // fatal_error(...), write_message(...)
#include "pugixml.hpp"  // XML
#include "openmc/geometry.h" // find_cell(...)
#include "openmc/hdf5_interface.h" // file_open, ...
#include "openmc/material.h" // std::vector<Material*> materials
#include "openmc/message_passing.h"
// nuclide_header ???
#include "openmc/random_lcg.h" // prn, prn_set_stream, set_particle_seed
#include "openmc/settings.h" // path_output
#include "openmc/capi.h"  // openmc_get_volumes() declaration
// timer_header ???


class VolumeCalculation {
public:
  // Constructors
  VolumeCalculation() = default;
  VolumeCalculation(pugi::xml_node node);

//  virtual ~VolumeCalculation();  // Needed??
  int domain_type;
  std::vector<int> domain_ids;
  std::vector<double> lower_left;
  std::vector<double> upper_right;
  int n_samples;

  /* ???
  std::map<int, std::vector<int> > hits;
  std::vector<int> volumes;
  std::vector<int> nuclide_vec;
  std::vector<double> atoms_vec;
  std::vector<double> uncertainty_vec;
  */

  // FIXME: This may ned to take more parameters...
  void check_hit(int domain, int material);

  void set_domain_type(int domain_type);
  void calculate_volumes(
      openmc::double_2dvec volumes, std::vector<int> i_nuclides,
      std::vector<double> n_atoms, std::vector<double> n_atoms_uncertainty);
  void write_volume(
    std::string filename, openmc::double_2dvec volumes, std::vector<int> i_nuclides,
    std::vector<double> n_atoms, std::vector<double> n_atoms_uncertainty);

};


#endif // VOLUME_CALC_H