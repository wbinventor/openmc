#include "openmc/volume_calc.h"

//==============================================================================
// Global variables
//==============================================================================

std::vector<VolumeCalculation> volume_calcs;


//==============================================================================
// VolumeCalculation implementation
//==============================================================================

VolumeCalculation::VolumeCalculation(pugi::xml_node node) {

  // NOTE: get_node_value(...) takes optional bool args lowercase and strip (False by default)

  // Read domain type (cell, material or universe)
  std::string domain_type = openmc::get_node_value(node, "domain_type");
  if (domain_type == "cell") {
    this->domain_type = openmc::FILTER_CELL;
  }
  else if (domain_type == "material") {
    this->domain_type = openmc::FILTER_MATERIAL;
  }
  else if (domain_type == "universe") {
    this->domain_type = openmc::FILTER_UNIVERSE;
  }
  else {
    openmc::fatal_error(std::string("Unrecognized domain type for stochastic "
           "volume calculation:" + domain_type));
  }

  // Read domain IDs, bounding corodinates and number of samples
  this->domain_ids = openmc::get_node_array<int>(node, "domain_ids");
  this->lower_left = openmc::get_node_array<double>(node, "lower_left");
  this->upper_right = openmc::get_node_array<double>(node, "upper_right");
  this->n_samples = std::stoi(openmc::get_node_value(node, "samples"));
}

//VolumeCalculation::~VolumeCalculation() {
  // FIXME: deallocate memory for arrays; i.e., volume_calcs
//}


  // FIXME: This may ned to take more parameters...
void VolumeCalculation::check_hit(int domain, int material) {
  return;
}


void VolumeCalculation::calculate_volumes(
    openmc::double_2dvec volumes, std::vector<int> i_nuclides,
    std::vector<double> n_atoms, std::vector<double> n_atoms_uncertainty) {
  return;
}

void VolumeCalculation::write_volume(
  std::string filename, openmc::double_2dvec volumes, std::vector<int> i_nuclides,
  std::vector<double> n_atoms, std::vector<double> n_atoms_uncertainty) {
  return;
}


//==============================================================================
// OPENMC_CALCULATE_VOLUMES runs each of the stochastic volume calculations
// that the user has specified and writes results to HDF5 files
//==============================================================================

int openmc_calculate_volumes() {

  int err = 0;
  int i, j, n_domains;
  openmc::double_2dvec volumes;  // volume mean/stdev in each domain
  std::string domain_type;
  std::string filename;    // filename for HDF5 file
  // FIXME: Initialize a timer...but is there one for C++ yet?
  std::vector<int> i_nuclides;   // indices in nuclides array
  std::vector<double> n_atoms;   // total # of atoms of each nuclide
  std::vector<double> n_atoms_uncertainty;  // uncertainty of total # of atoms

  // FIXME: Does this work for MPI in C++???
  if (openmc::mpi::master) {
    openmc::write_message("STOCHASTIC VOLUME CALCULATION", 3);
    // FIXME: Start the timer here....
  }

  // FIXME: Where is "volume_calcs" coming from???
  for (int i = 0; i < volume_calcs.size(); i++) {
    n_domains = volume_calcs[i].domain_ids.size();
    i_nuclides.resize(n_domains);
    n_atoms.resize(n_domains);
    n_atoms_uncertainty.resize(n_domains);

    // FIXME: This probably won't work for a nested std::vector
    volumes.resize(2);
    volumes[0].resize(n_domains);
    volumes[1].resize(n_domains);

    if (openmc::mpi::master) {
      std::stringstream msg;
      msg << "Running volume calculation " << i << "...";
      openmc::write_message(msg, 4);
    }

    // TODO: Implement this such that each MPI proc separatel computes the volume
    // and this function then reduces them at the end
    volume_calcs[i].calculate_volumes(
      volumes, i_nuclides, n_atoms, n_atoms_uncertainty);

    // TODO: Implement MPI reduction across volumes, atoms and nuclides

    if (openmc::mpi::master) {
      if (volume_calcs[i].domain_type == openmc::FILTER_CELL) {
        domain_type = "  Cell";
      }
      else if (volume_calcs[i].domain_type == openmc::FILTER_MATERIAL) {
        domain_type = "  Material";
      }
      else {
        domain_type = "  Universe";
      }

      // Display domain volumes
      for (int j = 0; j < volume_calcs[i].domain_ids.size(); j++) {
        std::stringstream msg;
        msg << domain_type << " " << volume_calcs[j].domain_ids[j]
            << volumes[0][j] << " +/- " << volumes[0][j] << " cm^3";
        openmc::write_message(msg.str(), 4);
      }
    
      // Write volumes to HDF5 file
      // FIXME: Can VolumeCalculation method take path and filename as args?
      // FIXME: Should this be a static class method, since it's only called by the master proc?
      std::stringstream filename;
      filename << openmc::settings::path_output << "volume_" << i << ".h5";
      volume_calcs[i].write_volume(
        filename.str(), volumes, i_nuclides, n_atoms, n_atoms_uncertainty);

    }

    // FIXME: Need to deallocate any std::vectors???

    // FIXME: Show elapsed time; is there a C++ timer yet???
  }

  err = 0;
  return err;
}



