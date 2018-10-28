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
    openmc::double_2dvec volumes, openmc::int_2dvec i_nuclides,
    openmc::double_2dvec n_atoms, openmc::double_2dvec n_atoms_uncertainty) {
  return;
}

void VolumeCalculation::write_volume(
  std::string filename, openmc::double_2dvec volumes, openmc::int_2dvec i_nuclides,
  openmc::double_2dvec n_atoms, openmc::double_2dvec n_atoms_uncertainty) {

  int i, j;
  int n;
  hid_t file_id;
  hid_t group_id;
  openmc::double_2dvec atom_data;  // mean/stdev of total # of atoms for each nuclide
  std::vector<std::string> nucnames;  // names of nuclides

  // Create HDF5 file
  file_id = file_open(filename, 'w');

  // Write header info
  write_attribute(fild_id, "filetype", "volume");
  write_attribute(file_id, "version", VERSION_VOLUME);
  write_attribute(file_id, "openmc_version", VERSION);
#ifdef GIT_SHA1
  write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  // FIXME: Will this work -- time_stamp is a FORTRAN routine in output.F90???
  write_attribute(file_id, "date_and_time", time_stamp());

  // Write basic metadata
  write_attribute(file_id, "n_samples", this->n_samples);
  write_attribute(file_id, "lower_left", this->lower_left);
  write_attribute(file_id, "upper_right", this->upper_right);
  if (this->domain_type == FILTER_CELL) {
    write_attribute(file_id, "domain_type", "cell");
  }
  else if (this->domain_type == FILTER_MATERIAL) {
    write_attribute(file_id, "domain_type", "material");
  }
  else if (this->domain_type == FILTER_UNIVERSE) {
    write_attribute(file_id, "domain_type", "universe");
  }

  for (int i=0; i < this->domain_ids.size(); i++)
  {
    group_id = create_group(file_id, "domain_" + std::to_string(this->domain_ids[i]));

    // Write volume for domain
    // FIXME: openmc::double_2dvec won't allow for index slicing; maybe use xtensor??
    write_dataset(group_id, "volume", volumes[:,i]);

    // Create array of nuclide names from the vector
    n = i_nuclides[i];
    if (n > 0) {
      nucnames.resize(n);
      for (int j=0; j < n; j++) {
        // NOTE: i_nuclides is a 2d array indexed by domain id and nuclide number
        // i.e., each nuclide in each domain for a specific volume calculation
        nucnames[j] = nuclides[i_nuclides[i]][j].name;
      }

      // Create array of total # of atoms with uncertainty for each nuclide
      atom_data.push_back(atoms_vec[i]);
      atom_data.push_back(uncertainty_vec[i]);

      // Write results
      write_dataset(group_id, "nuclides", nucnames);
      write_dataset(group_id, "atoms", atom_data);

      nucnames.clear();
      atom_data.clear();
    }

    close_group(group_id);
  }

  file_close(file_id);
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



