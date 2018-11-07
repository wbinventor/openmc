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
void VolumeCalculation::check_hit(
    int i_domain, int i_material, openmc::int_2dvec indices,
    openmc::int_2dvec hits, std::vector<int> n_mat) {

  // Check if this material was previously hit and if so, increment count
  bool already_hit = false;
  for (int j=0; j < n_mat.size(); j++) {
    if (indices[i_domain][j] == i_material) {
      hits[i_domain][j]++;
      already_hit = true;
    }
  }

  // If the material was not previously hit, append an entry to the material
  // indices and hits lists 
  if (!already_hit) {
    indices[i_domain].push_back(i_material);
    hits[i_domain].push_back(1);
  }
}

// Stochastically determin the volume of a set of domains along with the
// average number densities of nuclides within the domain
// @param volumes volume mean/stdev in each domain
// @param i_nuclides indices in nuclides array
// @param n_atoms total # of atoms of each nuclide
// @param n_atoms_uncertainty uncertainty of total # of atoms
void VolumeCalculation::calculate_volumes(
    openmc::double_2dvec volumes, openmc::int_2dvec i_nuclides,
    openmc::double_2dvec n_atoms, openmc::double_2dvec n_atoms_uncertainty) {

  // Variables that are private to each thread
  int i_material;   // index in materials array
  std::vector<int> n_mat;   // Number of materials for each domain
  openmc::int_2dvec indices;  // List of material indices for each domain
  openmc::int_2dvec hits;      // Number of hits for each material in each domain
  bool found_cell;
  openmc::Particle* p;

  // Shared variables
  long int i_start, i_end;    // Starting/ending sample for each process
  openmc::int_2dvec master_indices(this->domain_ids.size());
  openmc::int_2dvec master_hits(this->domain_ids.size());

  // Variables used outside of parallel region
  int i_nuclide;     // Index in nuclides array
  int total_hits;    // Total hits for a single domain (summed over materials)
  int min_samples;   // Minimum number of samples per process
  int remainder;     // Leftover samples from uneven divide
#ifdef OPENMC_MPI
  int m;    // Index over materials
  int n;    // Number of materials
  std::vector<int> data;   // Array used to send number of hits
#endif

  double f;       // Fraction of hits
  double var_f;   // Variance of fraction of hits
  double volume_sample;  // Total volume of sampled region

  openmc::double_2dvec atoms(2);   // FIXME: Unclear what this is or how to initialize it???
  atoms[0].resize(nuclides.size());
  atoms[1].resize(nuclides.size());

  // FIXME: Do I need to define an interface for MPI's send_int and recv_int routines???

  // Divide work over MPI processes
  min_samples = this->n_samples / openmc::mpi::n_procs;
  remainder = this->n_samples % openmc::mpi::n_procs;
  if (openmc::mpi::rank < remainder) {
    i_start = (min_samples + 1) * openmc::mpi::rank;
    i_end = i_start + min_samples;
  }
  else {
    i_start = (min_samples + 1) * remainder + (openmc::mpi::rank - remainder) * min_samples;
    i_end = i_start + min_samples - 1;
  }

  particle_initialize(p);

#pragma omp parallel private(i_material, found_cell, indices, hits, n_mat) firstprivate(p)

  // Create space for material indices and number of hits for each
  n_mat.resize(this->domain_ids.size(), 0);
  indices.resize(this->domain_ids.size());
  hits.resize(this->domain_ids.size());
  for (int i=0; i < this->domain_ids.size(); i++) {
    indices[i].resize(8);
    hits[i].resize(8);
  }

  openmc::prn_set_stream(openmc::STREAM_VOLUME);

  // Samples locations and count hits

#pragma omp for
  for (long int i=i_start; i < i_end; i++) {
    openmc::set_particle_seed(i);

    p->n_coord = 1;
    p->coord[0].xyz[0] = this->lower_left[0] + openmc::prn() * 
      (this->upper_right[0] - this->lower_left[0]);
    p->coord[0].xyz[1] = this->lower_left[1] + openmc::prn() * 
      (this->upper_right[1] - this->lower_left[1]);
    p->coord[0].xyz[2] = this->lower_left[2] + openmc::prn() * 
      (this->upper_right[2] - this->lower_left[2]);
    p->coord[0].uvw[0] = 0.5;
    p->coord[1].uvw[1] = 0.5;
    p->coord[2].uvw[2] = 0.5;

    // If this location is not in the geometry at all, move on to next block
    found_cell = find_cell(p, 0);
    if (!found_cell) {
      continue;
    }

    if (this->domain_type == openmc::FILTER_MATERIAL) {
      i_material = p->material;
      if (i_material != openmc::MATERIAL_VOID) {
        for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
          if (openmc::materials[i_material]->id_ == this->domain_ids[i_domain]) {
            check_hit(i_domain, i_material, indices, hits, n_mat);
          }
        }
      }
    }
    else if (this->domain_type == openmc::FILTER_CELL) {
      for (int level=0; level < p->n_coord; level++) {
        for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
          if (cells[p->coord[level].cell + 1].id == this->domain_ids[i_domain]) {
            i_material = p->material;
            check_hit(i_domain, i_material, indices, hits, n_mat);
          }
        }
      }
    }
    else if (this->domain_type == openmc::FILTER_UNIVERSE) {
      for (int level=0; level < p->n_coord; level++) {
        for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
          if (universe_id[p->coord[level]].universe == this->domain_ids[i_domain]) {
            i_material = p->material;
            check_hit(i_domain, i_material, indices, hits, n_mat);
          }
        }
      }
    }
  }

  // Reduce hits onto master thread
  // At this point, each thread has its own pair of index/hits lists and we now
  // need to reduce them. OpenMP is not nearly smart enough to do this on its own,
  // so we have to manually reduce them

#ifdef _OPENMP
#pragma omp for ordered schedule(static)
  for (int i=0; i < omp_get_num_threads(); i++) {
#pragma omp ordered
    for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
      for (int j=0; j < n_mat[i_domain]; j++) {
        // Check if this material has been added to the master list and if so,
        // accumulate the number of hits
        for (int k=0; k < master_indices[i_domain].size(); k++) {
          if (indices[i_domain][j] == master_indices[i_domain][k]) {
            master_hits[i_domain][k] += hits[i_domain][j];
          }
        }
        // If we made it here, the material hasn't yet been added to the master
        // list, so add entries to the master indices and master hits lists
        master_indices[i_domain].push_back(indices[i_domain][j]);
        master_hits[i_domain].push_back(hits[i_domain][j]);
      }
    }
  }
#else
  for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
    for (int j=0; j < n_mat[i_domain]) {
      master_indices[i_domain].push_back(indices[i_domain][j]);
      master_hits[i_domain].push_back(hits[i_domain][j]);
    }
  }
#endif

  openmc::prn_set_stream(openmc::STREAM_TRACKING);

  // Reduce hits onto master process

  volume_sample = 1.;
  for (int i=0; i < 3; i++) {
    volume_sample *= (this->upper_right[i] - this->lower_left[i]);
  }

  for (int i_domain=0; i_domain < this->domain_ids.size(); i_domain++) {
    std::fill(atoms[0].begin(), atoms[0].end(), 0.);
    std::fill(atoms[1].begin(), atoms[1].end(), 0.);
    total_hits = 0;

    if (openmc::mpi::master) {
#ifdef OPENMC_MPI
      for (int j=0; j < openmc::mpi::n_procs - 1; j++) { 
        recv_int(&n, 1, j, 0);
        data.resize(2 * n);
        recv_int(&data[0], 2 * n, j, 1);
        for (int k=0; k < n-1; k++) {
          for (int m=0, master_indices[i_domain].size(); m++) {
            if (data[2 * k] == master_indices[i_domain][m]) {
              master_hits[i_domain][m] = 
                master_hits[i_domain][m] + data[2 * k + 1];
            }
          }
        }
        data.clear();
      }
#endif

      for (int j=0; j < master_indices[i_domain].size(); j++) {
        total_hits += master_hits[i_domain][j];
        f = double(master_hits[i_domain][j]) / this->n_samples;
        var_f = f * (1. - f) / this->n_samples;

        i_material = master_indices[i_domain][j];
        if (i_material == openmc::MATERIAL_VOID) {
          continue;
        }

        openmc::Material* mat = openmc::materials[i_material];
        // FIXME: Material class does not have a collection of nuclides yet!!
        for (int k=0; k < mat->nuclides.size(); k++) {
          // Accumulate nuclide density
          i_nuclide = mat->nuclides[k];
          atoms[0][i_nuclide] += mat->atom_density[k] * f;
          atoms[1][i_nuclide] += pow(mat->atom_density[k], 2) * var_f;
        }
      }

      // Determine volume
      volumes[i_domain][0] = double(total_hits) / this->n_samples *
        volume_sample;
      volumes[i_domain][1] = sqrt(volumes[i_domain][0] * (volume_sample - 
        volumes[i_domain][0]) / this->n_samples);

      // Determine total number of atoms. At this point, we have values in
      // atoms/b-cm. To get to atoms we multiply by 10^24 V.
      for (int j=0; j < atoms[0].size(); j++) {
        atoms[0][j] *= 1.e24 * volume_sample;
        atoms[1][j] *= 1.e24 * volume_sample;
      }

      // Convert full arrays to vectors
      // FIXME: Is this necessary in C++???
      for (int j=0; nuclides.size(); j++) {
        if (atoms[0][j] > 0.) {
          i_nuclides[i_domain].push_back(j);
          n_atoms[i_domain].push_back(atoms[0][j]);
          n_atoms_uncertainty[i_domain].push_back(atoms[1][j]);
        }
      }
    }
    else {
#ifdef OPENMC_MPI
      n = master_indices[i_domain].size();
      data.resize(2 * n);
      for (int k=0; k < n-1; k++) {
        data[2 * k] = master_indices[i_domain][k];
        data[2 * k + 1] = master_hits[i_domain][k];
      }

      send_int(&n, 1, 0, 0);
      send_int(&data[0], 2 * n, 0, 1);
      data.clear();
#endif
    }
  }
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
  file_id = openmc::file_open(filename, 'w');

  // Write header info
  openmc::write_attribute(file_id, "filetype", "volume");
  openmc::write_attribute(file_id, "version", openmc::VERSION_VOLUME);
  openmc::write_attribute(file_id, "openmc_version", openmc::VERSION);
#ifdef GIT_SHA1
  openmc::write_attribute(file_id, "git_sha1", GIT_SHA1);
#endif

  // Write current date and time
  // FIXME: Will this work -- time_stamp is a FORTRAN routine in output.F90???
  openmc::write_attribute(file_id, "date_and_time", openmc:: time_stamp());

  // Write basic metadata
  openmc::write_attribute(file_id, "n_samples", this->n_samples);
  openmc::write_attribute(file_id, "lower_left", this->lower_left);
  openmc::write_attribute(file_id, "upper_right", this->upper_right);
  if (this->domain_type == openmc::FILTER_CELL) {
    openmc::write_attribute(file_id, "domain_type", "cell");
  }
  else if (this->domain_type == openmc::FILTER_MATERIAL) {
    openmc::write_attribute(file_id, "domain_type", "material");
  }
  else if (this->domain_type == openmc::FILTER_UNIVERSE) {
    openmc::write_attribute(file_id, "domain_type", "universe");
  }

  for (int i=0; i < this->domain_ids.size(); i++)
  {
    group_id = openmc::create_group(file_id, "domain_" + std::to_string(this->domain_ids[i]));

    // Write volume for domain
    // FIXME: openmc::double_2dvec won't allow for index slicing; maybe use xtensor??
    openmc::write_dataset(group_id, "volume", volumes[i]);

    // Create array of nuclide names from the vector
    n = i_nuclides[i].size();
    if (n > 0) {
      nucnames.resize(n);
      for (int j=0; j < n; j++) {
        // NOTE: i_nuclides is a 2d array indexed by domain id and nuclide number
        // i.e., each nuclide in each domain for a specific volume calculation
        nucnames[j] = nuclides[i_nuclides[i]][j].name;
      }

      // Create array of total # of atoms with uncertainty for each nuclide
      atom_data.push_back(n_atoms[i]);
      atom_data.push_back(n_atoms_uncertainty[i]);

      // Write results
      openmc::write_dataset(group_id, "nuclides", nucnames);
      openmc::write_dataset(group_id, "atoms", atom_data);

      nucnames.clear();
      atom_data.clear();
    }

    openmc::close_group(group_id);
  }

  openmc::file_close(file_id);
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
  openmc::int_2dvec i_nuclides;   // indices in nuclides array
  openmc::double_2dvec n_atoms;   // total # of atoms of each nuclide
  openmc::double_2dvec n_atoms_uncertainty;  // uncertainty of total # of atoms

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
    volumes.resize(n_domains);
    for (int j = 0; j < n_domains; j++) {
      volumes[j].resize(2);
    }

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
            << volumes[j][0] << " +/- " << volumes[j][1] << " cm^3";
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



