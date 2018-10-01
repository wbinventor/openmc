#include <openmc/volume_calc.h>


void VolumeCalculation::calculate_volume() {

}

void VolumeCalculation::write_volume(std::string filename, double volume,
                                     std::vector<int> nuclide_vec,
                                     std::vector<double> atoms_vec,
                                     std::vector<double> uncertainty_vec) {

}

void VolumeCalculation::check_hit(int domain, int material) {

}


void VolumeCalculation::from_xml(xml_node_struct* xml_node) {
  // TODO: Look at XML CPP API...

  // FIXME: Check domain type (cell, material or universe)
  std::string domain_type = get_node_value(xml_node, "domain_type");

  if (domain_type == "cell") {
    this.domain_type = FILTER_CELL;
  }
  else if (domain_type == "material") {
    this->domain_type = FILTER_MATERIAL;
  }
  else if (domain_type == "universe") {
    this->domain_type = FILTER_UNIVERSE;
  }
  else {
    std::stringstream err_msg;
    err_msg << "Unrecognized domain type for stochastic volume "
            << "calculation: \"" << domain_type << "\"";
    fatal_error(err_msg);
  }
}

void openmc_calculate_volumes() {

}