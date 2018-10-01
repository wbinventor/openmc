#ifndef VOLUME_CALC_H
#define VOLUME_CALC_H

#include <string>
#include <vector>
#include <map>

#ifdef _OPENMP
  #incude <omp.h>
#endif

#include <openmc/constants.h>
#include <openmc/xml_interface.h>


class VolumeCalculation {
private:
    int domain_type;
    std::vector<int> domain_id;
    double lower_left[3];
    double upper_right[3];
    int samples;
    std::map<int, std::vector<int> > hits;

    std::vector<int> volumes;
    std::vector<int> nuclide_vec;
    std::vector<double> atoms_vec;
    std::vector<double> uncertainty_vec;

    // FIXME: This may ned to take more parameters...
    void check_hit(int domain, int material);

public:
    void calculate_volume();
    void write_volume(std::string filename, double volume,
                      std::vector<int> nuclide_vec,
                      std::vector<double> atoms_vec,
                      std::vector<double> uncertainty_vec);
    void from_xml(xml_node_struct* xml_node);
}


// FIXME: Should I make this a static class method for VolumeCalculation?
VolumeCalculation* volume_from_xml(xml_node_struct* xml_node);

void openmc_calculate_volumes();

// FIXME:  Make volume_cals (collection of VolumeCalculation) a global variable???

#endif // VOLUME_CALC_H