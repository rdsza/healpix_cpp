#ifndef HEALPIX_SEARCH_H
#define HEALPIX_SEARCH_H

#include "src/Healpix_2.15a/healpix_base.h"
#include "src/Healpix_2.15a/pointing.h"
#include <vector>

// Forward declarations for compatibility with existing code
class ParameterMap;
class wxString;

class HealpixSearch {
public:
    // Public members - maintaining compatibility with EulerSearch interface
    int     number_of_search_positions;
    int     best_parameters_to_keep;
    float** list_of_search_parameters;  // [position][0] = phi, [1] = theta
    float** list_of_best_parameters;
    
    // HEALPix specific parameters
    int     nside;
    float   angular_step_size;
    float   resolution_limit;
    bool    test_mirror;
    float   theta_max;
    
    // Sampling method (0 = Fibonacci, 1 = Hopf - using Fibonacci as default)
    int     sampling_method;
    
    // Parameter map for compatibility
    ParameterMap* parameter_map;
    
    // Symmetry symbol for compatibility (not used in uniform coverage)
    wxString* symmetry_symbol;
    
    // Constructor and destructor
    HealpixSearch();
    HealpixSearch(const HealpixSearch& other_search);  // copy constructor
    ~HealpixSearch();
    
    // Assignment operators
    HealpixSearch& operator=(const HealpixSearch& other_search);
    HealpixSearch& operator=(const HealpixSearch* other_search);
    
    // Initialization methods - maintaining EulerSearch interface
    void Init(float wanted_resolution_limit, ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep);
    void InitGrid(const wxString& wanted_symmetry_symbol, float wanted_angular_step_size, float wanted_phi_start, 
                  float wanted_theta_start, float wanted_psi_max, float wanted_psi_step, float wanted_psi_start, 
                  float wanted_resolution_limit, ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep);
    
    // HEALPix specific initialization
    void InitHealpix(int wanted_nside, float wanted_angular_step_size, float wanted_resolution_limit, 
                     ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep, int wanted_sampling_method = 0);
    
    // Search position calculation methods
    void CalculateGridSearchPositions(bool random_start_angle = false);
    void CalculateHealpixSearchPositions();
    
    // Symmetry methods (for compatibility - not used in uniform coverage)
    void SetSymmetryLimits();

private:
    // Internal HEALPix generation methods
    std::vector<pointing> GenerateFibonacciPoints(int N);
    std::vector<pointing> GenerateHopfPoints(int N);
    std::vector<pointing> FilterPointsInPixel(const std::vector<pointing>& all_points, 
                                             const pointing& pixel_center, 
                                             double angular_step);
    double AngularDistance(const pointing& a, const pointing& b);
    std::vector<pointing> SelectPointsForPixelWithMinSpacing(const std::vector<pointing>& candidate_points,
                                                           const Healpix_Base& healpix_base,
                                                           long pixel_index,
                                                           double min_sep);
    
    // Memory management helpers
    void DeallocateArrays();
    void AllocateSearchArrays(int num_positions);
    void AllocateBestParametersArray();
};

#endif // HEALPIX_SEARCH_H
