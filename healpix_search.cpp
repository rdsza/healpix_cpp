#include "healpix_search.h"
#include <cmath>
#include <random>
#include <algorithm>
#include <iostream>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Constructor
HealpixSearch::HealpixSearch() {
    number_of_search_positions = 0;
    best_parameters_to_keep = 0;
    list_of_search_parameters = nullptr;
    list_of_best_parameters = nullptr;
    nside = 1;
    angular_step_size = 0.0f;
    resolution_limit = 0.0f;
    test_mirror = false;
    theta_max = 180.0f;
    sampling_method = 0; // Fibonacci as default
    parameter_map = nullptr;
    symmetry_symbol = nullptr;
}

// Copy constructor
HealpixSearch::HealpixSearch(const HealpixSearch& other_search) {
    list_of_search_parameters = nullptr;
    list_of_best_parameters = nullptr;
    parameter_map = nullptr;
    symmetry_symbol = nullptr;
    *this = other_search;
}

// Destructor
HealpixSearch::~HealpixSearch() {
    DeallocateArrays();
    if (parameter_map) {
        delete parameter_map;
        parameter_map = nullptr;
    }
    if (symmetry_symbol) {
        delete symmetry_symbol;
        symmetry_symbol = nullptr;
    }
}

// Assignment operators
HealpixSearch& HealpixSearch::operator=(const HealpixSearch& other_search) {
    return *this = &other_search;
}

HealpixSearch& HealpixSearch::operator=(const HealpixSearch* other_search) {
    if (this != other_search) {
        // Copy basic parameters
        number_of_search_positions = other_search->number_of_search_positions;
        best_parameters_to_keep = other_search->best_parameters_to_keep;
        nside = other_search->nside;
        angular_step_size = other_search->angular_step_size;
        resolution_limit = other_search->resolution_limit;
        test_mirror = other_search->test_mirror;
        theta_max = other_search->theta_max;
        sampling_method = other_search->sampling_method;
        
        // Handle parameter_map
        if (parameter_map) {
            delete parameter_map;
            parameter_map = nullptr;
        }
        if (other_search->parameter_map) {
            // We'll handle this in the implementation where we know the type
            parameter_map = other_search->parameter_map; // shallow copy for now
        }
        
        // Handle symmetry_symbol
        if (symmetry_symbol) {
            delete symmetry_symbol;
            symmetry_symbol = nullptr;
        }
        if (other_search->symmetry_symbol) {
            // We'll handle this in the implementation where we know the type
            symmetry_symbol = other_search->symmetry_symbol; // shallow copy for now
        }
        
        // Copy search parameters array
        DeallocateArrays();
        if (other_search->list_of_search_parameters && number_of_search_positions > 0) {
            AllocateSearchArrays(number_of_search_positions);
            for (int i = 0; i < number_of_search_positions; i++) {
                list_of_search_parameters[i][0] = other_search->list_of_search_parameters[i][0];
                list_of_search_parameters[i][1] = other_search->list_of_search_parameters[i][1];
            }
        }
        
        // Copy best parameters array
        if (other_search->list_of_best_parameters && best_parameters_to_keep > 0) {
            AllocateBestParametersArray();
            for (int i = 0; i <= best_parameters_to_keep; i++) {
                for (int j = 0; j < 6; j++) {
                    list_of_best_parameters[i][j] = other_search->list_of_best_parameters[i][j];
                }
            }
        }
    }
    return *this;
}

// Basic initialization
void HealpixSearch::Init(float wanted_resolution_limit, ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep) {
    resolution_limit = wanted_resolution_limit;
    parameter_map = &wanted_parameter_map;
    best_parameters_to_keep = wanted_parameters_to_keep;
    
    if (best_parameters_to_keep > 0) {
        AllocateBestParametersArray();
    }
}

// Grid initialization (for compatibility with EulerSearch interface)
void HealpixSearch::InitGrid(const wxString& wanted_symmetry_symbol, float wanted_angular_step_size, 
                           float wanted_phi_start, float wanted_theta_start, float wanted_psi_max, 
                           float wanted_psi_step, float wanted_psi_start, float wanted_resolution_limit, 
                           ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep) {
    
    if (number_of_search_positions == 0) {
        Init(wanted_resolution_limit, wanted_parameter_map, wanted_parameters_to_keep);
    }
    
    angular_step_size = wanted_angular_step_size;
    symmetry_symbol = const_cast<wxString*>(&wanted_symmetry_symbol);
    
    // Set basic parameters for compatibility
    test_mirror = false; // HEALPix generates uniform coverage, no mirroring needed
    theta_max = 180.0f;
    
    // For HEALPix, we'll calculate positions in CalculateHealpixSearchPositions
    // This method maintains interface compatibility but uses HEALPix internally
}

// HEALPix specific initialization
void HealpixSearch::InitHealpix(int wanted_nside, float wanted_angular_step_size, float wanted_resolution_limit, 
                               ParameterMap& wanted_parameter_map, int wanted_parameters_to_keep, int wanted_sampling_method) {
    nside = wanted_nside;
    angular_step_size = wanted_angular_step_size;
    resolution_limit = wanted_resolution_limit;
    parameter_map = &wanted_parameter_map;
    best_parameters_to_keep = wanted_parameters_to_keep;
    sampling_method = wanted_sampling_method;
    test_mirror = false;
    theta_max = 180.0f;
    
    if (best_parameters_to_keep > 0) {
        AllocateBestParametersArray();
    }
}

// Grid search position calculation (delegates to HEALPix for this implementation)
void HealpixSearch::CalculateGridSearchPositions(bool random_start_angle) {
    CalculateHealpixSearchPositions();
}

// Main HEALPix search position calculation
void HealpixSearch::CalculateHealpixSearchPositions() {
    if (nside <= 0) {
        std::cerr << "Error: nside must be positive for HEALPix search" << std::endl;
        return;
    }
    
    long npix = 12 * static_cast<long>(nside) * static_cast<long>(nside);
    Healpix_Base healpix_base(nside, RING, SET_NSIDE);
    
    // Convert angular step from degrees to radians
    double angular_step_rad = angular_step_size * (M_PI / 180.0);
    
    // Calculate area per point and points per pixel
    const double area_per_point = (std::sqrt(3.0) / 2.0) * angular_step_rad * angular_step_rad;
    const double pixel_area = (4.0 * M_PI) / static_cast<double>(npix);
    int points_per_pixel = std::max(1, static_cast<int>(std::floor(pixel_area / area_per_point)));
    
    // Generate candidate points
    int num_candidates = points_per_pixel * npix * 3; // oversample for safety
    std::vector<pointing> all_points;
    
    if (sampling_method == 0) {
        all_points = GenerateFibonacciPoints(num_candidates);
    } else {
        all_points = GenerateHopfPoints(num_candidates);
    }
    
    // Deterministic shuffle to avoid ordering bias
    const uint64_t seed = static_cast<uint64_t>(nside) * 1315423911ULL ^
                          static_cast<uint64_t>(sampling_method) * 2654435761ULL ^
                          static_cast<uint64_t>(std::llround(angular_step_rad * 1e9));
    std::mt19937_64 rng(seed);
    std::shuffle(all_points.begin(), all_points.end(), rng);
    
    // Assign points to pixels
    std::vector<std::vector<pointing>> pixel_points(npix);
    for (const auto& pt : all_points) {
        long pix = healpix_base.ang2pix(pt);
        if (pixel_points[pix].size() < static_cast<size_t>(points_per_pixel)) {
            // Enforce minimum separation within pixel
            bool far_enough = true;
            for (const auto& keep : pixel_points[pix]) {
                if (AngularDistance(pt, keep) < angular_step_rad) {
                    far_enough = false;
                    break;
                }
            }
            if (far_enough) {
                pixel_points[pix].push_back(pt);
            }
        }
    }
    
    // If any pixel has fewer than points_per_pixel, fill with closest available points
    for (long pix = 0; pix < npix; ++pix) {
        if (pixel_points[pix].size() < static_cast<size_t>(points_per_pixel)) {
            for (const auto& pt : all_points) {
                if (healpix_base.ang2pix(pt) == pix) {
                    pixel_points[pix].push_back(pt);
                    if (pixel_points[pix].size() == static_cast<size_t>(points_per_pixel)) break;
                }
            }
        }
    }
    
    // Count total positions
    int total_positions = 0;
    for (long pix = 0; pix < npix; ++pix) {
        total_positions += std::min(static_cast<int>(pixel_points[pix].size()), points_per_pixel);
    }
    
    // Allocate arrays
    DeallocateArrays();
    number_of_search_positions = total_positions;
    AllocateSearchArrays(number_of_search_positions);
    
    // Fill the search parameters array
    int position_index = 0;
    for (long pix = 0; pix < npix; ++pix) {
        const auto& samples = pixel_points[pix];
        for (size_t i = 0; i < samples.size() && i < static_cast<size_t>(points_per_pixel); ++i) {
            // Convert from radians to degrees for compatibility
            double phi_deg = samples[i].phi * (180.0 / M_PI);
            double theta_deg = samples[i].theta * (180.0 / M_PI);
            
            // Round to reasonable precision
            phi_deg = std::round(phi_deg * 10000.0) / 10000.0;
            theta_deg = std::round(theta_deg * 10000.0) / 10000.0;
            
            list_of_search_parameters[position_index][0] = static_cast<float>(phi_deg);
            list_of_search_parameters[position_index][1] = static_cast<float>(theta_deg);
            position_index++;
        }
    }
    
    // Adjust if we have fewer positions than expected
    number_of_search_positions = position_index;
    
    // Adjust best_parameters_to_keep if necessary
    if (number_of_search_positions < best_parameters_to_keep) {
        best_parameters_to_keep = number_of_search_positions;
        AllocateBestParametersArray(); // Reallocate with correct size
    }
}

// Symmetry limits (for compatibility - not used in uniform coverage)
void HealpixSearch::SetSymmetryLimits() {
    // HEALPix generates uniform coverage, so we don't limit based on symmetry
    // This method is kept for interface compatibility
    test_mirror = false;
    theta_max = 180.0f;
}

// Generate Fibonacci sphere sampling points
std::vector<pointing> HealpixSearch::GenerateFibonacciPoints(int N) {
    std::vector<pointing> points;
    double phi = M_PI * (3.0 - std::sqrt(5.0)); // Golden angle
    
    for (int i = 0; i < N; i++) {
        double y = 1.0 - (i / static_cast<double>(N - 1)) * 2.0;
        double radius = std::sqrt(1.0 - y * y);
        double theta = phi * i;
        
        double x = std::cos(theta) * radius;
        double z = std::sin(theta) * radius;
        
        // Convert to spherical coordinates
        double theta_sph = std::acos(z);
        double phi_sph = std::atan2(y, x);
        if (phi_sph < 0) phi_sph += 2 * M_PI;
        
        points.emplace_back(theta_sph, phi_sph);
    }
    return points;
}

// Generate Hopf fibration points
std::vector<pointing> HealpixSearch::GenerateHopfPoints(int N) {
    std::vector<pointing> points;
    
    // For Hopf fibration, we need to choose good parameters
    int n_phi = static_cast<int>(std::sqrt(N));
    int n_theta = N / n_phi;
    
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            double theta = M_PI * (i + 0.5) / n_theta;
            double phi = 2 * M_PI * j / n_phi;
            
            // Apply Hopf transformation
            double psi = 2 * M_PI * (i * n_phi + j) / N;
            
            // Convert to Cartesian coordinates on S3
            double x1 = std::cos(theta/2) * std::cos(psi);
            double x2 = std::cos(theta/2) * std::sin(psi);
            double x3 = std::sin(theta/2) * std::cos(phi + psi);
            double x4 = std::sin(theta/2) * std::sin(phi + psi);
            
            // Project to S2
            double theta_sph = std::acos(x3);
            double phi_sph = std::atan2(x4, std::sqrt(x1*x1 + x2*x2));
            if (phi_sph < 0) phi_sph += 2 * M_PI;
            
            points.emplace_back(theta_sph, phi_sph);
        }
    }
    
    // If we need more points, add some with different parameters
    while (points.size() < static_cast<size_t>(N)) {
        double theta = M_PI * (points.size() + 0.5) / N;
        double phi = 2 * M_PI * (points.size() % N) / N;
        points.emplace_back(theta, phi);
    }
    
    return points;
}

// Filter points within pixel boundaries (not used in current implementation but kept for completeness)
std::vector<pointing> HealpixSearch::FilterPointsInPixel(const std::vector<pointing>& all_points, 
                                                        const pointing& pixel_center, 
                                                        double angular_step) {
    std::vector<pointing> filtered_points;
    
    for (const auto& pt : all_points) {
        double angular_dist = AngularDistance(pt, pixel_center);
        if (angular_dist <= angular_step) {
            filtered_points.push_back(pt);
        }
    }
    
    return filtered_points;
}

// Compute great-circle angular distance between two directions (radians)
double HealpixSearch::AngularDistance(const pointing& a, const pointing& b) {
    const double cos_dist = std::cos(a.theta) * std::cos(b.theta) +
                            std::sin(a.theta) * std::sin(b.theta) *
                            std::cos(a.phi - b.phi);
    return std::acos(std::max(-1.0, std::min(1.0, cos_dist)));
}

// Select points that fall inside a given pixel and are at least min_sep apart (not used but kept for completeness)
std::vector<pointing> HealpixSearch::SelectPointsForPixelWithMinSpacing(
    const std::vector<pointing>& candidate_points,
    const Healpix_Base& healpix_base,
    long pixel_index,
    double min_sep) {
    
    std::vector<pointing> accepted;
    accepted.reserve(64);

    for (const pointing& pt : candidate_points) {
        // Keep only candidates that belong to this pixel
        if (healpix_base.ang2pix(pt) != pixel_index) continue;

        bool far_enough = true;
        for (const pointing& keep : accepted) {
            if (AngularDistance(pt, keep) < min_sep) {
                far_enough = false;
                break;
            }
        }
        if (far_enough) {
            accepted.push_back(pt);
        }
    }

    return accepted;
}

// Memory management helpers
void HealpixSearch::DeallocateArrays() {
    if (list_of_search_parameters) {
        for (int i = 0; i < number_of_search_positions; i++) {
            delete[] list_of_search_parameters[i];
        }
        delete[] list_of_search_parameters;
        list_of_search_parameters = nullptr;
    }
    
    if (list_of_best_parameters) {
        for (int i = 0; i <= best_parameters_to_keep; i++) {
            delete[] list_of_best_parameters[i];
        }
        delete[] list_of_best_parameters;
        list_of_best_parameters = nullptr;
    }
}

void HealpixSearch::AllocateSearchArrays(int num_positions) {
    if (num_positions <= 0) return;
    
    list_of_search_parameters = new float*[num_positions];
    for (int i = 0; i < num_positions; i++) {
        list_of_search_parameters[i] = new float[2];
    }
}

void HealpixSearch::AllocateBestParametersArray() {
    if (best_parameters_to_keep <= 0) return;
    
    if (list_of_best_parameters) {
        for (int i = 0; i <= best_parameters_to_keep; i++) {
            delete[] list_of_best_parameters[i];
        }
        delete[] list_of_best_parameters;
    }
    
    list_of_best_parameters = new float*[best_parameters_to_keep + 1];
    for (int i = 0; i <= best_parameters_to_keep; i++) {
        list_of_best_parameters[i] = new float[6];
        // Initialize to safe values
        for (int j = 0; j < 6; j++) {
            list_of_best_parameters[i][j] = 0.0f;
        }
    }
}
