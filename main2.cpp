#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <sstream>
#include <vector>
#include <algorithm>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include "src/Healpix_2.15a/healpix_base.h"
#include "src/Healpix_2.15a/pointing.h"

// Function to generate Fibonacci sphere sampling points
std::vector<pointing> generate_fibonacci_points(int N) {
    std::vector<pointing> points;
    double phi = M_PI * (3.0 - sqrt(5.0)); // Golden angle
    
    for (int i = 0; i < N; i++) {
        double y = 1.0 - (i / (double)(N - 1)) * 2.0;
        double radius = sqrt(1.0 - y * y);
        double theta = phi * i;
        
        double x = cos(theta) * radius;
        double z = sin(theta) * radius;
        
        // Convert to spherical coordinates
        double theta_sph = acos(z);
        double phi_sph = atan2(y, x);
        if (phi_sph < 0) phi_sph += 2 * M_PI;
        
        points.emplace_back(theta_sph, phi_sph);
    }
    return points;
}

// Function to generate Hopf fibration points
std::vector<pointing> generate_hopf_points(int N) {
    std::vector<pointing> points;
    
    // For Hopf fibration, we need to choose good parameters
    int n_phi = (int)sqrt(N);
    int n_theta = N / n_phi;
    
    for (int i = 0; i < n_theta; i++) {
        for (int j = 0; j < n_phi; j++) {
            double theta = M_PI * (i + 0.5) / n_theta;
            double phi = 2 * M_PI * j / n_phi;
            
            // Apply Hopf transformation
            double psi = 2 * M_PI * (i * n_phi + j) / N;
            
            // Convert to Cartesian coordinates on S3
            double x1 = cos(theta/2) * cos(psi);
            double x2 = cos(theta/2) * sin(psi);
            double x3 = sin(theta/2) * cos(phi + psi);
            double x4 = sin(theta/2) * sin(phi + psi);
            
            // Project to S2
            double theta_sph = acos(x3);
            double phi_sph = atan2(x4, sqrt(x1*x1 + x2*x2));
            if (phi_sph < 0) phi_sph += 2 * M_PI;
            
            points.emplace_back(theta_sph, phi_sph);
        }
    }
    
    // If we need more points, add some with different parameters
    while (points.size() < N) {
        double theta = M_PI * (points.size() + 0.5) / N;
        double phi = 2 * M_PI * (points.size() % N) / N;
        points.emplace_back(theta, phi);
    }
    
    return points;
}

// Function to filter points within pixel boundaries
std::vector<pointing> filter_points_in_pixel(const std::vector<pointing>& all_points, 
                                           const pointing& pixel_center, 
                                           double angular_step) {
    std::vector<pointing> filtered_points;
    
    for (const auto& pt : all_points) {
        // Calculate angular distance between point and pixel center
        double cos_dist = cos(pt.theta) * cos(pixel_center.theta) + 
                         sin(pt.theta) * sin(pixel_center.theta) * 
                         cos(pt.phi - pixel_center.phi);
        
        double angular_dist = acos(std::max(-1.0, std::min(1.0, cos_dist)));
        
        // Keep points within the angular step
        if (angular_dist <= angular_step) {
            filtered_points.push_back(pt);
        }
    }
    
    return filtered_points;
}

// Compute great-circle angular distance between two directions (radians)
static inline double angular_distance(const pointing& a, const pointing& b) {
    const double cos_dist = std::cos(a.theta) * std::cos(b.theta) +
                            std::sin(a.theta) * std::sin(b.theta) *
                            std::cos(a.phi - b.phi);
    return std::acos(std::max(-1.0, std::min(1.0, cos_dist)));
}

// Select points that fall inside a given pixel and are at least min_sep apart from each other
std::vector<pointing> select_points_for_pixel_with_min_spacing(
    const std::vector<pointing>& candidate_points,
    const Healpix_Base& healpix_base,
    long pixel_index,
    double min_sep)
{
    std::vector<pointing> accepted;
    accepted.reserve(64);

    for (const pointing& pt : candidate_points) {
        // Keep only candidates that belong to this pixel
        if (healpix_base.ang2pix(pt) != pixel_index) continue;

        bool far_enough = true;
        for (const pointing& keep : accepted) {
            if (angular_distance(pt, keep) < min_sep) {
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

// (best-candidate selection removed by request)

void write_healpix_deterministic_samples(int nside, double angular_step_rad, int sampling_method) {
    long npix = 12 * nside * nside;
    Healpix_Base healpix_base(nside, RING, SET_NSIDE);
    double pixel_size = std::sqrt(4 * M_PI / npix);  // Approximate size of each pixel in radians
    
    std::cout << "Pixel size (radians): " << pixel_size << std::endl;
    std::cout << "Pixel size (degrees): " << pixel_size * (180.0 / M_PI) << std::endl;
    
    // Create output directory
    std::string dir_name = "nside_" + std::to_string(nside);
    
    // Cross-platform directory creation
    #ifdef _WIN32
        _mkdir(dir_name.c_str());
    #else
        mkdir(dir_name.c_str(), 0755);
    #endif
    
    // Estimate candidate pool size from hexagonal packing with minimum separation angular_step_rad
    // Area per point ≈ (sqrt(3)/2) * d^2 for minimum distance d
    const double area_per_point = (std::sqrt(3.0) / 2.0) * angular_step_rad * angular_step_rad;
    // Use higher multiplier to reach near-saturation after greedy min-sep filtering
    const double candidate_density_multiplier = 6.0; // empirically good; deterministic
    int estimated_points = static_cast<int>(std::ceil(((4.0 * M_PI) / area_per_point) * candidate_density_multiplier));
    
    // Auto-densify loop: increase candidate density until per-pixel reaches target fraction of hex-pack
    const double pixel_area = (4.0 * M_PI) / static_cast<double>(npix);
    const double expected_per_pixel = pixel_area / area_per_point;
    const double target_fraction = 0.90; // 90% of hex-pack
    const int target_min_per_pixel = std::max(1, static_cast<int>(std::floor(expected_per_pixel * target_fraction)));
    std::cout << "Expected (hex-pack) points per pixel ≈ " << std::llround(expected_per_pixel)
              << ", target minimum per pixel = " << target_min_per_pixel << std::endl;

    std::vector<std::vector<pointing>> best_samples_per_pixel(npix);
    int best_min_count = -1;

    const int max_attempts = 6;
    double density_multiplier = 2.0; // start modestly
    for (int attempt = 0; attempt < max_attempts; ++attempt) {
        const int num_candidates = static_cast<int>(std::ceil(((4.0 * M_PI) / area_per_point) * density_multiplier));
        std::vector<pointing> all_points;
        all_points.reserve(num_candidates);
        if (sampling_method == 1) {
            std::cout << "Attempt " << (attempt + 1) << "/" << max_attempts << ": Fibonacci candidates = " << num_candidates << std::endl;
            all_points = generate_fibonacci_points(num_candidates);
        } else {
            std::cout << "Attempt " << (attempt + 1) << "/" << max_attempts << ": Hopf candidates = " << num_candidates << std::endl;
            all_points = generate_hopf_points(num_candidates);
        }

        // Deterministic shuffle per attempt to reduce ordering bias
        {
            const uint64_t seed = static_cast<uint64_t>(nside) * 1315423911ULL ^
                                  static_cast<uint64_t>(sampling_method) * 2654435761ULL ^
                                  static_cast<uint64_t>(std::llround(angular_step_rad * 1e9)) ^
                                  static_cast<uint64_t>(attempt);
            std::mt19937_64 rng(seed);
            std::shuffle(all_points.begin(), all_points.end(), rng);
        }

        // Per-pixel greedy min-separation selection
        std::vector<std::vector<pointing>> samples_per_pixel(npix);
        samples_per_pixel.shrink_to_fit();
        int min_count = std::numeric_limits<int>::max();
        long total_count = 0;
        for (long pix = 0; pix < npix; ++pix) {
            auto samples = select_points_for_pixel_with_min_spacing(all_points, healpix_base, pix, angular_step_rad);
            if (pix == 0 && attempt == 0) {
                std::cout << "First pixel (attempt 1) generated " << samples.size() << " points" << std::endl;
            }
            min_count = std::min(min_count, static_cast<int>(samples.size()));
            total_count += static_cast<long>(samples.size());
            samples_per_pixel[pix] = std::move(samples);
        }

        const double avg_count = static_cast<double>(total_count) / static_cast<double>(npix);
        std::cout << "Attempt " << (attempt + 1) << ": min/avg per pixel = " << min_count << "/" << avg_count << std::endl;

        if (min_count > best_min_count) {
            best_min_count = min_count;
            best_samples_per_pixel = std::move(samples_per_pixel);
        }

        if (min_count >= target_min_per_pixel) {
            std::cout << "Target reached. Proceeding with current selection." << std::endl;
            break;
        }

        density_multiplier *= 1.8; // densify
    }

    // Write out using best achieved samples
    for (long pix = 0; pix < npix; ++pix) {
        const auto& samples = best_samples_per_pixel[pix];
        std::stringstream filename;
        filename << dir_name << "/healpix_data_" << pix << ".txt";
        std::ofstream file(filename.str());
        if (!file.is_open()) {
            std::cerr << "Failed to open file: " << filename.str() << std::endl;
            continue;
        }

        //file << "Phi[Degrees] Theta[Degrees]" << std::endl;
        for (const auto& sample : samples) {
            double theta_deg = sample.theta * (180.0 / M_PI);
            double phi_deg = sample.phi * (180.0 / M_PI);
            theta_deg = round(theta_deg * 10000.0) / 10000.0;
            phi_deg = round(phi_deg * 10000.0) / 10000.0;
            file << phi_deg << " " << theta_deg << std::endl;
        }
        file.close();

        if ((pix + 1) % 100 == 0) {
            std::cout << "Processed " << (pix + 1) << "/" << npix << " pixels" << std::endl;
        }
    }
    
    // Calculate and display final statistics
    long total_points = 0;
    for (long pix = 0; pix < npix; ++pix) {
        std::string filename = dir_name + "/healpix_data_" + std::to_string(pix) + ".txt";
        std::ifstream file(filename);
        if (file.is_open()) {
            std::string line;
            int line_count = 0;
            while (std::getline(file, line)) {
                line_count++;
            }
            total_points += (line_count - 1); // Subtract 1 for header
            file.close();
        }
    }
    
    std::cout << "\n=== SAMPLING STATISTICS ===" << std::endl;
    std::cout << "nside: " << nside << std::endl;
    std::cout << "Total pixels: " << npix << std::endl;
    std::cout << "Angular step size: " << (angular_step_rad * 180.0 / M_PI) << "°" << std::endl;
    std::cout << "Total sampling points: " << total_points << std::endl;
    std::cout << "Average points per pixel: " << (double)total_points / npix << std::endl;
    std::cout << "===========================" << std::endl;
}

//int main() {
//    int nside = 1;  // Example nside
//    write_healpix_random_samples(nside);
//    return 0;
//}

int main() {
    int nside;
    float ang_step_size;
    
    std::cout << "Enter the value of nside (e.g., 1, 2, 3, 4, 5, 6, 8, 9, 12, 16, etc.): ";
    std::cin >> nside;
    
    // Validate nside is positive
    if (nside <= 0) {
        std::cerr << "Error: nside must be a positive integer" << std::endl;
        return 1;
    }
    
    std::cout << "Enter the angular step size in degrees (e.g., 5.7): ";
    std::cin >> ang_step_size;
    
    std::cout << "Choose sampling method:" << std::endl;
    std::cout << "1. Fibonacci sphere sampling" << std::endl;
    std::cout << "2. Hopf fibration sampling" << std::endl;
    std::cout << "Enter choice (1 or 2): ";
    
    int sampling_method;
    std::cin >> sampling_method;
    
    if (sampling_method != 1 && sampling_method != 2) {
        std::cerr << "Error: Invalid choice. Using Fibonacci sampling." << std::endl;
        sampling_method = 1;
    }
    
    std::cout << "The number of HEALPix pixels created: " << 12 * nside * nside << std::endl;
    std::cout << "Note: Using " << ((nside & (nside - 1)) == 0 ? "hierarchical" : "extended") << " HEALPix scheme" << std::endl;
    std::cout << "Generating deterministic samples using " << (sampling_method == 1 ? "Fibonacci" : "Hopf") << " method..." << std::endl;
    
    // Convert angular step size from degrees to radians for internal calculations
    double ang_step_rad = ang_step_size * (M_PI / 180.0);
    std::cout << "Angular step size: " << ang_step_size << "° (" << ang_step_rad << " radians)" << std::endl;
    std::cout << "Note: Internal calculations use radians, output files show degrees" << std::endl;
    std::cout << "Number of points per pixel will be determined automatically based on spacing..." << std::endl;
    
    write_healpix_deterministic_samples(nside, ang_step_rad, sampling_method);
    
    std::cout << "Sampling complete! Files saved in nside_" << nside << " folder." << std::endl;
    return 0;
}
