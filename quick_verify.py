#!/usr/bin/env python3
"""
Quick verification script to check if HEALPix sampling points are properly bounded.
This is a lightweight version for fast verification.
"""

import numpy as np
import matplotlib.pyplot as plt
import healpy as hp
import os
import glob

def quick_verify_and_plot(nside, data_dir=None, show_plot=True):
    """
    Quick verification of HEALPix pixel containment with simple visualization.
    
    Parameters:
    -----------
    nside : int
        HEALPix nside parameter
    data_dir : str, optional
        Directory containing the data files
    show_plot : bool
        Whether to show the plot
    
    Returns:
    --------
    dict : Verification results
    """
    if data_dir is None:
        data_dir = f"nside_{nside}"
    
    if not os.path.exists(data_dir):
        print(f"Error: Directory {data_dir} not found. Run the C++ code first.")
        return None
    
    print(f"Verifying HEALPix sampling for nside={nside}")
    print("-" * 50)
    
    # Load and verify all data
    pattern = os.path.join(data_dir, "healpix_data_*.txt")
    files = glob.glob(pattern)
    
    total_points = 0
    correct_assignments = 0
    wrong_assignments = []
    pixel_counts = {}
    all_points = []
    
    for file_path in files:
        # Extract pixel ID
        filename = os.path.basename(file_path)
        assigned_pixel = int(filename.replace("healpix_data_", "").replace(".txt", ""))
        
        try:
            # Load data
            data = np.loadtxt(file_path)
            if data.size == 0:
                pixel_counts[assigned_pixel] = 0
                continue
                
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            # Convert degrees to radians
            phi_deg, theta_deg = data[:, 0], data[:, 1]
            phi_rad = np.deg2rad(phi_deg)
            theta_rad = np.deg2rad(theta_deg)
            
            pixel_counts[assigned_pixel] = len(phi_rad)
            
            # Verify each point
            for phi, theta in zip(phi_rad, theta_rad):
                computed_pixel = hp.ang2pix(nside, theta, phi)
                total_points += 1
                
                if computed_pixel == assigned_pixel:
                    correct_assignments += 1
                else:
                    wrong_assignments.append((assigned_pixel, computed_pixel, phi, theta))
                
                all_points.append((phi, theta, assigned_pixel, computed_pixel))
                
        except Exception as e:
            print(f"Warning: Could not process {file_path}: {e}")
    
    # Calculate statistics
    npix = 12 * nside * nside
    accuracy = correct_assignments / total_points if total_points > 0 else 0
    
    print(f"Total pixels: {npix}")
    print(f"Pixels with data: {len(pixel_counts)}")
    print(f"Total sample points: {total_points}")
    print(f"Correct assignments: {correct_assignments}")
    print(f"Wrong assignments: {len(wrong_assignments)}")
    print(f"Accuracy: {accuracy:.4f} ({accuracy*100:.2f}%)")
    
    # Points per pixel statistics
    if pixel_counts:
        counts = list(pixel_counts.values())
        print(f"\nPoints per pixel statistics:")
        print(f"  Mean: {np.mean(counts):.2f}")
        print(f"  Min: {np.min(counts)}")
        print(f"  Max: {np.max(counts)}")
        print(f"  Std: {np.std(counts):.2f}")
    
    # Show first few wrong assignments
    if wrong_assignments:
        print(f"\nFirst 5 misassigned points:")
        for i, (assigned, computed, phi, theta) in enumerate(wrong_assignments[:5]):
            phi_deg = np.rad2deg(phi)
            theta_deg = np.rad2deg(theta)
            print(f"  {i+1}: Assigned to pixel {assigned}, computed as pixel {computed}")
            print(f"      Position: φ={phi_deg:.3f}°, θ={theta_deg:.3f}°")
    
    # Create simple visualization if requested
    if show_plot and all_points:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
        
        # Plot 1: All points colored by assigned pixel
        all_phi, all_theta, assigned_pix, computed_pix = zip(*all_points)
        all_phi = np.array(all_phi)
        all_theta = np.array(all_theta)
        assigned_pix = np.array(assigned_pix)
        computed_pix = np.array(computed_pix)
        
        # Convert to lat/lon for plotting
        lon = all_phi - np.pi
        lat = np.pi/2 - all_theta
        
        scatter1 = ax1.scatter(np.rad2deg(lon), np.rad2deg(lat), c=assigned_pix, 
                             s=10, alpha=0.7, cmap='tab20')
        ax1.set_title(f'Points colored by assigned pixel (nside={nside})')
        ax1.set_xlabel('Longitude (degrees)')
        ax1.set_ylabel('Latitude (degrees)')
        ax1.grid(True, alpha=0.3)
        plt.colorbar(scatter1, ax=ax1, label='Pixel ID')
        
        # Plot 2: Show misassigned points
        if wrong_assignments:
            # Correct points in blue
            correct_mask = assigned_pix == computed_pix
            ax2.scatter(np.rad2deg(lon[correct_mask]), np.rad2deg(lat[correct_mask]), 
                       c='blue', s=5, alpha=0.5, label='Correct')
            
            # Wrong points in red
            wrong_mask = ~correct_mask
            ax2.scatter(np.rad2deg(lon[wrong_mask]), np.rad2deg(lat[wrong_mask]), 
                       c='red', s=20, alpha=0.8, label='Misassigned')
            
            ax2.set_title(f'Correct vs Misassigned Points\n({len(wrong_assignments)} errors)')
        else:
            ax2.scatter(np.rad2deg(lon), np.rad2deg(lat), c='green', s=5, alpha=0.7)
            ax2.set_title('All Points Correctly Assigned!')
        
        ax2.set_xlabel('Longitude (degrees)')
        ax2.set_ylabel('Latitude (degrees)')
        ax2.grid(True, alpha=0.3)
        ax2.legend()
        
        plt.tight_layout()
        plt.show()
    
    # Return results
    return {
        'nside': nside,
        'total_pixels': npix,
        'pixels_with_data': len(pixel_counts),
        'total_points': total_points,
        'correct_assignments': correct_assignments,
        'wrong_assignments': len(wrong_assignments),
        'accuracy': accuracy,
        'pixel_counts': pixel_counts,
        'wrong_details': wrong_assignments[:10]
    }

def main():
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: python quick_verify.py <nside> [data_dir]")
        print("Example: python quick_verify.py 1")
        return
    
    nside = int(sys.argv[1])
    data_dir = sys.argv[2] if len(sys.argv) > 2 else None
    
    result = quick_verify_and_plot(nside, data_dir, show_plot=True)
    
    if result and result['accuracy'] >= 0.99:
        print("\n✅ VERIFICATION PASSED: Points are properly bounded within HEALPix pixels!")
    elif result:
        print(f"\n⚠️  VERIFICATION ISSUES: Only {result['accuracy']*100:.2f}% accuracy")
    else:
        print("\n❌ VERIFICATION FAILED: Could not load data")

if __name__ == "__main__":
    main()