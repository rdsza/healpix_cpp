#!/usr/bin/env python3
"""
Visualize HEALPix sampling output from the C++ code to verify proper pixel boundaries.
This script reads the output files and plots both the sampled points and HEALPix pixel boundaries on S2.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import healpy as hp
import os
import glob
import argparse
from mpl_toolkits.mplot3d import Axes3D

def load_healpix_samples(nside, data_dir=None):
    """
    Load all HEALPix sample files for a given nside.
    
    Parameters:
    -----------
    nside : int
        The HEALPix nside parameter
    data_dir : str, optional
        Directory containing the data files. If None, uses 'nside_<nside>'
    
    Returns:
    --------
    dict : Dictionary with pixel_id as key and list of (phi, theta) tuples as values
    """
    if data_dir is None:
        data_dir = f"nside_{nside}"
    
    if not os.path.exists(data_dir):
        raise FileNotFoundError(f"Directory {data_dir} not found. Run the C++ code first.")
    
    pixel_data = {}
    pattern = os.path.join(data_dir, "healpix_data_*.txt")
    files = glob.glob(pattern)
    
    print(f"Found {len(files)} data files in {data_dir}")
    
    for file_path in files:
        # Extract pixel ID from filename
        filename = os.path.basename(file_path)
        pixel_id = int(filename.replace("healpix_data_", "").replace(".txt", ""))
        
        # Load the data
        try:
            data = np.loadtxt(file_path)
            if data.size == 0:
                continue
                
            if data.ndim == 1:
                data = data.reshape(1, -1)
            
            # Convert from degrees to radians for internal calculations
            phi_deg, theta_deg = data[:, 0], data[:, 1]
            phi_rad = np.deg2rad(phi_deg)
            theta_rad = np.deg2rad(theta_deg)
            
            pixel_data[pixel_id] = list(zip(phi_rad, theta_rad))
        except Exception as e:
            print(f"Warning: Could not load {file_path}: {e}")
    
    return pixel_data

def get_healpix_pixel_boundaries(nside, pixel_id, resolution=100):
    """
    Get the boundary of a HEALPix pixel.
    
    Parameters:
    -----------
    nside : int
        HEALPix nside parameter
    pixel_id : int
        Pixel index
    resolution : int
        Number of points along the boundary
    
    Returns:
    --------
    tuple : (phi_boundary, theta_boundary) arrays in radians
    """
    # Get pixel boundaries using healpy
    boundary = hp.boundaries(nside, pixel_id, step=resolution)
    
    # Convert from vectors to spherical coordinates
    theta_boundary = np.arccos(boundary[2])
    phi_boundary = np.arctan2(boundary[1], boundary[0])
    phi_boundary[phi_boundary < 0] += 2 * np.pi
    
    return phi_boundary, theta_boundary

def plot_mollweide_projection(pixel_data, nside, show_boundaries=True, 
                             colors_per_pixel=True, figsize=(15, 8)):
    """
    Plot the HEALPix sampling on a Mollweide projection.
    
    Parameters:
    -----------
    pixel_data : dict
        Dictionary with pixel data
    nside : int
        HEALPix nside parameter
    show_boundaries : bool
        Whether to show pixel boundaries
    colors_per_pixel : bool
        Whether to use different colors for each pixel
    figsize : tuple
        Figure size
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, 
                                   subplot_kw={'projection': 'mollweide'})
    
    # Generate colormap for pixels
    npix = 12 * nside * nside
    colors = plt.cm.tab20(np.linspace(0, 1, min(20, npix)))
    if npix > 20:
        colors = plt.cm.viridis(np.linspace(0, 1, npix))
    
    # Plot 1: Sample points only
    ax1.set_title(f'HEALPix Sampling Points (nside={nside})')
    
    total_points = 0
    for pixel_id, points in pixel_data.items():
        if not points:
            continue
            
        phi_vals, theta_vals = zip(*points)
        phi_vals = np.array(phi_vals)
        theta_vals = np.array(theta_vals)
        
        # Convert to longitude/latitude for Mollweide
        lon = phi_vals - np.pi  # Convert to [-π, π]
        lat = np.pi/2 - theta_vals  # Convert to [-π/2, π/2]
        
        color = colors[pixel_id % len(colors)] if colors_per_pixel else 'blue'
        ax1.scatter(lon, lat, c=[color], s=10, alpha=0.7, 
                   label=f'Pixel {pixel_id}' if pixel_id < 5 else "")
        total_points += len(points)
    
    ax1.grid(True, alpha=0.3)
    if npix <= 5:
        ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Plot 2: Sample points with pixel boundaries
    ax2.set_title(f'HEALPix Sampling with Pixel Boundaries (nside={nside})')
    
    # Plot boundaries first (so they appear behind points)
    if show_boundaries:
        print("Computing pixel boundaries...")
        for pixel_id in range(npix):
            try:
                phi_bound, theta_bound = get_healpix_pixel_boundaries(nside, pixel_id, resolution=50)
                
                # Convert to longitude/latitude for Mollweide
                lon_bound = phi_bound - np.pi
                lat_bound = np.pi/2 - theta_bound
                
                # Handle wrap-around at boundaries
                if np.max(lon_bound) - np.min(lon_bound) > np.pi:
                    # Split the boundary at the wrap-around point
                    wrap_idx = np.where(np.abs(np.diff(lon_bound)) > np.pi)[0]
                    if len(wrap_idx) > 0:
                        for i in range(len(wrap_idx) + 1):
                            start = wrap_idx[i-1] + 1 if i > 0 else 0
                            end = wrap_idx[i] + 1 if i < len(wrap_idx) else len(lon_bound)
                            ax2.plot(lon_bound[start:end], lat_bound[start:end], 
                                   'k-', linewidth=0.5, alpha=0.6)
                    else:
                        ax2.plot(lon_bound, lat_bound, 'k-', linewidth=0.5, alpha=0.6)
                else:
                    ax2.plot(lon_bound, lat_bound, 'k-', linewidth=0.5, alpha=0.6)
                    
            except Exception as e:
                print(f"Warning: Could not compute boundary for pixel {pixel_id}: {e}")
    
    # Plot sample points on top
    for pixel_id, points in pixel_data.items():
        if not points:
            continue
            
        phi_vals, theta_vals = zip(*points)
        phi_vals = np.array(phi_vals)
        theta_vals = np.array(theta_vals)
        
        # Convert to longitude/latitude for Mollweide
        lon = phi_vals - np.pi
        lat = np.pi/2 - theta_vals
        
        color = colors[pixel_id % len(colors)] if colors_per_pixel else 'red'
        ax2.scatter(lon, lat, c=[color], s=15, alpha=0.8, edgecolors='white', linewidth=0.5)
    
    ax2.grid(True, alpha=0.3)
    
    # Add statistics
    fig.suptitle(f'Total Points: {total_points}, Pixels: {len(pixel_data)}/{npix}', 
                 fontsize=14, y=0.95)
    
    plt.tight_layout()
    return fig

def plot_3d_sphere(pixel_data, nside, show_boundaries=False, figsize=(12, 8)):
    """
    Plot the HEALPix sampling on a 3D sphere.
    
    Parameters:
    -----------
    pixel_data : dict
        Dictionary with pixel data
    nside : int
        HEALPix nside parameter
    show_boundaries : bool
        Whether to show pixel boundaries (can be slow)
    figsize : tuple
        Figure size
    """
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection='3d')
    
    # Generate colormap
    npix = 12 * nside * nside
    colors = plt.cm.tab20(np.linspace(0, 1, min(20, npix)))
    if npix > 20:
        colors = plt.cm.viridis(np.linspace(0, 1, npix))
    
    # Plot sample points
    total_points = 0
    for pixel_id, points in pixel_data.items():
        if not points:
            continue
            
        phi_vals, theta_vals = zip(*points)
        phi_vals = np.array(phi_vals)
        theta_vals = np.array(theta_vals)
        
        # Convert to Cartesian coordinates
        x = np.sin(theta_vals) * np.cos(phi_vals)
        y = np.sin(theta_vals) * np.sin(phi_vals)
        z = np.cos(theta_vals)
        
        color = colors[pixel_id % len(colors)]
        ax.scatter(x, y, z, c=[color], s=20, alpha=0.8)
        total_points += len(points)
    
    # Draw wireframe sphere for reference
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    x_sphere = np.outer(np.sin(v), np.cos(u))
    y_sphere = np.outer(np.sin(v), np.sin(u))
    z_sphere = np.outer(np.cos(v), np.ones(np.size(u)))
    ax.plot_wireframe(x_sphere, y_sphere, z_sphere, alpha=0.1, color='gray')
    
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title(f'3D Sphere View (nside={nside}, {total_points} points)')
    
    # Set equal aspect ratio
    ax.set_box_aspect([1,1,1])
    
    return fig

def verify_pixel_containment(pixel_data, nside):
    """
    Verify that all points are actually contained within their assigned pixels.
    
    Parameters:
    -----------
    pixel_data : dict
        Dictionary with pixel data
    nside : int
        HEALPix nside parameter
    
    Returns:
    --------
    dict : Statistics about containment verification
    """
    print("\nVerifying pixel containment...")
    
    total_points = 0
    correct_assignments = 0
    wrong_assignments = []
    
    for pixel_id, points in pixel_data.items():
        for phi, theta in points:
            # Convert to HEALPix pixel
            computed_pixel = hp.ang2pix(nside, theta, phi)
            
            total_points += 1
            if computed_pixel == pixel_id:
                correct_assignments += 1
            else:
                wrong_assignments.append((pixel_id, computed_pixel, phi, theta))
    
    stats = {
        'total_points': total_points,
        'correct_assignments': correct_assignments,
        'wrong_assignments': len(wrong_assignments),
        'accuracy': correct_assignments / total_points if total_points > 0 else 0,
        'wrong_details': wrong_assignments[:10]  # First 10 wrong assignments
    }
    
    print(f"Total points checked: {total_points}")
    print(f"Correct assignments: {correct_assignments}")
    print(f"Wrong assignments: {len(wrong_assignments)}")
    print(f"Accuracy: {stats['accuracy']:.2%}")
    
    if wrong_assignments:
        print("\nFirst few wrong assignments:")
        for assigned_pix, computed_pix, phi, theta in wrong_assignments[:5]:
            print(f"  Point assigned to pixel {assigned_pix}, but belongs to pixel {computed_pix}")
    
    return stats

def main():
    parser = argparse.ArgumentParser(description='Visualize HEALPix sampling output')
    parser.add_argument('nside', type=int, help='HEALPix nside parameter')
    parser.add_argument('--data-dir', type=str, default=None, 
                       help='Directory containing data files (default: nside_<nside>)')
    parser.add_argument('--no-boundaries', action='store_true', 
                       help='Do not show pixel boundaries (faster)')
    parser.add_argument('--3d', action='store_true', 
                       help='Show 3D sphere plot')
    parser.add_argument('--verify-only', action='store_true', 
                       help='Only verify containment, do not plot')
    parser.add_argument('--save', type=str, default=None, 
                       help='Save plots to file (without extension)')
    
    args = parser.parse_args()
    
    try:
        # Load the data
        print(f"Loading HEALPix samples for nside={args.nside}...")
        pixel_data = load_healpix_samples(args.nside, args.data_dir)
        
        if not pixel_data:
            print("No data found! Make sure you've run the C++ code first.")
            return
        
        # Verify containment
        stats = verify_pixel_containment(pixel_data, args.nside)
        
        if args.verify_only:
            return
        
        # Create plots
        if not args.no_boundaries or args.nside <= 4:
            print("Creating Mollweide projection plot...")
            fig1 = plot_mollweide_projection(pixel_data, args.nside, 
                                           show_boundaries=not args.no_boundaries)
            if args.save:
                fig1.savefig(f"{args.save}_mollweide.png", dpi=300, bbox_inches='tight')
                print(f"Saved Mollweide plot to {args.save}_mollweide.png")
        
        if args._3d:
            print("Creating 3D sphere plot...")
            fig2 = plot_3d_sphere(pixel_data, args.nside, 
                                show_boundaries=False)  # Boundaries too slow in 3D
            if args.save:
                fig2.savefig(f"{args.save}_3d.png", dpi=300, bbox_inches='tight')
                print(f"Saved 3D plot to {args.save}_3d.png")
        
        if not args.save:
            plt.show()
        
        print(f"\nVisualization complete!")
        print(f"Containment accuracy: {stats['accuracy']:.2%}")
        
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    main()