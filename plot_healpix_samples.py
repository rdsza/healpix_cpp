#!/usr/bin/env python3
"""
HEALPix Sampling Visualization Script

This script reads the HEALPix sampling data files and plots them on a sphere (S2)
using:
- 3D scatter with optional HEALPix pixel boundary overlay (recommended for figures)
- Color-coding by pixel index

Requirements:
- matplotlib
- numpy
- healpy (optional, needed for pixel boundary overlay)
- glob (built-in)
- os (built-in)
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from matplotlib.colors import Normalize
import matplotlib.patches as mpatches
try:
    import healpy as hp
    HAVE_HEALPY = True
except Exception:
    HAVE_HEALPY = False

def read_healpix_data(nside):
    """
    Read all HEALPix data files for a given nside.
    
    Args:
        nside (int): The nside parameter
        
    Returns:
        tuple: (theta_deg, phi_deg, pixel_indices) where all are numpy arrays
    """
    folder_name = f"nside_{nside}"
    if not os.path.exists(folder_name):
        print(f"Error: Folder {folder_name} not found!")
        return None, None, None
    
    all_theta = []
    all_phi = []
    pixel_indices = []
    
    # Find all data files
    pattern = os.path.join(folder_name, "healpix_data_*.txt")
    files = glob.glob(pattern)
    
    if not files:
        print(f"Error: No data files found in {folder_name}")
        return None, None, None
    
    print(f"Found {len(files)} data files in {folder_name}")
    
    for file_path in sorted(files):
        # Extract pixel index from filename
        filename = os.path.basename(file_path)
        pixel_idx = int(filename.split('_')[-1].split('.')[0])
        
        try:
            data = np.loadtxt(file_path, skiprows=1)  # Skip header
            if len(data.shape) == 2 and data.shape[1] == 2:
                # Files now store columns as: Phi[Degrees], Theta[Degrees]
                phi_deg = data[:, 0]
                theta_deg = data[:, 1]
                
                all_theta.extend(theta_deg)
                all_phi.extend(phi_deg)
                pixel_indices.extend([pixel_idx] * len(theta_deg))
                
        except Exception as e:
            print(f"Warning: Could not read {file_path}: {e}")
            continue
    
    if not all_theta:
        print("Error: No valid data found!")
        return None, None, None
    
    return np.array(all_theta), np.array(all_phi), np.array(pixel_indices)

def plot_healpix_samples(nside, save_plot=True):
    """
    Plot HEALPix sampling points on a sphere.
    
    Args:
        nside (int): The nside parameter
        save_plot (bool): Whether to save the plot
    """
    # Read data
    theta_deg, phi_deg, pixel_indices = read_healpix_data(nside)
    
    if theta_deg is None:
        return
    
    # Convert to radians for calculations
    theta_rad = np.radians(theta_deg)
    phi_rad = np.radians(phi_deg)
    
    # Create figure with spherical projection
    fig = plt.figure(figsize=(12, 8))
    
    # Create 3D projection
    ax = fig.add_subplot(111, projection='3d')
    
    # Convert spherical to Cartesian coordinates
    x = np.sin(theta_rad) * np.cos(phi_rad)
    y = np.sin(theta_rad) * np.sin(phi_rad)
    z = np.cos(theta_rad)
    
    # Create colormap based on pixel indices
    norm = Normalize(vmin=pixel_indices.min(), vmax=pixel_indices.max())
    colors = plt.cm.jet(norm(pixel_indices))
    
    # Plot points
    scatter = ax.scatter(x, y, z, c=pixel_indices, cmap='jet', 
                         s=20, alpha=0.7, edgecolors='none')
    
    # Add colorbar
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20)
    cbar.set_label('Pixel Index', fontsize=12)
    
    # Set labels and title
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.set_title(f'HEALPix Sampling (nside={nside})\n{len(theta_deg)} total points, {len(set(pixel_indices))} pixels', 
                 fontsize=14, fontweight='bold')
    
    # Set equal aspect ratio
    ax.set_box_aspect([1, 1, 1])
    
    # Set view angle for better visualization
    ax.view_init(elev=20, azim=45)
    
    # Add grid
    ax.grid(True, alpha=0.3)
    
    # Add sphere outline
    u = np.linspace(0, 2 * np.pi, 100)
    v = np.linspace(0, np.pi, 100)
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones(np.size(u)), np.cos(v))
    
    ax.plot_wireframe(sphere_x, sphere_y, sphere_z, color='gray', alpha=0.2, linewidth=0.8)
    
    plt.tight_layout()
    
    if save_plot:
        plot_filename = f"healpix_visualization_nside_{nside}.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {plot_filename}")
    
    plt.show()

def plot_healpix_samples_with_boundaries(nside, save_plot=True, fill_faces=False, edge_alpha=0.5):
    """
    Plot points and overlay HEALPix pixel boundaries on a 3D sphere.

    Args:
        nside (int): HEALPix nside
        save_plot (bool): Save plot to PNG
        fill_faces (bool): If True, lightly fill pixel faces (slower)
        edge_alpha (float): Edge transparency for boundaries
    """
    if not HAVE_HEALPY:
        print("healpy is required for boundary overlay. Install with: pip install healpy")
        return

    theta_deg, phi_deg, pixel_indices = read_healpix_data(nside)
    if theta_deg is None:
        return

    theta_rad = np.radians(theta_deg)
    phi_rad = np.radians(phi_deg)

    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')

    # Points
    x = np.sin(theta_rad) * np.cos(phi_rad)
    y = np.sin(theta_rad) * np.sin(phi_rad)
    z = np.cos(theta_rad)
    norm = Normalize(vmin=pixel_indices.min(), vmax=pixel_indices.max())
    scatter = ax.scatter(x, y, z, c=pixel_indices, cmap='jet', s=10, alpha=0.8, edgecolors='none')
    cbar = plt.colorbar(scatter, ax=ax, shrink=0.6, aspect=20)
    cbar.set_label('Pixel Index', fontsize=12)

    # Pixel boundaries
    from mpl_toolkits.mplot3d.art3d import Poly3DCollection
    npix = 12 * nside * nside
    for ipix in range(npix):
        # hp.boundaries returns (3, M) array of unit vectors outlining the pixel
        verts = hp.boundaries(nside, ipix, step=1)  # (3, M)
        bx, by, bz = verts[0], verts[1], verts[2]
        # Draw boundary polyline
        ax.plot(bx, by, bz, color='k', linewidth=1.5, alpha=edge_alpha)
        if fill_faces:
            poly = Poly3DCollection([list(zip(bx, by, bz))], facecolors=(0.8, 0.8, 0.8, 0.05), edgecolors='none')
            ax.add_collection3d(poly)

    # Labels and cosmetics
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_zlabel('Z', fontsize=12)
    ax.set_title(f'HEALPix Sampling with Boundaries (nside={nside})\n{len(theta_deg)} total points, {int(npix)} pixels',
                 fontsize=14, fontweight='bold')
    ax.set_box_aspect([1, 1, 1])
    ax.view_init(elev=20, azim=45)
    ax.grid(True, alpha=0.2)

    # Sphere wireframe
    u = np.linspace(0, 2 * np.pi, 80)
    v = np.linspace(0, np.pi, 40)
    sphere_x = np.outer(np.cos(u), np.sin(v))
    sphere_y = np.outer(np.sin(u), np.sin(v))
    sphere_z = np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(sphere_x, sphere_y, sphere_z, color='gray', alpha=0.15, linewidth=0.6)

    plt.tight_layout()
    if save_plot:
        plot_filename = f"healpix_visualization_boundaries_nside_{nside}.png"
        plt.savefig(plot_filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved as {plot_filename}")
    plt.show()

def plot_multiple_nside(nside_values, save_plots=True):
    """
    Plot multiple nside values for comparison.
    
    Args:
        nside_values (list): List of nside values to plot
        save_plots (bool): Whether to save the plots
    """
    for nside in nside_values:
        print(f"\n{'='*50}")
        print(f"Processing nside = {nside}")
        print(f"{'='*50}")
        
        try:
            plot_healpix_samples(nside, save_plots)
        except Exception as e:
            print(f"Error plotting nside {nside}: {e}")
            continue

def main():
    """Main function to run the visualization."""
    print("HEALPix Sampling Visualization Tool")
    print("="*40)
    
    # Check if any nside folders exist
    nside_folders = [d for d in os.listdir('.') if d.startswith('nside_') and os.path.isdir(d)]
    
    if not nside_folders:
        print("No nside folders found in current directory.")
        print("Please run the C++ program first to generate data files.")
        return
    
    # Extract nside values
    available_nsides = []
    for folder in nside_folders:
        try:
            nside = int(folder.split('_')[1])
            available_nsides.append(nside)
        except:
            continue
    
    available_nsides.sort()
    print(f"Available nside values: {available_nsides}")
    
    if len(available_nsides) == 1:
        print(f"Plotting nside = {available_nsides[0]}")
        plot_healpix_samples(available_nsides[0])
    else:
        print("\nChoose visualization option:")
        print("1. Plot single nside (scatter only)")
        print("2. Plot all available nsides (scatter only)")
        print("3. Plot single nside with HEALPix pixel boundaries (3D)")
        
        try:
            choice = input("Enter choice (1 or 2): ").strip()
            
            if choice == "1":
                nside = int(input(f"Enter nside value ({', '.join(map(str, available_nsides))}): "))
                if nside in available_nsides:
                    plot_healpix_samples(nside)
                else:
                    print(f"Invalid nside. Available values: {available_nsides}")
            elif choice == "2":
                plot_multiple_nside(available_nsides)
            elif choice == "3":
                nside = int(input(f"Enter nside value ({', '.join(map(str, available_nsides))}): "))
                if nside in available_nsides:
                    plot_healpix_samples_with_boundaries(nside)
                else:
                    print(f"Invalid nside. Available values: {available_nsides}")
            else:
                print("Invalid choice. Plotting first available nside.")
                plot_healpix_samples(available_nsides[0])
                
        except (ValueError, KeyboardInterrupt):
            print("\nPlotting first available nside.")
            plot_healpix_samples(available_nsides[0])

if __name__ == "__main__":
    main()
