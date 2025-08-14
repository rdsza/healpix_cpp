# HEALPix Deterministic Sampling System

This system generates deterministic sampling points within HEALPix pixels using either Fibonacci sphere sampling or Hopf fibration, ensuring equal representation across all pixels for statistical analysis.

## Features

- **Two Sampling Algorithms**: Choose between Fibonacci sphere sampling and Hopf fibration
- **Deterministic Results**: Same parameters always produce the same sampling pattern
- **Equal Pixel Representation**: Each pixel gets the same number of sampling points
- **Flexible nside Support**: Works with any positive integer nside value
- **Organized Output**: Files saved in organized folders with consistent naming
- **Visualization**: Python script for 3D visualization of sampling patterns

## C++ Program (main2.cpp)

### Compilation
```bash
g++ -std=c++17 -O2 -o healpix_sampler main2.cpp -I./src/Healpix_2.15a
```

### Usage
```bash
./healpix_sampler
```

The program will prompt you for:
1. **nside**: HEALPix resolution parameter (e.g., 1, 2, 3, 4, 5, 6, 8, 9, 12, 16)
2. **Angular step size**: Maximum angular distance from pixel center in degrees (e.g., 5.7)
3. **Sampling method**: 1 for Fibonacci, 2 for Hopf fibration

### Output
- Creates a folder `nside_X/` for each nside value
- Each pixel gets its own file: `healpix_data_PIXEL_INDEX.txt`
- Files contain two columns: Theta[Degrees] Phi[Degrees]
- Angles are rounded to 4 significant digits
- **Automatic point calculation**: Number of points per pixel is determined by angular spacing
- **Consistent results**: Same nside + spacing always produces the same number of points

## Python Visualization (plot_healpix_samples.py)

### Installation
```bash
pip install -r requirements.txt
```

### Usage
```bash
python plot_healpix_samples.py
```

The script will:
1. Automatically detect available nside folders
2. Offer to plot single or multiple nside values
3. Create 3D spherical visualizations with colormaps
4. Save high-resolution plots as PNG files

### Visualization Features
- **3D Sphere Projection**: Points plotted on a unit sphere
- **Color-coded by Pixel**: Each pixel gets a unique color
- **Interactive Viewing**: Rotate and zoom the 3D plot
- **Statistics Display**: Shows total points and pixel count
- **High-resolution Output**: 300 DPI PNG files

## Sampling Algorithms

### Fibonacci Sphere Sampling
- **Advantages**: Excellent uniformity, works with any number of points
- **Best for**: General purpose sampling, variable point counts
- **Mathematical Basis**: Uses golden ratio for optimal distribution

### Hopf Fibration
- **Advantages**: Mathematically elegant, perfect for spherical geometry
- **Best for**: When maximum uniformity is required
- **Mathematical Basis**: Maps S³ to S² using quaternion structure

## File Structure Example

```
nside_3/
├── healpix_data_0.txt
├── healpix_data_1.txt
├── healpix_data_2.txt
...
└── healpix_data_107.txt

nside_6/
├── healpix_data_0.txt
├── healpix_data_1.txt
...
└── healpix_data_431.txt
```

## Statistical Analysis Use Case

This system is designed for:
- **Monte Carlo Integration**: Equal representation ensures unbiased results
- **Statistical Comparison**: Same sampling pattern across different quantities
- **Reproducible Research**: Deterministic sampling for consistent results
- **Multi-resolution Analysis**: Compare results across different nside values

## Performance Notes

- **Memory Usage**: Scales with total number of points
- **File I/O**: Creates many small files (consider combining if needed)
- **Visualization**: 3D plotting can be slow for very high nside values
- **Storage**: Each pixel file contains header + data rows

## Troubleshooting

### Common Issues
1. **No data files found**: Run the C++ program first
2. **Compilation errors**: Ensure C++17 support and HEALPix headers
3. **Python import errors**: Install required packages with pip
4. **Large nside values**: May take significant time and memory

### Recommendations
- Start with small nside values (1-4) for testing
- Use moderate angular step sizes (0.1-0.5 radians)
- Consider combining output files for very high nside values
- Use the visualization script to verify sampling quality

## Mathematical Details

- **Coordinate System**: HEALPix uses colatitude (θ) and longitude (φ)
- **Units**: Internal calculations in radians, output in degrees
- **Pixel Boundaries**: Angular step size constrains points within pixel area
- **Equal Area**: All HEALPix pixels have equal surface area
