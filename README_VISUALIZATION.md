# HEALPix Sampling Boundary Verification

This directory contains tools to verify that the C++ HEALPix sampling code properly bounds sampling points within HEALPix pixel boundaries.

## Files

- `main2.cpp` - C++ code for HEALPix sampling using Fibonacci or Hopf fibration methods
- `visualize_healpix_sampling.py` - Full-featured visualization script with Mollweide projections and 3D plots
- `quick_verify.py` - Lightweight verification script for fast boundary checking
- `test_boundaries.py` - Interactive tool that runs the full pipeline (C++ + verification)
- `requirements.txt` - Python package dependencies

## Quick Start

### 1. Install Dependencies

```bash
pip install -r requirements.txt
```

### 2. Compile the C++ Code

```bash
g++ -std=c++11 -O2 main2.cpp src/Healpix_2.15a/*.cc -o main2
```

### 3. Run Complete Verification

```bash
python test_boundaries.py
```

This interactive script will:
- Run the C++ sampling code with your chosen parameters
- Verify that all points are bounded within their HEALPix pixels
- Display visualizations showing the results

### 4. Manual Usage

#### Generate HEALPix Samples
```bash
./main2
# Enter nside (e.g., 1)
# Enter angular step size (e.g., 5.7)
# Choose sampling method (1=Fibonacci, 2=Hopf)
```

#### Quick Verification
```bash
python quick_verify.py 1
```

#### Full Visualization
```bash
python visualize_healpix_sampling.py 1
```

## Script Details

### `visualize_healpix_sampling.py`

Full-featured visualization with multiple options:

```bash
python visualize_healpix_sampling.py <nside> [options]

Options:
  --data-dir DIR        Specify data directory (default: nside_<nside>)
  --no-boundaries       Skip pixel boundaries (faster for large nside)
  --3d                  Show 3D sphere plot
  --verify-only         Only verify containment, no plots
  --save PREFIX         Save plots to files
```

Examples:
```bash
# Basic visualization with boundaries
python visualize_healpix_sampling.py 1

# Large nside without boundaries (faster)
python visualize_healpix_sampling.py 8 --no-boundaries

# 3D visualization
python visualize_healpix_sampling.py 2 --3d

# Save plots to files
python visualize_healpix_sampling.py 1 --save healpix_nside1

# Just verify containment
python visualize_healpix_sampling.py 4 --verify-only
```

### `quick_verify.py`

Lightweight verification for fast checking:

```bash
python quick_verify.py <nside> [data_dir]
```

This script:
- Loads all sampling data for the given nside
- Verifies each point belongs to its assigned pixel using `healpy.ang2pix()`
- Reports accuracy statistics
- Shows simple scatter plots of correct vs misassigned points

### Output Interpretation

#### Verification Results

The scripts will report:
- **Total points**: Number of sampling points generated
- **Correct assignments**: Points that belong to their assigned pixel
- **Wrong assignments**: Points incorrectly assigned (should be 0 for proper bounded sampling)
- **Accuracy**: Percentage of correctly assigned points (should be 100% for proper bounds)

#### Visual Output

1. **Mollweide Projection**: Shows sampling points on a world map projection
   - Left panel: Points only, colored by pixel
   - Right panel: Points with HEALPix pixel boundaries

2. **3D Sphere**: Interactive 3D view of the sampling on the unit sphere

3. **Verification Plot**: 
   - Correct points (blue)
   - Misassigned points (red) - should be none for proper bounds

## Expected Results

For properly bounded sampling:
- **Accuracy should be 100%** - all points belong to their assigned pixels
- **Visual verification**: Points should lie within their pixel boundaries
- **No red points** in the verification scatter plot

## Troubleshooting

### Common Issues

1. **"Directory nside_X not found"**
   - Run the C++ code first to generate sampling data

2. **"Missing required package"**
   - Install dependencies: `pip install -r requirements.txt`

3. **"No compiled executable found"**
   - Compile the C++ code: `g++ -std=c++11 -O2 main2.cpp src/Healpix_2.15a/*.cc -o main2`

4. **Slow visualization for large nside**
   - Use `--no-boundaries` flag
   - Try `--verify-only` for just numerical verification

### Performance Notes

- Boundary visualization can be slow for nside > 4
- Use `quick_verify.py` for fast verification of large datasets
- 3D plots are always computed without boundaries for performance

## Technical Details

### Coordinate Systems

- **C++ code**: Uses standard spherical coordinates (θ, φ) in radians
  - θ ∈ [0, π] (colatitude, 0 at north pole)
  - φ ∈ [0, 2π] (azimuth)

- **HEALPix**: Uses the same convention
  - Pixels numbered 0 to 12×nside²-1

- **Visualizations**: Convert to longitude/latitude for plotting
  - Longitude = φ - π ∈ [-π, π]
  - Latitude = π/2 - θ ∈ [-π/2, π/2]

### Verification Method

The verification uses `healpy.ang2pix(nside, theta, phi)` to:
1. Load each point from the C++ output files
2. Compute which HEALPix pixel it should belong to
3. Compare with the pixel it was assigned to in the filename
4. Report any mismatches as boundary violations

A properly bounded sampling should have 100% accuracy in this test.