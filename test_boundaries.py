#!/usr/bin/env python3
"""
Simple usage example for HEALPix boundary verification.
Run this after generating data with the C++ code.
"""

import subprocess
import sys
import os

def check_dependencies():
    """Check if required Python packages are installed."""
    try:
        import numpy
        import matplotlib
        import healpy
        print("✅ All required packages are installed")
        return True
    except ImportError as e:
        print(f"❌ Missing required package: {e}")
        print("Please install requirements: pip install -r requirements.txt")
        return False

def run_cpp_sampling(nside=1, angular_step=5.7, method=1):
    """
    Run the C++ sampling code with default parameters.
    
    Parameters:
    -----------
    nside : int
        HEALPix nside parameter
    angular_step : float
        Angular step size in degrees
    method : int
        Sampling method (1=Fibonacci, 2=Hopf)
    """
    print(f"Running C++ code with nside={nside}, angular_step={angular_step}°, method={method}")
    
    # Check if executable exists
    if os.path.exists("./main2"):
        executable = "./main2"
    elif os.path.exists("./main"):
        executable = "./main"
    else:
        print("❌ No compiled executable found. Please compile the C++ code first:")
        print("   g++ -std=c++11 -O2 main2.cpp src/Healpix_2.15a/*.cc -o main2")
        return False
    
    try:
        # Run with input
        input_str = f"{nside}\n{angular_step}\n{method}\n"
        result = subprocess.run([executable], input=input_str, text=True, 
                              capture_output=True, timeout=30)
        
        if result.returncode == 0:
            print("✅ C++ sampling completed successfully")
            print(result.stdout)
            return True
        else:
            print("❌ C++ sampling failed:")
            print(result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("❌ C++ sampling timed out")
        return False
    except Exception as e:
        print(f"❌ Error running C++ code: {e}")
        return False

def run_verification(nside):
    """Run the Python verification script."""
    print(f"\nRunning verification for nside={nside}...")
    
    try:
        # Import and run the quick verification
        sys.path.insert(0, '.')
        from quick_verify import quick_verify_and_plot
        
        result = quick_verify_and_plot(nside, show_plot=True)
        
        if result and result['accuracy'] >= 0.99:
            print("\n🎉 SUCCESS: HEALPix boundary verification passed!")
            print(f"   - {result['total_points']} points sampled across {result['pixels_with_data']} pixels")
            print(f"   - {result['accuracy']*100:.2f}% of points are correctly bounded")
            return True
        elif result:
            print(f"\n⚠️  WARNING: Only {result['accuracy']*100:.2f}% accuracy in boundary containment")
            return False
        else:
            print("\n❌ FAILED: Could not verify boundary containment")
            return False
            
    except Exception as e:
        print(f"❌ Error during verification: {e}")
        return False

def main():
    print("HEALPix Boundary Verification Tool")
    print("=" * 40)
    
    # Check dependencies
    if not check_dependencies():
        return 1
    
    # Get parameters
    try:
        nside = int(input("Enter nside (default: 1): ") or "1")
        angular_step = float(input("Enter angular step in degrees (default: 5.7): ") or "5.7")
        method = int(input("Enter method 1=Fibonacci, 2=Hopf (default: 1): ") or "1")
    except ValueError:
        print("❌ Invalid input. Using defaults: nside=1, angular_step=5.7°, method=1")
        nside, angular_step, method = 1, 5.7, 1
    
    print(f"\nConfiguration:")
    print(f"  - nside: {nside}")
    print(f"  - angular step: {angular_step}°")
    print(f"  - method: {'Fibonacci' if method == 1 else 'Hopf'}")
    print()
    
    # Step 1: Generate data (optional, if not already exists)
    data_dir = f"nside_{nside}"
    if not os.path.exists(data_dir):
        print("Data directory not found. Running C++ sampling...")
        if not run_cpp_sampling(nside, angular_step, method):
            return 1
    else:
        print(f"✅ Found existing data in {data_dir}")
    
    # Step 2: Verify boundaries
    if run_verification(nside):
        print("\n✨ Verification completed successfully!")
        print("The visualization window should show the sampling points and boundary verification.")
        return 0
    else:
        return 1

if __name__ == "__main__":
    sys.exit(main())