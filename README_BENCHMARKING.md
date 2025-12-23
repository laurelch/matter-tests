# Snow Benchmark System - Quick Start

## ✅ What Was Done

Your Matter simulation has been updated to benchmark **only the snow scenario** with comprehensive CSV output tracking:

-   End frame count
-   Grid size (total nodes)
-   Number of particles
-   Steps from frame 1 to 20
-   Total simulation steps
-   Total execution time

Results are saved to **timestamped CSV files** in format `YYMMDD_HHMM.csv` with **append-only mode** (no file recreation).

## 🚀 Quick Start

### 1. Run the Benchmark Suite

```bash
cd /Users/chengguizihan/Documents/Research/2024/webassembly/23-matter-tests
./src/mpm
```

This will:

-   Test 12 particle counts (1k through 100k)
-   Test 8 thread counts (1-8)
-   Run 96 total simulations
-   Save results to `YYMMDD_HHMM.csv`

### 2. Check Your Results

```bash
# See the CSV filename that was created
ls -lt *.csv | head -1

# View the results
cat 241223_1430.csv  # Replace with your actual filename
```

### 3. CSV Output Format

```csv
end_frame,grid_size,num_particles,steps_frame_1_to_20,total_steps,total_time_ms
20,32768,1000,450,473,2719
20,32768,1000,438,461,1405
...
```

## 📊 What Changed

### Removed

-   ❌ `sand_collision()` benchmark
-   ❌ `bouncy_cube()` benchmark
-   ❌ `no_plasticity()` benchmark
-   ❌ Old CSV format with just time values

### Added

-   ✅ Timestamped CSV files (`YYMMDD_HHMM.csv`)
-   ✅ Comprehensive metrics in CSV
-   ✅ Append-only mode (never recreates files)
-   ✅ Progress indicators in console
-   ✅ Public getters for simulation metrics

### Optimized

-   ✅ Disabled PLY file output during benchmarking
-   ✅ Disabled grid data saving
-   ✅ Reduced console verbosity
-   ✅ Focus on timing accuracy

## 📖 Documentation

Three detailed guides have been created:

1. **`BENCHMARK_GUIDE.md`** - Complete usage guide

    - How to run benchmarks
    - How to customize tests
    - How to analyze results with Python/CLI
    - Troubleshooting tips

2. **`CHANGES_SUMMARY.md`** - Technical changes

    - All modified files
    - Code changes explained
    - API changes
    - Build notes

3. **`README_BENCHMARKING.md`** (this file) - Quick start

## ⚡ Quick Examples

### Run Fewer Tests (Faster)

Edit `src/mpm.cpp` main():

```cpp
int main() {
    // Test only 1k particles with different thread counts
    for (int threads = 1; threads <= 8; threads++) {
        record_to_csv(snow(1000, 20, threads));
    }
    return 0;
}
```

Then rebuild and run:

```bash
make -j 8
./src/mpm
```

### Analyze Results with Python

```python
import pandas as pd

# Load results
df = pd.read_csv('241223_1430.csv',  # Your actual filename
                 names=['end_frame', 'grid_size', 'num_particles',
                        'steps_1_to_20', 'total_steps', 'time_ms'])

# Calculate speedup
df['threads'] = df.index % 8 + 1
for p in df['num_particles'].unique():
    subset = df[df['num_particles'] == p]
    baseline = subset[subset['threads'] == 1]['time_ms'].values[0]
    print(f"{p} particles: {baseline/subset['time_ms'].mean():.2f}x speedup")
```

## 🔧 If You Need to Rebuild

```bash
cd /Users/chengguizihan/Documents/Research/2024/webassembly/23-matter-tests

# Clean
rm -rf CMakeCache.txt CMakeFiles src/CMakeFiles tests/CMakeFiles

# Configure (with OpenMP for macOS)
cmake -DCMAKE_BUILD_TYPE=Release -DUSE_VDB=OFF \
  -DOpenMP_CXX_FLAGS="-Xclang -fopenmp" \
  -DOpenMP_CXX_LIB_NAMES=libomp \
  -DOpenMP_libomp_LIBRARY=/usr/local/Cellar/libomp/19.1.7/lib/libomp.dylib \
  -DOpenMP_CXX_INCLUDE_DIR=/usr/local/Cellar/libomp/19.1.7/include \
  .

# Build
make -j 8
```

## 📝 CSV Columns Explained

| Column                | Description                           | Example |
| --------------------- | ------------------------------------- | ------- |
| `end_frame`           | Number of frames simulated            | 20      |
| `grid_size`           | Total grid nodes (Nx×Ny×Nz)           | 32768   |
| `num_particles`       | Particle count                        | 1000    |
| `steps_frame_1_to_20` | Approximate steps between frames 1-20 | 450     |
| `total_steps`         | Total simulation timesteps            | 473     |
| `total_time_ms`       | Execution time in milliseconds        | 2719    |

## ⏱️ Expected Runtime

Approximate times for full benchmark suite (96 tests):

-   **Apple M1/M2**: ~2-4 hours
-   **Intel i7/i9**: ~4-8 hours
-   **Older systems**: 8+ hours

Start with a smaller test first to estimate!

## 🐛 Troubleshooting

### "Failed to open CSV"

-   Check write permissions in the current directory
-   Make sure you're running from the project root

### Build errors

-   Ensure OpenMP is installed: `brew list libomp`
-   Try the full rebuild commands above
-   Check `BENCHMARK_GUIDE.md` for detailed troubleshooting

### Very slow execution

-   This is normal for large particle counts
-   Consider testing fewer particle counts first
-   Check CPU temperature/throttling

## 📧 Notes

-   The benchmark executable is: `./src/mpm`
-   Results append to the same CSV file for the entire session
-   Each time you start `./src/mpm`, a new timestamped CSV is created
-   Grid size may vary slightly between runs due to adaptive meshing
-   The `steps_frame_1_to_20` is an approximation: `(total_steps × 19) / end_frame`

## ✨ You're All Set!

Everything is configured and ready to go. Just run:

```bash
./src/mpm
```

Results will be saved to a timestamped CSV file. The filename will be displayed when the program starts.

Happy benchmarking! 🎉
