# CorsikaParallelism.jl - Testing Instructions

## Setup

First do this: 

export LD_LIBRARY_PATH=/n/home09/tkrishnan/.julia/juliaup/julia-1.12.6+0.x64.linux.gnu/lib/julia:$LD_LIBRARY_PATH 
export JULIA_DEPOT_PATH=/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/.julia_cache 

### 1. Activate the package in Julia

```julia
using Pkg
Pkg.activate("/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/CorsikaParallelism.jl")
Pkg.instantiate()  # Install dependencies
```

### 2. Load the package

```julia
using CorsikaParallelism
using DataFrames, Arrow, Distributed
```

### 3. Add worker processes (optional but recommended)

```julia
addprocs(4)  # Add 4 worker processes for parallel execution
```

---

## Test 1: Run Primary Shower with First-Interaction Stopping

```julia
# Setup
corsika_binary = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-build/applications/c8_air_shower"
output_dir = "/tmp/corsika_test_primary"
mkpath(output_dir)

# Primary config
config = Dict(
    "pdg" => 2212,         # proton (strongly interacting, good for producing secondaries)
    "energy" => 1000.0,    # 10 TeV
    "zenith" => 0.0,       # straight down
    "azimuth" => 0.0,
    "nevent" => 1,         # run a few events to reliably get non-empty secondaries
    "seed" => 12345,       # fixed seed for reproducible test runs
    "eslope" => -1.0,
    "emcut" => 0.3         # 300 MeV cutoff
)

# Keep these in global scope for later tests in REPL sessions
output_file = ""
secondaries = DataFrame()

# Run with first-interaction stopping
try
    global output_file = CorsikaParallelism.run_corsika_with_stopping(
        corsika_binary, 
        config, 
        output_dir;
        generation=0  # Generation 0 = primary splits
    )
    println("✓ Primary shower complete: $output_file")
    
    # Read secondaries
    global secondaries = CorsikaParallelism.read_secondary_particles(output_file)
    println("✓ Read $(size(secondaries, 1)) secondary particles")
    if size(secondaries, 1) > 0 && hasproperty(secondaries, :energy)
        println("  Energy range: $(minimum(secondaries.energy)) - $(maximum(secondaries.energy)) GeV")
    else
        error("No usable secondaries produced. Increase nevent (e.g. 10-20) or raise primary energy.")
    end
    
catch e
    println("✗ Primary shower failed: $e")
end
```

---

## Test 2: Test All Four Schedulers

```julia
# Use secondaries from Test 1
n_workers = 4
size(secondaries, 1) == 0 && error("No secondaries found from Test 1. Re-run Test 1 first.")

schedulers = [:naive, :threshold, :binpack, :workstealing]

for sched in schedulers
    println("\n=== Testing :$sched scheduler ===")
    
    try
        batches = CorsikaParallelism.schedule_secondaries(
            secondaries, 
            sched, 
            n_workers
        )
        
        # Show distribution
        batch_sizes = [size(b, 1) for b in batches]
        batch_energies = [sum(skipmissing(b.energy)) for b in batches]
        
        println("✓ Batch sizes: $batch_sizes")
        println("✓ Total energy per worker: $(round.(batch_energies; digits=1)) GeV")
        
        # Compute load balance metric (lower is better)
        if sum(batch_sizes) > 0
            avg_size = sum(batch_sizes) / length(batch_sizes)
            imbalance = maximum(batch_sizes) / avg_size
            println("✓ Load imbalance ratio: $(round(imbalance; digits=2))")
        end
        
    catch e
        println("✗ Scheduler failed: $e")
    end
end
```

---

## Test 3: Power-Law Runtime Estimation

```julia
using Statistics

# Example: profile runtimes at different energies (if you have profiling data)
# For now, just test the functions exist:

energy = 100.0  # GeV
alpha = 2.0     # Power-law exponent
baseline = 10.0 # 10 seconds at 1 GeV

# Estimate runtime for a 100 GeV particle
est_time = CorsikaParallelism.estimate_runtime(energy, alpha, baseline)
println("Estimated runtime for $energy GeV particle: $(round(est_time; digits=1)) seconds")

# Test power-law fitting (with synthetic data)
E_test = [10.0, 50.0, 100.0, 500.0, 1000.0]
t_test = 0.1 .* (E_test .^ 1.8)  # Synthetic: t = 0.1 * E^1.8

alpha_fit, baseline_fit = CorsikaParallelism.fit_powerlaw(E_test, t_test)
println("✓ Fitted power-law: t = $(round(baseline_fit; digits=3)) * E^$(round(alpha_fit; digits=2))")
```

---

## Test 4: Full Orchestration (Parallel Shower Splitting)

**IMPORTANT:** This runs CORSIKA multiple times in parallel. Takes significant time and compute.

```julia
# Only run if you have sufficient compute resources
# Recommended: 8+ cores, 32GB+ memory
@everywhere using CorsikaParallelism

if nprocs() < 5
    addprocs(4)
end

try
    # Run full parallel orchestration
    results = CorsikaParallelism.distribute_and_run_secondaries(
        secondaries,
        corsika_binary,
        config,
        :naive;  # Use best-performing scheduler
        nworkers=Distributed.nworkers()
    )
    
    println("✓ Parallel execution complete")
    println("✓ Aggregated results from $(length(results)) worker batches")
    
catch e
    println("✗ Orchestration failed: $e")
    println(stacktrace(catch_backtrace()))
end
```

---

## Test 5: Physics Validation (Serial vs Parallel)

```julia
# After running both serial and parallel versions
# Read their observation-plane outputs

serial_output = DataFrame(Arrow.Table("path/to/serial_output.parquet"))
parallel_output = DataFrame(Arrow.Table("path/to/parallel_output.parquet"))

# Validate
passed = CorsikaParallelism.validate_physics(
    serial_output,
    parallel_output;
    energy_tol=0.10,        # 10% tolerance
    spectrum_pval_min=0.05  # 5% significance
)

if passed
    println("✓ Physics validation PASSED")
else
    println("✗ Physics validation FAILED - check output above")
end
```

---

## Quick Start Script

Save this as `test_package.jl`:

```julia
using Pkg
Pkg.activate("/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/CorsikaParallelism.jl")
Pkg.instantiate()

using CorsikaParallelism, DataFrames, Distributed

corsika_binary = "/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-build/applications/c8_air_shower"
output_dir = "/tmp/corsika_test"
mkpath(output_dir)

config = Dict(
    "pdg" => 16,
    "energy" => 1000.0,
    "zenith" => 0.0,
    "azimuth" => 0.0,
    "nevent" => 1
)

println("=== Test 1: Primary Shower ===")
output_file = ""
secondaries = DataFrame()
try
    global output_file = CorsikaParallelism.run_corsika_with_stopping(
        corsika_binary, config, output_dir; generation=0
    )
    global secondaries = CorsikaParallelism.read_secondary_particles(output_file)
    println("✓ Success: $(size(secondaries, 1)) secondaries")
    
    println("\n=== Test 2: Schedulers ===")
    for sched in [:naive, :threshold, :binpack, :workstealing]
        batches = CorsikaParallelism.schedule_secondaries(secondaries, sched, 4)
        sizes = [size(b, 1) for b in batches]
        println("  $sched: $sizes particles")
    end
    
catch e
    println("✗ Error: $e")
    stacktrace(catch_backtrace())
end
```

Run with:
```bash
julia test_package.jl
```

---

## Environment Notes

- **CORSIKA binary paths:**
  - Modified: `/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-build/applications/c8_air_shower`
  - Baseline: `/n/holylfs05/LABS/arguelles_delgado_lab/Everyone/tkrishnan/corsika-baseline-build/applications/c8_air_shower`

- **Julia environment:** Use SLURM for heavy tests:
  ```bash
  srun --account=arguelles_delgado_lab --partition=sapphire,shared --mem=32000 --time=120 julia test_package.jl
  ```

- **Expected output files:** Parquet format in the `output_dir` specified

---

## Troubleshooting

| Problem | Solution |
|---------|----------|
| `libstdc++.so.6` version error | Run from compute node with GCC 12.2.0 environment set |
| CORSIKA binary not found | Check path is correct and binary is executable |
| Out of memory | Reduce `nevent` or use fewer parallel workers |
| Parquet file not found | Check CORSIKA ran successfully and `--output` dir exists |
| `--filename: Path already exists` | Delete old output (`rm -rf /tmp/corsika_test_primary/shower_gen0_e1000.0.parquet`) or rerun with updated package code that auto-cleans old outputs |
| Package not loading | Run `Pkg.instantiate()` to install dependencies |
