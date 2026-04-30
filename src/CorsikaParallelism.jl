"""
    CorsikaParallelism

Sub-shower parallelism for PeV tau air showers using CORSIKA8 cascade splitting.

Provides:
- Orchestration of CORSIKA8 with first-interaction stopping
- Scheduling strategies for distributing secondaries to parallel workers
- Runtime profiling and load balancing
- Physics validation for parallel vs serial execution
"""
module CorsikaParallelism

using Arrow, DataFrames, Distributed, Parquet2, Statistics, StatsBase, YAML, Unitful

export
    run_corsika_with_stopping,
    read_secondary_particles,
    schedule_secondaries,
    estimate_runtime,
    validate_energy_conservation,
    compare_spectra

include("orchestrator.jl")
include("schedulers.jl")
include("profiling.jl")
include("validation.jl")

end # module CorsikaParallelism
