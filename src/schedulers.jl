"""
    schedulers.jl

Four scheduling strategies for distributing secondaries to parallel workers.
Each returns a Vector{DataFrame}, one batch per worker.
"""

"""
    schedule_secondaries_naive(particles::DataFrame, nworkers::Int)

Naive strategy: Round-robin distribution.
Simple but may be imbalanced due to energy variation.
"""
function schedule_secondaries_naive(particles::DataFrame, nworkers::Int)
    worker_batches = [DataFrame() for _ in 1:nworkers]

    for (i, row) in enumerate(eachrow(particles))
        worker_id = (i - 1) % nworkers + 1
        push!(worker_batches[worker_id], row)
    end

    return worker_batches
end

"""
    schedule_secondaries_threshold(particles::DataFrame, nworkers::Int)

Threshold strategy: Bin by energy, distribute bins to workers.
Groups similar-energy particles together, which may improve cache locality.

Energy thresholds are adaptive based on particle energy distribution.
"""
function schedule_secondaries_threshold(particles::DataFrame, nworkers::Int)
    # Sort by energy
    sorted_df = sort(particles, :energy)

    # Compute energy quantiles for bin boundaries
    quantiles = range(0, 1, length=nworkers+1)
    thresholds = quantile(sorted_df.energy, quantiles)

    # Assign particles to workers based on energy threshold
    worker_batches = [DataFrame() for _ in 1:nworkers]

    for row in eachrow(sorted_df)
        E = row.energy
        # Find which bin this energy falls into
        bin_idx = searchsortedlast(thresholds, E)
        bin_idx = clamp(bin_idx, 1, nworkers)
        push!(worker_batches[bin_idx], row)
    end

    return worker_batches
end

"""
    schedule_secondaries_binpack(particles::DataFrame, nworkers::Int; alpha::Float64=2.0)

Bin-packing strategy: Greedily pack particles to balance estimated runtime.

Assumes runtime ∝ E^α (power-law). Uses greedy algorithm: assign each particle
to the worker with smallest current estimated runtime.

# Arguments
- `particles`: Secondary particles
- `nworkers`: Number of workers
- `alpha`: Power-law exponent for runtime estimation (default 2.0)
"""
function schedule_secondaries_binpack(particles::DataFrame, nworkers::Int; alpha::Float64=2.0)
    # Track estimated runtime per worker (in arbitrary units)
    worker_loads = zeros(Float64, nworkers)
    worker_batches = [DataFrame() for _ in 1:nworkers]

    # Sort by energy (descending) for better packing
    sorted_df = sort(particles, :energy, rev=true)

    # Greedy assignment: assign each particle to least-loaded worker
    for row in eachrow(sorted_df)
        # Estimated runtime ∝ E^α
        particle_runtime = row.energy ^ alpha

        # Find worker with minimum current load
        min_load, min_worker = findmin(worker_loads)

        # Assign to that worker
        push!(worker_batches[min_worker], row)
        worker_loads[min_worker] += particle_runtime
    end

    return worker_batches
end

"""
    schedule_secondaries_workstealing(particles::DataFrame, nworkers::Int; alpha::Float64=2.0)

Work-stealing strategy: Initial balanced distribution, with mechanism for dynamic stealing.

Uses greedy bin-packing for initial assignment, but marks particles as "stealable"
so a faster-finishing worker can steal particles from slower workers.

Returns the same initial allocation as :binpack.
Dynamic work-stealing is implemented at runtime in orchestrator.jl.
"""
function schedule_secondaries_workstealing(particles::DataFrame, nworkers::Int; alpha::Float64=2.0)
    # Use bin-packing as initial distribution
    worker_batches = schedule_secondaries_binpack(particles, nworkers; alpha=alpha)

    return worker_batches
end

"""
    schedule_secondaries(
        particles::DataFrame,
        scheduler::Symbol,
        nworkers::Int;
        kwargs...
    )::Vector{DataFrame}

Main scheduling dispatcher. Routes to appropriate strategy.

# Arguments
- `particles`: DataFrame of secondary particles
- `scheduler`: Strategy symbol (:naive, :threshold, :binpack, :workstealing)
- `nworkers`: Number of parallel workers
- `kwargs`: Additional arguments (e.g., alpha for power-law strategies)

# Returns
- Vector of DataFrames, one batch per worker
"""
function schedule_secondaries(
    particles::DataFrame,
    scheduler::Symbol,
    nworkers::Int;
    kwargs...
)::Vector{DataFrame}

    if isempty(particles)
        return [DataFrame() for _ in 1:nworkers]
    end

    # Validate scheduler choice
    valid_schedulers = (:naive, :threshold, :binpack, :workstealing)
    !(scheduler in valid_schedulers) && error(
        "Unknown scheduler :$scheduler. Must be one of $valid_schedulers"
    )

    @info "Scheduling $(nrow(particles)) particles using :$scheduler strategy"

    result = if scheduler == :naive
        schedule_secondaries_naive(particles, nworkers)
    elseif scheduler == :threshold
        schedule_secondaries_threshold(particles, nworkers)
    elseif scheduler == :binpack
        schedule_secondaries_binpack(particles, nworkers; kwargs...)
    elseif scheduler == :workstealing
        schedule_secondaries_workstealing(particles, nworkers; kwargs...)
    end

    # Log distribution
    batch_sizes = [nrow(batch) for batch in result]
    total_energy = [nrow(batch) > 0 ? sum(skipmissing(batch.energy)) : 0.0 for batch in result]
    @info "Distribution: $batch_sizes particles, $(round.(total_energy; digits=1)) GeV per worker"

    return result
end

