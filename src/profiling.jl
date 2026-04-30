"""
    profiling.jl

Runtime profiling and empirical power-law fitting for secondaries.
Fit empirical data to: t ∝ E^α
"""

"""
    estimate_runtime(energy::Float64, alpha::Float64, baseline::Float64)

Estimate secondary shower runtime from energy using power-law model.

# Arguments
- `energy`: Secondary particle energy (GeV)
- `alpha`: Power-law exponent (typically 1.5-2.0 for cascade showers)
- `baseline`: Reference runtime at 1 GeV (seconds)

# Returns
- Estimated runtime in seconds
"""
function estimate_runtime(energy::Float64, alpha::Float64, baseline::Float64)
    return baseline * (energy ^ alpha)
end

"""
    fit_powerlaw(energies::Vector, runtimes::Vector)

Fit power-law model t = baseline * E^α to empirical profiling data.

Uses log-log linear regression: log(t) = log(baseline) + α*log(E)

# Arguments
- `energies`: Vector of secondary particle energies (GeV)
- `runtimes`: Corresponding measured runtimes (seconds)

# Returns
- `(alpha, baseline)`: Power-law parameters
  - `alpha`: Exponent (slope in log-log space)
  - `baseline`: Reference runtime at 1 GeV
"""
function fit_powerlaw(energies::Vector, runtimes::Vector)
    # Validate inputs
    length(energies) != length(runtimes) && error(
        "energies and runtimes must have same length"
    )
    isempty(energies) && error("Cannot fit power-law to empty data")

    # Remove zeros/negatives (log undefined)
    valid_idx = (energies .> 0) .& (runtimes .> 0)
    E_valid = energies[valid_idx]
    t_valid = runtimes[valid_idx]

    if length(E_valid) < 2
        @warn "Insufficient valid data points for power-law fit, using defaults"
        return (2.0, 1.0)  # Default α=2.0, baseline=1s
    end

    # Log-log regression: log(t) = a + b*log(E)
    # where a = log(baseline), b = alpha
    log_E = log.(E_valid)
    log_t = log.(t_valid)

    # Least-squares fit
    n = length(E_valid)
    mean_log_E = mean(log_E)
    mean_log_t = mean(log_t)

    numerator = sum((log_E .- mean_log_E) .* (log_t .- mean_log_t))
    denominator = sum((log_E .- mean_log_E) .^ 2)

    if denominator ≈ 0
        @warn "Cannot fit power-law: energies are constant"
        return (1.0, mean(t_valid))
    end

    alpha = numerator / denominator
    log_baseline = mean_log_t - alpha * mean_log_E
    baseline = exp(log_baseline)

    # Compute R² for goodness of fit
    y_pred = baseline .* (E_valid .^ alpha)
    ss_res = sum((t_valid .- y_pred) .^ 2)
    ss_tot = sum((t_valid .- mean(t_valid)) .^ 2)
    r_squared = 1 - (ss_res / ss_tot)

    @info "Power-law fit: t = $(round(baseline; digits=3)) * E^$(round(alpha; digits=2))" R²=round(r_squared; digits=3)

    return (alpha, baseline)
end

"""
    profile_shower_runtime(
        energy_range::Vector,
        corsika_binary::String,
        config::Dict,
        n_trials::Int=3
    )

Empirically profile shower runtimes across an energy range.

Runs CORSIKA multiple times at different energies and fits power-law model.

# Arguments
- `energy_range`: Vector of energies (GeV) to profile
- `corsika_binary`: Path to CORSIKA binary
- `config`: Base configuration dict
- `n_trials`: Number of trials per energy (for averaging)

# Returns
- `(alpha, baseline)`: Fitted power-law parameters
"""
function profile_shower_runtime(
    energy_range::Vector,
    corsika_binary::String,
    config::Dict,
    n_trials::Int=3
)
    energies = Float64[]
    runtimes = Float64[]

    @info "Profiling shower runtime across $(length(energy_range)) energies, $n_trials trials each"

    for E in energy_range
        trial_times = Float64[]

        for trial in 1:n_trials
            @info "Profiling E=$E GeV (trial $trial/$n_trials)"

            # Create temporary output directory
            output_dir = mktempdir()

            try
                # Measure runtime
                start_time = time()

                trial_config = Dict{String,Any}(config)
                trial_config["energy"] = E
                trial_config["nevent"] = get(config, "nevent", 1)
                trial_config["stop_after_first_interaction"] = false
                trial_config["output_tag"] = "profile_E$(E)_trial$(trial)"

                run_corsika_with_stopping(corsika_binary, trial_config, output_dir; generation=1)

                elapsed_time = time() - start_time
                push!(trial_times, elapsed_time)

            finally
                rm(output_dir, recursive=true, force=true)
            end
        end

        # Average over trials
        mean_time = mean(trial_times)
        push!(energies, E)
        push!(runtimes, mean_time)

        @info "E=$E GeV: $(round(mean_time; digits=2))s ($(length(trial_times)) trials)"
    end

    # Fit power-law
    alpha, baseline = fit_powerlaw(energies, runtimes)

    return (alpha, baseline)
end

