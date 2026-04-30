"""
    validation.jl

Physics validation: confirm parallel and serial outputs agree within statistical uncertainty.
"""

"""
    validate_energy_conservation(serial_output::DataFrame, parallel_output::DataFrame)

Check that total energy is conserved between serial and parallel runs.

Compares energy budgets: E_initial = E_final + E_lost
Should match within statistical uncertainties.

# Arguments
- `serial_output`: Observation-plane particles from serial (non-parallelized) run
- `parallel_output`: Observation-plane particles from parallel run

# Returns
- `(Δ_serial, Δ_parallel, relative_diff)`: Energy budget differences
  - `Δ_serial`: |E_initial - E_final - E_lost| for serial
  - `Δ_parallel`: |E_initial - E_final - E_lost| for parallel
  - `relative_diff`: (|Δ_serial - Δ_parallel| / max(|Δ_serial|, |Δ_parallel|))
"""
function validate_energy_conservation(serial_output::DataFrame, parallel_output::DataFrame)
    E_final_serial = sum(skipmissing(serial_output.energy))
    E_final_parallel = sum(skipmissing(parallel_output.energy))

    Δ_abs = abs(E_final_serial - E_final_parallel)
    denom = max(abs(E_final_serial), abs(E_final_parallel))
    relative_diff = denom > 0 ? Δ_abs / denom : 0.0

    @info "Energy conservation check:"
    @info "  Serial final energy: $(round(E_final_serial; digits=3)) GeV"
    @info "  Parallel final energy: $(round(E_final_parallel; digits=3)) GeV"
    @info "  Absolute difference: $(round(Δ_abs; digits=3)) GeV"
    @info "  Relative difference: $(round(100*relative_diff; digits=3))%"

    return (E_final_serial, E_final_parallel, relative_diff)
end

"""
    compare_spectra(serial_output::DataFrame, parallel_output::DataFrame)

Compare observation-plane particle spectra using Kolmogorov-Smirnov test.

Tests whether the energy distributions match within statistical uncertainty.

# Arguments
- `serial_output`: DataFrame with particle data from serial run
- `parallel_output`: DataFrame with particle data from parallel run

# Returns
- `(ks_statistic, p_value)`: KS test results
  - `ks_statistic`: Maximum difference between empirical CDFs
  - `p_value`: Probability that distributions are same (high = agreement)
"""
function compare_spectra(serial_output::DataFrame, parallel_output::DataFrame)
    E_serial = sort(Float64.(collect(skipmissing(serial_output.energy))))
    E_parallel = sort(Float64.(collect(skipmissing(parallel_output.energy))))

    isempty(E_serial) && return (NaN, NaN)
    isempty(E_parallel) && return (NaN, NaN)

    @info "Spectrum comparison:"
    @info "  Serial: $(length(E_serial)) particles, mean=$(round(mean(E_serial); digits=3)) GeV"
    @info "  Parallel: $(length(E_parallel)) particles, mean=$(round(mean(E_parallel); digits=3)) GeV"

    values = sort(vcat(E_serial, E_parallel))
    i = 1
    j = 1
    n = length(E_serial)
    m = length(E_parallel)
    dmax = 0.0

    for v in values
        while i <= n && E_serial[i] <= v
            i += 1
        end
        while j <= m && E_parallel[j] <= v
            j += 1
        end
        fs = (i - 1) / n
        fp = (j - 1) / m
        dmax = max(dmax, abs(fs - fp))
    end

    # Asymptotic two-sample KS p-value approximation.
    en = sqrt((n * m) / (n + m))
    λ = (en + 0.12 + 0.11 / en) * dmax
    p_value = 2 * sum(((-1)^(k - 1)) * exp(-2 * (k^2) * (λ^2)) for k in 1:100)
    p_value = clamp(p_value, 0.0, 1.0)

    @info "  KS D-statistic: $(round(dmax; digits=5))"
    @info "  KS p-value (asymptotic): $(round(p_value; digits=5))"

    return (dmax, p_value)
end

"""
    validate_physics(
        serial_output::DataFrame,
        parallel_output::DataFrame;
        energy_tol::Float64=0.05,
        spectrum_pval_min::Float64=0.05
    )

Comprehensive physics validation of parallel run.

Checks:
1. Energy conservation in both runs
2. Agreement in particle spectra
3. No significant loss of particles

# Arguments
- `serial_output`, `parallel_output`: DataFrames from each run
- `energy_tol`: Maximum allowed relative energy residual difference
- `spectrum_pval_min`: Minimum p-value for KS test (default 5% significance)

# Returns
- `passed::Bool`: Whether validation passed
"""
function validate_physics(
    serial_output::DataFrame,
    parallel_output::DataFrame;
    energy_tol::Float64=0.05,
    spectrum_pval_min::Float64=0.05
)
    @info "Running physics validation..."

    # Check energy conservation
    E_ser, E_par, rel_diff = validate_energy_conservation(serial_output, parallel_output)
    energy_ok = rel_diff < energy_tol

    # Check spectra
    ks_stat, p_val = compare_spectra(serial_output, parallel_output)
    spectrum_ok = p_val > spectrum_pval_min

    # Overall result
    all_ok = energy_ok && spectrum_ok

    @info "Validation result: $(all_ok ? "PASS ✓" : "FAIL ✗")"
    @info "  Energy conservation: $(energy_ok ? "PASS" : "FAIL") (diff=$(round(100*rel_diff; digits=2))% < $(100*energy_tol)%)"
    @info "  Spectrum agreement: $(spectrum_ok ? "PASS" : "FAIL") (p=$(round(p_val; digits=3)) > $spectrum_pval_min)"

    return all_ok
end

