"""
    orchestrator.jl

Main orchestration layer for sub-shower parallelism.
Runs CORSIKA8 with first-interaction stopping, reads secondaries, and distributes to workers.
"""

function uses_geometry_interface(config::Dict)
    geometry_keys = (
        "xpos", "ypos", "zpos",
        "xdir", "ydir", "zdir",
        "x-intercept", "y-intercept", "z-intercept",
    )
    return any(haskey(config, key) for key in geometry_keys)
end

function observation_plane_basis(config::Dict)
    if haskey(config, "xdir") || haskey(config, "ydir") || haskey(config, "zdir")
        nx = Float64(get(config, "xdir", 0.0))
        ny = Float64(get(config, "ydir", 0.0))
        nz = Float64(get(config, "zdir", 1.0))
    else
        nx, ny, nz = 0.0, 0.0, 1.0
    end

    normal_norm = sqrt(nx^2 + ny^2 + nz^2)
    normal_norm == 0.0 && error("Observation plane normal cannot be zero")
    normal = (nx / normal_norm, ny / normal_norm, nz / normal_norm)

    xaxis = if abs(ny) > 0.0 || abs(nz) > 0.0
        denom = sqrt(ny^2 + nz^2)
        (0.0, -nz / denom, ny / denom)
    else
        (0.0, 1.0, 0.0)
    end

    yaxis = (
        normal[2] * xaxis[3] - normal[3] * xaxis[2],
        normal[3] * xaxis[1] - normal[1] * xaxis[3],
        normal[1] * xaxis[2] - normal[2] * xaxis[1],
    )

    return xaxis, yaxis, normal
end

function observation_plane_center_km(config::Dict)
    if uses_geometry_interface(config)
        return (
            Float64(get(config, "x-intercept", 0.0)),
            Float64(get(config, "y-intercept", 0.0)),
            Float64(get(config, "z-intercept", 0.0)),
        )
    end

    return (
        0.0,
        0.0,
        Float64(get(config, "observation-level", 0.0)) / 1000.0,
    )
end

function local_to_surface_coords_km(local_xyz_m, config::Dict)
    xaxis, yaxis, normal = observation_plane_basis(config)
    center = observation_plane_center_km(config)
    scale = 1 / 1000.0

    return (
        center[1] + scale * (local_xyz_m[1] * xaxis[1] + local_xyz_m[2] * yaxis[1] + local_xyz_m[3] * normal[1]),
        center[2] + scale * (local_xyz_m[1] * xaxis[2] + local_xyz_m[2] * yaxis[2] + local_xyz_m[3] * normal[2]),
        center[3] + scale * (local_xyz_m[1] * xaxis[3] + local_xyz_m[2] * yaxis[3] + local_xyz_m[3] * normal[3]),
    )
end

function local_to_global_direction(local_xyz, config::Dict)
    xaxis, yaxis, normal = observation_plane_basis(config)
    return (
        local_xyz[1] * xaxis[1] + local_xyz[2] * yaxis[1] + local_xyz[3] * normal[1],
        local_xyz[1] * xaxis[2] + local_xyz[2] * yaxis[2] + local_xyz[3] * normal[2],
        local_xyz[1] * xaxis[3] + local_xyz[2] * yaxis[3] + local_xyz[3] * normal[3],
    )
end

"""
    run_corsika_with_stopping(
        corsika_binary::String,
        config::Dict,
        output_dir::String;
        generation::Int=0
    )

Run the modified CORSIKA8 binary with first-interaction stopping enabled.

# Arguments
- `corsika_binary`: Path to the modified c8_air_shower executable
- `config`: Configuration dict with primary particle settings (energy, theta, etc.)
- `output_dir`: Directory where CORSIKA outputs (Parquet, logs) will be written
- `generation`: Generation depth (0 = primary splits, 1+ = from interactions)

# Returns
- `output_file`: Path to the primary shower Parquet output file
"""
function run_corsika_with_stopping(
    corsika_binary::String,
    config::Dict,
    output_dir::String;
    generation::Int=0
)
    config = Dict{String,Any}(config)
    # Validate inputs
    !isfile(corsika_binary) && error("CORSIKA binary not found: $corsika_binary")
    !isdir(output_dir) && mkpath(output_dir)

    # Build command line arguments
    cmd_args = String[]

    # Primary particle settings
    push!(cmd_args, "-p", string(config["pdg"]))
    push!(cmd_args, "-E", string(config["energy"]))
    push!(cmd_args, "-N", string(config["nevent"]))

    if uses_geometry_interface(config)
        for key in (
            "xpos", "ypos", "zpos",
            "xdir", "ydir", "zdir",
            "x-intercept", "y-intercept", "z-intercept",
        )
            haskey(config, key) && push!(cmd_args, "--$key", string(config[key]))
        end
    else
        push!(cmd_args, "-z", string(config["zenith"]))
        push!(cmd_args, "-a", string(config["azimuth"]))
    end

    # Output filename (CORSIKA requires --filename, not --output dir)
    base_name = haskey(config, "output_tag") ? string(config["output_tag"]) : "shower_gen$(generation)_e$(config["energy"])"
    output_file = joinpath(output_dir, "$(base_name).parquet")
    if ispath(output_file)
        output_file = joinpath(output_dir, "$(base_name)_$(time_ns()).parquet")
    end
    push!(cmd_args, "--filename", output_file)

    # Generation control
    push!(cmd_args, "--generation", string(generation))
    if get(config, "stop_after_first_interaction", generation == 0)
        push!(cmd_args, "--stop-after-first-interaction")
    end

    # Additional config options
    if haskey(config, "eslope")
        push!(cmd_args, "--eslope", string(config["eslope"]))
    end
    if haskey(config, "emcut")
        push!(cmd_args, "--emcut", string(config["emcut"]))
    end
    if haskey(config, "hadcut")
        push!(cmd_args, "--hadcut", string(config["hadcut"]))
    end
    if haskey(config, "mucut")
        push!(cmd_args, "--mucut", string(config["mucut"]))
    end
    if haskey(config, "taucut")
        push!(cmd_args, "--taucut", string(config["taucut"]))
    end
    if haskey(config, "observation-level")
        push!(cmd_args, "--observation-level", string(config["observation-level"]))
    end
    if haskey(config, "injection-height")
        push!(cmd_args, "--injection-height", string(config["injection-height"]))
    end
    if haskey(config, "seed")
        push!(cmd_args, "--seed", string(config["seed"]))
    end
    if get(config, "force_interaction", false)
        push!(cmd_args, "--force-interaction")
    end
    if get(config, "force_decay", false)
        push!(cmd_args, "--force-decay")
    end
    if haskey(config, "photoncut")
        push!(cmd_args, "--photoncut", string(config["photoncut"]))
    end

    @info "Running CORSIKA8 with first-interaction stopping" generation
    @info "Command: $corsika_binary $(join(cmd_args, " "))"

    # Run CORSIKA
    run(Cmd([corsika_binary; cmd_args]))

    # CORSIKA outputs a directory structure
    !isdir(output_file) && error("CORSIKA output directory not found: $output_file")
    @info "Primary shower output: $output_file"

    return output_file
end

"""
    read_secondary_particles(output_dir::String; generation::Int=0)

Read secondary particle data from CORSIKA output directory.
For generation < 2: reads interactions (split products).
For generation >= 2: reads observation-level particles.

# Returns
- DataFrame with columns: energy, x, y, z, pdg_id, (other relevant fields)
"""
function read_secondary_particles(output_dir::String; generation::Int=0)
    parquet_file = if generation == 0
        joinpath(output_dir, "interactions", "interactions.parquet")
    else
        joinpath(output_dir, "particles", "particles.parquet")
    end

    !isfile(parquet_file) && error("Parquet file not found: $parquet_file")

    @info "Reading secondary particles from $parquet_file"

    table = Parquet2.Dataset(parquet_file)
    df = DataFrame(table)
    @info "Read $(nrow(df)) secondary particles"

    # Log energy range if column exists
    if hasproperty(df, :energy) && nrow(df) > 0
        @info "Energy range: $(minimum(df.energy)) - $(maximum(df.energy)) GeV"
    end

    return df
end

function _canonical_seed(x)
    # Keep seed in CORSIKA-friendly positive Int32 range.
    s = try
        Int(round(Float64(x)))
    catch
        1
    end
    return clamp(abs(s), 1, 2_147_483_646)
end

function derive_subshower_seed(base_seed, shower_id::Int, idx::Int, pdg::Int)
    b = UInt64(_canonical_seed(base_seed))
    # SplitMix-style mixing for reproducible distinct seeds.
    z = b + UInt64(0x9e3779b97f4a7c15)
    z ⊻= UInt64(shower_id) * UInt64(0xbf58476d1ce4e5b9)
    z ⊻= UInt64(idx) * UInt64(0x94d049bb133111eb)
    z ⊻= UInt64(abs(pdg) + 11) * UInt64(0x369dea0f31a53f85)
    z ⊻= (z >> 30)
    z *= UInt64(0xbf58476d1ce4e5b9)
    z ⊻= (z >> 27)
    z *= UInt64(0x94d049bb133111eb)
    z ⊻= (z >> 31)
    return Int(mod(z, UInt64(2_147_483_646)) + 1)
end

function _secondary_task_from_row(row)
    return (
        pdg = Int(row.pdg),
        energy = Float64(row.energy),
        px = Float64(row.px),
        py = Float64(row.py),
        pz = Float64(row.pz),
        x = Float64(row.x),
        y = Float64(row.y),
        z = Float64(row.z),
        shower = Int(getproperty(row, :shower)),
    )
end

function _worksteal_worker_loop(task_queue::RemoteChannel, corsika_binary::String, config::Dict)
    combined = DataFrame()
    have_combined = false

    while true
        task = take!(task_queue)
        if haskey(task, :stop) && task.stop
            break
        end

        particle = DataFrame(
            pdg=[task.pdg],
            energy=[task.energy],
            px=[task.px],
            py=[task.py],
            pz=[task.pz],
            x=[task.x],
            y=[task.y],
            z=[task.z],
            shower=[task.shower],
        )

        output = run_secondary_batch(particle, corsika_binary, config)
        if !isempty(output)
            if !have_combined
                combined = DataFrame(output, copycols=true)
                have_combined = true
            else
                append!(combined, output; promote=true)
            end
        end
        GC.gc()
    end

    return have_combined ? combined : DataFrame()
end

function run_secondaries_workstealing(
    secondaries::DataFrame,
    corsika_binary::String,
    config::Dict;
    nworkers::Int=nworkers()
)
    isempty(secondaries) && return DataFrame[]

    np = nworkers
    avail = workers()
    length(avail) < np && error("Requested $np workers but only $(length(avail)) Distributed workers are available")
    worker_pids = avail[1:np]

    tasks = [_secondary_task_from_row(row) for row in eachrow(secondaries)]
    queue = RemoteChannel(() -> Channel{NamedTuple}(length(tasks) + np))
    for t in tasks
        put!(queue, t)
    end
    for _ in 1:np
        put!(queue, (stop=true,))
    end

    futures = [@spawnat pid _worksteal_worker_loop(queue, corsika_binary, config) for pid in worker_pids]
    return map(fetch, futures)
end

"""
    distribute_and_run_secondaries(
        secondaries::DataFrame,
        corsika_binary::String,
        config::Dict,
        scheduler::Symbol=:naive;
        nworkers::Int=nprocs()-1
    )

Main orchestration function: distribute secondaries to workers and run in parallel.

# Arguments
- `secondaries`: DataFrame of secondary particles from first interaction
- `corsika_binary`: Path to modified CORSIKA binary
- `config`: Configuration dict (reuses some settings from primary run)
- `scheduler`: Scheduling strategy (:naive, :threshold, :binpack, :workstealing)
- `nworkers`: Number of parallel workers to use

# Returns
- `results`: Aggregated shower results from all secondaries
"""
function distribute_and_run_secondaries(
    secondaries::DataFrame,
    corsika_binary::String,
    config::Dict,
    scheduler::Symbol=:naive;
    nworkers::Int=nworkers()
)
    config = Dict{String,Any}(config)
    n_secondaries = nrow(secondaries)
    @info "Distributing $n_secondaries secondaries to $nworkers workers using :$scheduler scheduler"

    if scheduler == :workstealing
        results = run_secondaries_workstealing(secondaries, corsika_binary, config; nworkers=nworkers)
        @info "All secondary showers complete"
        return results
    end

    worker_batches = schedule_secondaries(secondaries, scheduler, nworkers)
    batch_sizes = [nrow(b) for b in worker_batches]
    @info "Created $(length(worker_batches)) batches: $batch_sizes particles"

    futures = []
    worker_pids = workers()[1:length(worker_batches)]
    for (i, batch) in enumerate(worker_batches)
        isempty(batch) && continue
        future = @spawnat worker_pids[i] run_secondary_batch(batch, corsika_binary, config)
        push!(futures, future)
    end

    results = map(fetch, futures)
    @info "All secondary showers complete"
    return results
end

"""
    run_secondary_batch(
        particles::DataFrame,
        corsika_binary::String,
        config::Dict
    )

Run CORSIKA for a batch of secondary particles (called on worker process).

# Arguments
- `particles`: DataFrame of secondary particles for this worker
- `corsika_binary`: Path to CORSIKA binary
- `config`: Configuration (inherited from primary run)

# Returns
- `output`: Aggregated particle output at observation plane
"""
function run_secondary_batch(
    particles::DataFrame,
    corsika_binary::String,
    config::Dict
)
    config = Dict{String,Any}(config)
    # Create temporary directory for this batch's output
    batch_dir = mktempdir()

    try
        # Persist each sub-shower output immediately, then append incrementally.
        results_dir = mkpath(joinpath(batch_dir, "subshower_outputs"))
        combined = DataFrame()
        have_combined = false

        for (i, row) in enumerate(eachrow(particles))
            energy = Float64(row.energy)
            px = Float64(row.px)
            py = Float64(row.py)
            pz = Float64(row.pz)

            momentum_norm = sqrt(px^2 + py^2 + pz^2)
            if !isfinite(momentum_norm) || momentum_norm == 0.0
                @warn "Skipping secondary with invalid momentum" px py pz
                continue
            end

            injection_pos = local_to_surface_coords_km(
                (Float64(row.x), Float64(row.y), Float64(row.z)),
                config,
            )
            direction = local_to_global_direction(
                (px / momentum_norm, py / momentum_norm, pz / momentum_norm),
                config,
            )
            plane_center = observation_plane_center_km(config)
            _, _, normal = observation_plane_basis(config)
            denom = direction[1] * normal[1] + direction[2] * normal[2] + direction[3] * normal[3]
            if !isfinite(denom) || abs(denom) < 1e-12
                @warn "Skipping secondary with direction parallel to observation plane" direction
                continue
            end

            delta_to_plane = (
                plane_center[1] - injection_pos[1],
                plane_center[2] - injection_pos[2],
                plane_center[3] - injection_pos[3],
            )
            t = (delta_to_plane[1] * normal[1] + delta_to_plane[2] * normal[2] + delta_to_plane[3] * normal[3]) / denom
            intercept = (
                injection_pos[1] + t * direction[1],
                injection_pos[2] + t * direction[2],
                injection_pos[3] + t * direction[3],
            )

            # Create batch-specific config
            secondary_config = Dict{String,Any}(config)
            secondary_config["pdg"] = Int(row.pdg)
            secondary_config["energy"] = energy
            secondary_config["xpos"] = injection_pos[1]
            secondary_config["ypos"] = injection_pos[2]
            secondary_config["zpos"] = injection_pos[3]
            secondary_config["xdir"] = normal[1]
            secondary_config["ydir"] = normal[2]
            secondary_config["zdir"] = normal[3]
            secondary_config["x-intercept"] = intercept[1]
            secondary_config["y-intercept"] = intercept[2]
            secondary_config["z-intercept"] = intercept[3]
            secondary_config["nevent"] = 1
            secondary_config["stop_after_first_interaction"] = false
            secondary_config["force_interaction"] = false
            secondary_config["output_tag"] = "subshower_$(getproperty(row, :shower))_$(i)"
            secondary_config["seed"] = derive_subshower_seed(get(config, "seed", 1), Int(getproperty(row, :shower)), i, Int(row.pdg))

            # Run CORSIKA for this secondary
            output_dir = run_corsika_with_stopping(
                corsika_binary,
                secondary_config,
                batch_dir;
                generation=1
            )

            # Sub-showers should emit observation-level particles.
            if isdir(output_dir)
                output = read_secondary_particles(output_dir; generation=1)
            else
                output = DataFrame()  # Empty if no particles reached observation level
            end

            output_file = joinpath(results_dir, string(secondary_config["output_tag"], ".parquet"))
            Parquet2.writefile(output_file, output)

            if !have_combined
                combined = DataFrame(output, copycols=true)
                have_combined = true
            else
                append!(combined, output; promote=true)
            end

            GC.gc()
        end
        return have_combined ? combined : DataFrame()

    catch err
        rethrow(err)
    finally
        # Clean up temporary directory
        rm(batch_dir, recursive=true, force=true)
    end
end
