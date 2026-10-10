# Oceananigans twin of the NEMO CANAL barotropic jet with NEMO's exact geometry and initial noise:
#   julia canal_twin5.jl <Δt seconds> <scheme: weno5|weno3|een (vorticity only) | wvi5|wvi3 (full WENOVectorInvariant)> <output prefix>
#   env: SUBSTEPS (default 90), STOPDAYS (default 60), OUTHOURS (default 6), NOISE (default canal_noise.bin next to this script)
using Oceananigans
using Oceananigans.Units
using Printf, JLD2

Δt_seconds = parse(Float64, ARGS[1])
vort = ARGS[2]
prefix = ARGS[3]

# NEMO CANAL (rn_dx = 10 km): 201 periodic T columns centred on x = -1000…1000 km,
# 101 T rows centred on y = -500…500 km of which the first and last are land -> 99 wet rows, walls at y = ±495 km
Nx, Ny, Δ, H = 201, 99, 10kilometers, 1000
grid = RectilinearGrid(CPU(), size=(Nx, Ny, 1), x=(-Nx*Δ/2, Nx*Δ/2), y=(-Ny*Δ/2, Ny*Δ/2), z=(-H, 0),
                       halo=(5, 5, 3), topology=(Periodic, Bounded, Bounded))

g  = 9.80665
f₀ = 2 * 7.292115e-5 * sind(38.5)
if startswith(vort, "wvi")
    momentum_advection = WENOVectorInvariant(order=parse(Int, vort[4:end]))
else
    vorticity_scheme = vort == "weno5" ? WENO(order=5) : vort == "weno3" ? WENO(order=3) : Oceananigans.Advection.EnergyConserving()
    momentum_advection = VectorInvariant(; vorticity_scheme)
end
println(momentum_advection)

model = HydrostaticFreeSurfaceModel(grid; momentum_advection,
                                    coriolis = FPlane(f=f₀),
                                    free_surface = SplitExplicitFreeSurface(grid; substeps=parse(Int, get(ENV, "SUBSTEPS", "90")), gravitational_acceleration=g),
                                    tracers = :c, buoyancy = nothing,
                                    tracer_advection = WENO(order=5),
                                    timestepper = :SplitRungeKutta3)

# NEMO CANAL nn_initcase = 1: top-hat geostrophic jet, dye 35 inside / 30 outside (T-point test |y| <= 200 km)
Uj, yj = 1.0, 200kilometers + 1     # +1 m: NEMO's T rows at y = ±200 km are inside the jet
uᵢ(x, y, z) = abs(y) <= yj ? Uj : 0.0
ηᵢ(x, y, z) = abs(y) <= yj ? -f₀ * Uj * y / g : -f₀ * Uj * sign(y) * 200kilometers / g
cᵢ(x, y, z) = abs(y) <= yj ? 35.0 : 30.0
set!(model, u=uᵢ, η=ηᵢ, c=cᵢ)

# ln_sshnoise: pv(ji,jj) += 0.5 (zrandom(ji,jj) - 0.5), zrandom(201, 101) from gfortran RANDOM_NUMBER with seed 1234567.
# NEMO global V row jj sits at y = -495 + 10 (jj-1) km, the same as Oceananigans v face j = jj; faces 2…99 are wet.
noise = Array{Float64}(undef, Nx, Ny + 2)
read!(get(ENV, "NOISE", joinpath(@__DIR__, "canal_noise.bin")), noise)
interior(model.velocities.v, :, 2:Ny, 1) .+= 0.5 .* (noise[:, 2:Ny] .- 0.5)

u, v, w = model.velocities
ζ = Field(∂x(v) - ∂y(u))
simulation = Simulation(model; Δt=Δt_seconds, stop_time=parse(Float64, get(ENV, "STOPDAYS", "60")) * days)

times, KE, Z, umax = Float64[], Float64[], Float64[], Float64[]
function record(sim)
    compute!(ζ)
    push!(times, time(sim))
    push!(KE, 0.5 * (sum(interior(u).^2) + sum(interior(v).^2)))
    push!(Z, 0.5 * sum(interior(ζ).^2))
    push!(umax, max(maximum(abs, interior(u)), maximum(abs, interior(v))))
    @printf("day %5.2f  max|U| %.3f  KE/KE₀ %.3f  Z/Z₀ %.3f\n", time(sim)/day, umax[end], KE[end]/KE[1], Z[end]/Z[1])
    isfinite(umax[end]) && umax[end] < 10 || (sim.running = false)
end
out = TimeInterval(parse(Float64, get(ENV, "OUTHOURS", "6")) * hours)
add_callback!(simulation, record, out)
simulation.output_writers[:fields] = JLD2Writer(model, (; ζ, c=model.tracers.c, u, v); filename=prefix * "_fields.jld2",
                                                 schedule=out, overwrite_files=true)
run!(simulation)
jldsave(prefix * "_stats.jld2"; times, KE, Z, umax)
