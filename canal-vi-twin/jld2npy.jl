# dump JLD2Writer timeseries (u, v, c, t) to raw float64 files: julia jld2npy.jl <fields.jld2> <out prefix>
using JLD2
f = jldopen(ARGS[1]); its = sort(parse.(Int, keys(f["timeseries/t"])))
for n in ("u", "v", "c")
    A = cat((f["timeseries/$n/$k"][:, :, 4] for k in its)...; dims=3)
    write(ARGS[2] * "_$n.bin", A); println(n, " ", size(A))
end
write(ARGS[2] * "_t.bin", [f["timeseries/t/$k"] for k in its]); println("t ", length(its))
