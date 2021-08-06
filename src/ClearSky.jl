module ClearSky

using Base: tail
using Base.Threads: @threads, nthreads, threadid
using BasicInterpolators
using Cubature: hquadrature
using QuadGK: gauss
using ScalarRadau
using ProgressMeter: Progress, next!

#order matters
include("constants.jl")
include("util.jl")
include("par.jl")
include("faddeyeva.jl")
include("line_shapes.jl")
include("molparam.jl")
include("gases.jl")
include("collision_induced_absorption.jl")
include("atmospherics.jl")
include("radiation.jl")
include("modeling_core.jl")
include("modeling.jl")
include("orbits.jl")
include("insolation.jl")

using .Faddeyeva

end
