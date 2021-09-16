module ClearSky

using Base: tail
using Base.Threads: @threads, @spawn, Task, fetch
using BasicInterpolators
using Cubature: hquadrature
using Faddeyeva985
using QuadGK: gauss
using ScalarRadau
using ProgressMeter: Progress, next!

#order matters
include("constants.jl")
include("util.jl")
include("par.jl")
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

end
