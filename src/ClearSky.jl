module ClearSky

using Base: tail
using Base.Threads: @threads
using QuadGK: gauss
using Cubature
using BasicInterpolators
using ScalarRadau

#order matters
include("constants.jl")
include("util.jl")
include("molparam.jl")
include("par.jl")
include("faddeyeva.jl")
include("line_shapes.jl")
include("gases.jl")
include("collision_induced_absorption.jl")
include("atmospherics.jl")
include("radiation.jl")
include("modeling_core.jl")
include("modeling.jl")
include("orbits.jl")
include("insolation.jl")

using .Faddeyeva
export faddeyeva

end
