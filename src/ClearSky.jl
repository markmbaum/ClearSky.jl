module ClearSky

using Base: tail
using Base.Threads: @threads, @spawn, Task, fetch
using BasicInterpolators
using Cubature: hquadrature
using Faddeyeva985
using FastGaussQuadrature: gausslegendre, gausslobatto
using Memoize: @memoize
using ProgressMeter: Progress, next!
using ScalarRadau
using UnPack: @unpack

#order shouldn't matter for these
include("constants.jl")
include("util.jl")
include("radiation.jl")
include("orbits.jl")
include("insolation.jl")

#reading HITRAN data and molecular properties
include(joinpath("hitran", "par.jl"))
include(joinpath("hitran", "molparam.jl"))

#computing absorption cross-sections & coefficients
include(joinpath("absorption", "line_shapes.jl"))
include(joinpath("absorption", "gases.jl"))
include(joinpath("absorption", "collision_induced_absorption.jl"))
include(joinpath("absorption", "absorbers.jl"))

#adiabatic profiles, lapse rates, hydrostatics, etc.
include("atmospherics.jl")

#core modeling functions/numerics
include(joinpath("core", "shared.jl"))
include(joinpath("core", "radau.jl"))
include(joinpath("core", "discretized.jl"))

#wrappers of core functions for one shot calculations
include("fluxes.jl")
#highest-level modeling stuff
include("radiative_convective.jl")

end
