#-------------------------------------------------------------------------------
# general checks and stuff

function checkazimuth(θ)
    @assert 0 <= θ < π/2 "azimuth angle θ must be ∈ [0,π/2)"
end

function checkstreams(n)
    n < 4 && @warn "careful! using nstream < 4 is likely to be inaccurate!" maxlog=1 
end

#handling of inputs for temperature, molar mass, heat capacity
formprofile(P, x::AbstractVector) = AtmosphericProfile(P, x)
formprofile(::Any, x::Real) = (::Any...) -> x
formprofile(::Any, x) = x
formprofiles(P, X...) = map(x -> formprofile(P,x), X)

#------------------------------------------------------------------------------
export opticaldepth

"""
    opticaldepth(P₁, P₂, g, 𝒻T, 𝒻μ, θ, absorbers...; tol=1e-5)

Compute monochromatic [optical depths](https://en.wikipedia.org/wiki/Optical_depth#Atmospheric_sciences) (``\\tau``) between two pressure levels

# Arguments

* `P₁`: first pressure level [Pa]
* `P₂`: second pressure level [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `𝒻T`: temperature [K] as a function of pressure [Pa], `𝒻T(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `𝒻μ`: mean molar mass as a function of temperature [K] and pressure [Pa], `𝒻μ(T,P)`
* `θ`: angle [radians] of path, must be ∈ [0,π/2), where 0 is straight up/down
* `absorbers`: at least one gas object and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

Returns a vector of optical depths across all wavenumbers stored in gas objects. The `tol` keyword argument adjusts integrator error tolerance.
"""
function opticaldepth(P₁::Real,
                      P₂::Real,
                      g::Real,
                      𝒻T::Q,
                      𝒻μ::R,
                      θ::Real,
                      absorbers...;
                      tol::Float64=1e-5) where {Q,R}
    #parse absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #get pressure limits from high to low
    P₁, P₂ = max(P₁, P₂), min(P₁, P₂)
    #converted pressure coordinates
    ω₁, ω₂ = P2ω(P₁), P2ω(P₂)
    #zenith angle factor    
    𝓂 = 1/cos(θ)
    #final checks
    checkpressures(𝒜, P₁, P₂)
    checkazimuth(θ)

    #spawn asynchronous parallel integrations using radau function
    tasks = Vector{Task}(undef, nν)
    for i ∈ eachindex(ν)
        tasks[i] = @spawn 𝓇depth(dτdω, ω₁, ω₂, 𝒜, i, g, 𝓂, 𝒻T, 𝒻μ, tol)
    end
    #fetch and return the results
    [fetch(task) for task ∈ tasks]
end

function opticaldepth(P::AbstractVector{<:Real},
                      g::Real,
                      T::Q,
                      μ::R,
                      θ::Real,
                      absorbers...;
                      nlobatto::Int=4) where {Q,R}
    #parse absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #parse profiles
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #let pressure be ascending, shouldn't matter here
    P = sort(P)
    #angle factor 
    𝓂 = 1/cos(θ)
    #precalculated constant
    C = 1e-4*𝐍𝐚/g
    #pre-evaluate T and μ at intra-layer lobatto nodes
    T, μ = lobattoevaluations(P, 𝒻T, 𝒻μ, nlobatto)
    #final checks
    checkpressures(𝒜, P[end], P[1])
    checkazimuth(θ)
    
    #integrate in parallel
    τ = Vector{eltype(T)}(undef,nν)
    @threads for i ∈ eachindex(ν)
        τ[i] = 𝒹depth(P, T, μ, 𝒜, i, C, 𝓂, nlobatto)
    end
    return τ
end

#------------------------------------------------------------------------------
#already exported

"""
    transmittance(P₁, P₂, g, 𝒻T, 𝒻μ, θ, absorbers...; tol=1e-5)

Compute monochromatic [transmittances](https://en.wikipedia.org/wiki/Transmittance). between two pressure levels

Accepts the same arguments as [`opticaldepth`](@ref) and returns a vector of transmittances across all wavenumbers stored in gas objects.
"""
transmittance(args...; kwargs...) = exp.(-opticaldepth(args...; kwargs...))

#------------------------------------------------------------------------------
export outgoing

"""
    outgoing(Pₛ, g, 𝒻T, 𝒻μ, absorbers; Ptop=1.0, nstream=5, tol=1e-5)

Compute outgoing monochromatic radiative fluxes [W/m``^2``/cm``^{-1}``], line-by-line. Integrates the [`schwarzschild`](@ref) equation from `Pₛ` to `Ptop` at each wavenumber in the provided gas object(s) using any number of streams/angles. Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function. This function does not include reflected stellar radiation, which is accounted for in [`topfluxes`](@ref) and [`topimbalance`](@ref)

# Arguments
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `𝒻T`: temperature [K] as a function of pressure [Pa], `𝒻T(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `𝒻μ`: mean molar mass as a function of temperature [K] and pressure [Pa], `𝒻μ(T,P)`
* `absorbers`: at least one [gas object](gas_objects.md) and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

# Keywords
* `Ptop`: top of atmopshere pressure [Pa]
* `nstream`: the number of atmospheric streams (radiation angles) to calculate and integrate
* `tol`: integrator error tolerance

The keyword argument `nstream` specifies how many independent streams, or beam angles through the atmosphere, to integrate. The keyword argument `tol` is a numerical error tolerance passed to the [`radau`](https://github.com/markmbaum/ScalarRadau.jl) integrator.
"""
function outgoing(Pₛ::Real,
                  g::Real,
                  𝒻T::Q,
                  𝒻μ::R,
                  absorbers...;
                  Ptop::Real=1.0, #top of atmosphere pressure, 1 Pa by default
                  nstream::Int=5, #number of streams to use
                  tol::Real=1e-5) where {Q,R}
    #parse the absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #convert pressure coordinates
    ω₁, ω₂ = P2ω(Pₛ, Ptop)
    #surface temperature
    Tₛ = 𝒻T(Pₛ)
    #final checks
    checkpressures(𝒜, Pₛ, Ptop)
    checkstreams(nstream)

    #integrate in parallel, dynamic schedule
    tasks = Vector{Task}(undef, nν)
    for (i,νᵢ) ∈ enumerate(ν)
        I₀ = planck(νᵢ, Tₛ)
        tasks[i] = @spawn 𝓇streams(dIdω, I₀, ω₁, ω₂, 𝒜, i, g, 𝒻T, 𝒻μ, nstream, tol)
    end
    [fetch(task) for task ∈ tasks]
end

function outgoing(P::AbstractVector{<:Real},
                  g::Real,
                  T::Q,
                  μ::R,
                  absorbers...;
                  nstream::Int=5,
                  nlobatto::Int=3) where {Q,R}
    #parse the absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #parse profiles
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #pressure should be descending for outgoing calculation
    P = sort(P, rev=true)
    #surface temperature [K]
    Tₛ = 𝒻T(P[1])
    #constant
    C = 1e-4*𝐍𝐚/g
    #pre-evaluate T and μ at intra-layer lobatto nodes
    T, μ = lobattoevaluations(P, 𝒻T, 𝒻μ, nlobatto)
    #final checks
    checkpressures(𝒜, P[1], P[end])
    checkstreams(nstream)

    #no need for asynchronous when fully-discretized
    olr = zeros(typeof(Tₛ), nν)
    @threads for i ∈ eachindex(ν)
        #initial irradiance
        I₀ = planck(ν[i], Tₛ)
        #monochromatic outgoing flux
        olr[i] = 𝒹streams(P, T, μ, 𝒜, i, I₀, C, nstream, nlobatto)
    end
    return olr
end

#------------------------------------------------------------------------------
export monochromaticfluxes

function monochromaticfluxes!(M⁺::AbstractMatrix,
                              M⁻::AbstractMatrix,
                              τ::AbstractMatrix,
                              core::Radau,
                              P::AbstractVector{<:Real},
                              g::Real,
                              T::Q,
                              μ::R,
                              𝒻S::U,
                              𝒻a::V,
                              absorbers...;
                              θₛ::Real=0.841)::Nothing where {Q,R,U,V}
    #parse absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #parse profiles
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #integration instructions
    @unpack nstream, tol = core
    #check ascending order of pressure coordinates
    @assert issorted(P) "pressure coordinates must be in ascending order (sorted)"
    #transformed coordinates
    ω, ι = P2ω.(P), P2ι.(P)
    reverse!(ω)
    #final checks
    checkpressures(𝒜, P[end], P[1])
    checkstreams(nstream)
    checkazimuth(θₛ)
    #make it clear that optical depth is not used
    τ[:] .*= NaN

    #asynchronous, parallel integrations
    tasks = Vector{Task}(undef, nν)
    for j ∈ eachindex(ν)
        Mⱼ⁺ = @view M⁺[:,j]
        Mⱼ⁻ = @view M⁻[:,j]
        tasks[j] = @spawn 𝓇monoflux!(Mⱼ⁺, Mⱼ⁻, P, ω, ι, 𝒜, j, g, 𝒻T, 𝒻μ, 𝒻S, 𝒻a, θₛ, nstream, tol)
    end
    [fetch(task) for task ∈ tasks]
    return nothing
end

function monochromaticfluxes!(M⁺::AbstractMatrix,
                              M⁻::AbstractMatrix,
                              τ::AbstractMatrix,
                              core::Discretized,
                              P::AbstractVector{<:Real},
                              g::Real,
                              T::Q,
                              μ::R,
                              𝒻S::U,
                              𝒻a::V,
                              absorbers...;
                              θₛ::Real=0.841)::Nothing where {Q,R,U,V}
    #parse absorber(s)
    𝒜, ν, _ = unifyabsorbers(absorbers)
    #parse profiles
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #integration instructions
    @unpack nstream, nlobatto = core
    #check ascending order of pressure coordinates
    @assert issorted(P) "pressure coordinates must be in ascending order (sorted)"
    #precalculated constant
    C = 1e-4*𝐍𝐚/g
    #temperature and molar mass at lobatto nodes within layers
    T, μ = lobattoevaluations(P, 𝒻T, 𝒻μ, nlobatto)
    #pre-evaluate Planck at pressure coordinates
    B = planckevaluations(P, ν, 𝒻T)
    #final checks
    checkpressures(𝒜, P[end], P[1])
    checkstreams(nstream)
    checkazimuth(θₛ)

    #no need for asynchronous when fully discretized
    @threads for j ∈ eachindex(ν)
        τⱼ = @view τ[:,j]
        𝒹depth!(τⱼ, P, T, μ, 𝒜, j, C, nlobatto)
        Mⱼ⁺ = @view M⁺[:,j]
        Mⱼ⁻ = @view M⁻[:,j]
        Bⱼ = @view B[:,j]
        𝒹monoflux!(Mⱼ⁺, Mⱼ⁻, τⱼ, P, Bⱼ, ν[j], 𝒻S, 𝒻a, θₛ, nstream)
    end
    return nothing
end

function monochromaticfluxes(P::AbstractVector{<:Real},
                             g::Real,
                             T::Q,
                             μ::R,
                             𝒻S::U,
                             𝒻a::V,
                             absorbers...;
                             core::AbstractNumericalCore=Discretized(),
                             θₛ::Real=0.841) where {Q,R,U,V}
    #parse absorber(s)
    𝒜, _, nν = unifyabsorbers(absorbers)
    #parse profile
    𝒻T = formprofile(P, T)
    #output type
    W = typeof(𝒻T(P[1]))
    #big blocks of flux
    np = length(P)
    M⁺ = zeros(W, np, nν) #monochromatic upward fluxes
    M⁻ = zeros(W, np, nν) #monochromatic downward
    τ = zeros(W, np-1, nν)

    #main calculation of all fluxes
    monochromaticfluxes!(M⁺, M⁻, τ, core, P, g, T, μ, 𝒻S, 𝒻a, 𝒜; θₛ=θₛ)

    return M⁺, M⁻
end

#------------------------------------------------------------------------------
export fluxes, netfluxes

function fluxes(P::AbstractVector{<:Real},
                g::Real,
                T::Q,
                μ::R,
                𝒻S::U,
                𝒻a::V,
                absorbers...;
                core::AbstractNumericalCore=Discretized(),
                θₛ::Real=0.841) where {Q,R,U,V} θₛ
    #parse absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #parse profiles
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #coordinate count
    np = length(P)
    #output type
    W = typeof(𝒻T(P[1]))

    #big blocks of flux and depth
    M⁺ = zeros(W, np, nν) #monochromatic upward fluxes
    M⁻ = zeros(W, np, nν) #monochromatic downward
    τ = zeros(W, np-1, nν)
    #big blocks of monochromatic fluxes
    monochromaticfluxes!(M⁺, M⁻, τ, core, P, g, 𝒻T, 𝒻μ, 𝒻S, 𝒻a, 𝒜; θₛ=θₛ)
    #integrate across wavenumbers
    F⁺ = similar(M⁺, np)
    F⁻ = similar(M⁻, np)
    ∫F!(F⁺, F⁻, M⁺, M⁻, ν)
    return F⁺, F⁻
end

function netfluxes(P::AbstractVector{<:Real},
                   g::Real,
                   T::Q,
                   μ::R,
                   𝒻S::U,
                   𝒻a::V,
                   absorbers...;
                   kwargs...) where {Q,R,U,V}
    F⁺, F⁻ = fluxes(P, g, T, μ, 𝒻S, 𝒻a, absorbers...; kwargs...)
    return F⁺ .- F⁻
end

#------------------------------------------------------------------------------
export radiate!, radiate

function radiate!(F::FluxPack,
                  core::AbstractNumericalCore,
                  P::AbstractVector{<:Real},
                  g::Real,
                  T::Q,
                  μ::R,
                  𝒻S::U,
                  𝒻a::V,
                  absorbers...;
                  kwargs...)::Nothing where {Q,R,U,V}
    #parse absorber(s)
    𝒜, ν, nν = unifyabsorbers(absorbers)
    #form profile functions
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #unpack the model arrays
    @unpack τ, M⁺, M⁻, F⁺, F⁻, Fnet = F
    #check sizes
    np = length(P)
    @assert size(F) == (np, nν) "size of FluxPack does not match number of pressure or wavenumber coordinates"
    #compute monochromatic fluxes in place
    monochromaticfluxes!(M⁺, M⁻, τ, core, P, g, 𝒻T, 𝒻μ, 𝒻S, 𝒻a, 𝒜; kwargs...)
    #the other stuff is easy
    ∫F!(F⁺, F⁻, M⁺, M⁻, ν)
    @. Fnet = F⁺ - F⁻
    #nothing to return
    return nothing
end

function radiate(P::AbstractVector{<:Real},
                 g::Real,
                 T::Q,
                 μ::R,
                 𝒻S::U,
                 𝒻a::V,
                 absorbers...;
                 core::AbstractNumericalCore=Discretized(),
                 θₛ::Real=0.841
                 )::FluxPack where {Q,R,U,V}
    #parse absorber(s)
    𝒜, ν, _ = unifyabsorbers(absorbers)
    #form profile functions
    𝒻T, 𝒻μ = formprofiles(P, T, μ)
    #allocate arrays of the proper type
    F = FluxPack(P, ν, typeof(𝒻T(P[1])))
    #compute
    radiate!(F, core, P, g, T, 𝒻μ, 𝒻S, 𝒻a, 𝒜; θₛ=θₛ)
    return F
end