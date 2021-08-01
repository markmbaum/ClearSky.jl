#-------------------------------------------------------------------------------
export opticaldepth

"""
    opticaldepth(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [optical depths](https://en.wikipedia.org/wiki/Optical_depth#Atmospheric_sciences) (``\\tau``) between two pressure levels

# Arguments

* `P₁`: first pressure level [Pa]
* `P₂`: second pressure level [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `θ`: angle [radians] of path, must be ∈ [0,π/2), where 0 is straight up/down
* `absorbers`: at least one gas object and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

Returns a vector of optical depths across all wavenumbers stored in gas objects. The `tol` keyword argument adjusts integrator error tolerance.
"""
function opticaldepth(P₁::Real,
                      P₂::Real,
                      g::Real,
                      fT::Q,
                      fμ::R,
                      θ::Real,
                      absorbers...;
                      tol::Float64=1e-5
                      )::Vector{Float64} where {Q,R}
    #initialization
    P₁, P₂ = max(P₁, P₂), min(P₁, P₂)
    A = unifyabsorbers(absorbers)
    checkpressures(A, P₁, P₂)
    checkazimuth(θ)
    #integrate wavenumbers in parallel
    τ = zeros(Float64, A.nν)
    m = 1/cos(θ)
    @threads for i ∈ eachindex(A.ν)
        #integrate
        τ[i] = depth(dτdω, P2ω(P₁), P2ω(P₂), A, i, g, m, fT, fμ, tol)
    end
    return τ
end

#-------------------------------------------------------------------------------
#the basic transmittance method exp(-τ) is already exported

"""
    transmittance(P₁, P₂, g, fT, fμ, θ, absorbers...; tol=1e-5)

Compute monochromatic [transmittances](https://en.wikipedia.org/wiki/Transmittance). between two pressure levels

Accepts the same arguments as [`opticaldepth`](@ref) and returns a vector of transmittances across all wavenumbers stored in gas objects.
"""
transmittance(args...; kwargs...) = transmittance.(opticaldepth(args...; kwargs...))

#-------------------------------------------------------------------------------
export outgoing

"""
    outgoing(Pₛ, g, fT, fμ, absorbers; Pₜ=1.0, nstream=5, tol=1e-5)

Compute outgoing monochromatic radiative fluxes [W/m``^2``/cm``^{-1}``], line-by-line. Integrates the [`schwarzschild`](@ref) equation from `Pₛ` to `Pₜ` at each wavenumber in the provided gas object(s) using any number of streams/angles. Total flux [W/m``^2``] can be evaluated with the [`trapz`](@ref) function. This function does not include reflected stellar radiation, which is accounted for in [`topfluxes`](@ref) and [`topimbalance`](@ref)

# Arguments
* `Pₛ`: surface pressure [Pa]
* `g`: gravitational acceleration [m/s``^2``]
* `fT`: temperature [K] as a function of pressure [Pa], `fT(P)`. This may be any callable object, like [`MoistAdiabat`](@ref), for example.
* `fμ`: mean molar mass as a function of temperature [K] and pressure [Pa], `fμ(T,P)`
* `absorbers`: at least one [gas object](gas_objects.md) and any number of [`CIATables`](@ref) and functions in the form σ(ν, T, P)

# Keywords
* `Pₜ`: top of atmopshere pressure [Pa]
* `nstream`: the number of atmospheric streams (radiation angles) to calculate and integrate
* `tol`: integrator error tolerance

The keyword argument `nstream` specifies how many independent streams, or beam angles through the atmosphere, to integrate. The keyword argument `tol` is a numerical error tolerance passed to the [`radau`](https://github.com/markmbaum/ScalarRadau.jl) integrator.
"""
function outgoing(Pₛ::Real,
                  g::Real,
                  fT::Q,
                  fμ::R,
                  absorbers...;
                  Pₜ::Float64=1.0, #top of atmosphere pressure, 1 Pa by default
                  nstream::Int64=5, #number of streams to use
                  tol::Float64=1e-5 #integrator tolerance
                  )::Vector{Float64} where {Q,R}
    #initialization
    A = unifyabsorbers(absorbers)
    checkpressures(A, Pₛ, Pₜ)
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    #surface temperature
    T₀ = fT(Pₛ)
    #integrate wavenumbers in parallel
    F = zeros(Float64, A.nν)
    @threads for i ∈ eachindex(A.ν)
        I₀ = planck(A.ν[i], T₀)
        F[i] = streams(dIdω, I₀, ω₁, ω₂, A, i, g, nstream, fT, fμ, tol)
    end
    return F
end

#-------------------------------------------------------------------------------
export topflux, surfaceflux

function topflux(Pₛ::Real,
                 g::Real,
                 fT::Q, # fT(P)
                 fμ::R, # fμ(T,P)
                 fS::S, # fS(ν) [W/m^2]
                 fα::U, #albedo as function of wavenumber
                 absorbers...;
                 Pₜ::Float64=1.0, #top of atmosphere pressure, 1 Pa by default
                 nstream::Int64=5, #number of streams to use
                 θₛ::Float64=0.841, #corresponds to cos(θ) = 2/3
                 tol::Float64=1e-4
                 )::NTuple{2,Vector{Float64}} where {Q,R,S,U}
    #setup
    A = unifyabsorbers(absorbers)
    checkpressures(A, Pₛ, Pₜ)
    checkazimuth(θₛ)
    ν, nν = A.ν, A.nν
    ω₁, ω₂ = P2ω(Pₛ, Pₜ)
    ι₁, ι₂ = P2ι(Pₜ, Pₛ)

    Tₛ = fT(Pₛ) #surface temperature
    Fₜ⁻ = fS.(ν)*cos(θₛ) #incoming stellar flux at TOA
    Fₛ⁻ = zeros(nν) #downward stellar flux at the ground/surface
    Fₜ⁺ = zeros(nν) #outgoing flux at TOA
    m = 1/cos(θₛ)

    @threads for i ∈ eachindex(ν)
        #downward stellar flux at surface
        τ = depth(dτdι, ι₁, ι₂, A, i, g, m, fT, fμ, tol)
        Fₛ⁻[i] = Fₜ⁻[i]*exp(-τ)
        #some of it gets reflected back upward at the surface
        Iₛ⁺ = Fₛ⁻[i]*fα(ν[i])/π #Lambertian reflection
        #then radiation streams out of the atmosphere
        I₀ = Iₛ⁺ + planck(ν[i], Tₛ)
        Fₜ⁺[i] = streams(dIdω, I₀, ω₁, ω₂, A, i, g, nstream, fT, fμ, tol)
    end

    return Fₜ⁻, Fₜ⁺
end

function surfaceflux(Pₛ::Real,
                     g::Real,
                     fT::Q, # fT(P)
                     fμ::R, # fμ(T,P)
                     fS::S, # fS(ν) [W/m^2]
                     fα::U, #albedo as function of wavenumber
                     absorbers...;
                     Pₜ::Float64=1.0, #top of atmosphere pressure, 1 Pa by default
                     nstream::Int64=5, #number of streams to use
                     θₛ::Float64=0.841, #corresponds to cos(θ) = 2/3
                     tol::Float64=1e-4
                     )::NTuple{2,Vector{Float64}} where {Q,R,S,U}
    #setup
    A = unifyabsorbers(absorbers)
    checkpressures(A, Pₛ, Pₜ)
    checkazimuth(θₛ)
    ν, nν = A.ν, A.nν
    ι₁, ι₂ = P2ι(Pₜ, Pₛ)

    Tₛ = fT(Pₛ) #surface temperature
    Fₜ⁻ = fS.(ν)*cos(θₛ) #incoming stellar flux at TOA
    Fₛ⁻ = zeros(nν) #downward flux at the ground/surface
    Fₛ⁺ = zeros(nν) #upward flux at surface
    m = 1/cos(θₛ)
    
    @threads for i ∈ eachindex(ν)
        #incoming stellar flux at surface
        τ = depth(dτdι, ι₁, ι₂, A, i, g, m, fT, fμ, tol)
        Fₛ⁻[i] += Fₜ⁻[i]*exp(-τ)
        #some of it gets reflected back upward at the surface
        Fₛ⁺[i] = Fₛ⁻[i]*fα(ν[i])
        #and the surface emits some radiation
        Fₛ⁺[i] += π*planck(ν[i], Tₛ)
        #the atmosphere emits downward too
        Fₛ⁻[i] += streams(dIdι, 0.0, ι₁, ι₂, A, i, g, nstream, fT, fμ, tol)
    end

    return Fₛ⁻, Fₛ⁺
end

#-------------------------------------------------------------------------------
export netflux

function netflux(P::AbstractVector{<:Real},
                 g::Real,
                 fT::Q,
                 fμ::R,
                 fS::S,
                 fα::U,
                 absorbers...;
                 nstream::Int64=5,
                 θₛ::Float64=0.841,
                 tol::Float64=1e-5
                 )::Vector{Float64} where {Q,R,S,U}
    #setup
    A = unifyabsorbers(absorbers)
    @assert P[end] < P[1] "pressure coordinates must be in descending order"
    checkpressures(A, P[1], P[end])
    P = reverse(collect(Float64, P))
    nP = length(P)
    checkazimuth(θₛ)
    ν, nν = A.ν, A.nν
    ω = reverse(P2ω.(P))
    ι = P2ι.(P)

    #number of threads
    nthr = nthreads()
    #temporary storage
    X = zeros(Float64, nP, nthr)
    #surface temperature
    Tₛ = fT(P[1]) #surface temperature
    #big blocks of flux
    M⁺ = zeros(Float64, nP, nν) #monochromatic upward fluxes
    M⁻ = zeros(Float64, nP, nν) #monochromatic downward fluxes
    #F⁺ = zeros(Float64, nP) #total updward fluxes
    #F⁻ = zeros(Float64, nP) #total downward fluxes
    Fₙ = zeros(Float64, nP) #net fluxes
    
    @threads for i ∈ eachindex(ν)
        #temporary storage
        x = view(X, :, threadid())
        #downward stellar irradiance
        I₀ = fS(ν[i])
        #downward stellar irradiance through atmosphere
        M = view(M⁻, :, i)
        stream!(dIdι, I₀, M, ι, ι[1], ι[end], A, i, g, 1/cos(θₛ), fT, fμ, tol)
        M .*= cos(θₛ) #convert to flux
        #some of the stellar flux is reflected
        Iₛ⁺ = M[end]*fα(ν[i])/π #Lambertian reflection
        #and the surface emits radiation
        I₀ = Iₛ⁺ + planck(ν[i], Tₛ)
        #total upward radiation streams
        M = view(M⁺, nP:-1:1, i) #reversed view, from surface to TOA
        streams!(dIdω, I₀, x, M, ω, ω[1], ω[end], A, i, g, nstream, fT, fμ, tol)
        #downward atmospheric contribution
        M = view(M⁻, :, i)
        streams!(dIdι, 0.0, x, M, ι, ι[1], ι[end], A, i, g, nstream, fT, fμ, tol)
    end

    @threads for i ∈ eachindex(P)
        Fₙ[i] = trapz(ν, view(M⁺, i, :)) - trapz(ν, view(M⁻, i, :))
    end

    return reverse(Fₙ)

end