#-------------------------------------------------------------------------------
# core differential equations

function dτdP(P, ::Any, param::Tuple)
    #unpack parameters
    𝒜, idx, g, 𝓂, fT, fμ = param
    #temperature from given profile [K]
    T = fT(P)
    #mean molar mass [mole/kg]
    μ = fμ(T, P)
    #sum of all cross-sections [cm^2/molecule]
    σ = Σ(𝒜, idx, T, P)
    #dτ/dlnP, scaled by the angle m = 1/cos(θ)
    𝓂*dτdP(σ, g, μ) #no Planck emission
end

function dIdP(P, I, param::Tuple)
    #unpack parameters
    𝒜, idx, g, 𝓂, fT, fμ = param
    #temperature from given profile [K]
    T = fT(P)
    #mean molar mass [mole/kg]
    μ = fμ(T, P)
    #sum of all cross-sections [cm^2/molecule]
    σ = Σ(𝒜, idx, T, P)
    #pull out wavenumber
    ν = @inbounds 𝒜.ν[idx]
    #dI/dlnP, scaled by the angle m = 1/cos(θ)
    𝓂*schwarzschild(I, ν, σ, g, μ, T)
end

function dJdP(P, I, param::Tuple)
    #unpack parameters
    𝒜, idx, g, 𝓂, fT, fμ = param
    #temperature from given profile [K]
    T = fT(P)
    #mean molar mass [mole/kg]
    μ = fμ(T, P)
    #sum of all cross-sections [cm^2/molecule]
    σ = Σ(𝒜, idx, T, P)
    #dI/dlnP, scaled by the angle m = 1/cos(θ), without emission
    𝓂*absorption(I, σ, g, μ)
end

#-------------------------------------------------------------------------------
# wrappers for log pressure coordinates

function dτdι(ι, τ, param::Tuple)
    P = ι2P(ι)
    dιfac(P)*dτdP(P, τ, param)
end

function dτdω(ω, τ, param::Tuple)
    P = ω2P(ω)
    dωfac(P)*dτdP(P, τ, param)
end

function dIdω(ω, I, param::Tuple)
    P = ω2P(ω)
    dωfac(P)*dIdP(P, I, param)
end

function dIdι(ι, I, param::Tuple)
    P = ι2P(ι)
    dιfac(P)*dIdP(P, I, param)
end

function dJdι(ι, I, param::Tuple)
    P = ι2P(ι)
    dιfac(P)*dJdP(P, I, param)
end

#-------------------------------------------------------------------------------
# functions for optical depth paths

function 𝓇depth(dτdx::Q,
                x₁::Real,
                x₂::Real,
                𝒜::R,
                idx::Int,
                g::Real,
                𝓂::Real, # 1/cos(θ)
                fT::S,
                fμ::U,
                tol::Float64
                ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters, copying the functions to avoid race conditions
    param = (𝒜, idx, g, 𝓂, copyprofile(fT), copyprofile(fμ))
    #integrate with the ODE solver (appears to be faster than quadrature)
    radau(dτdx, zero(x₁), x₁, x₂, param, atol=tol, rtol=tol)
end

function 𝓇depth!(τ,
                 x,
                 dτdx::Q,
                 x₁::Real,
                 x₂::Real,
                 𝒜::R,
                 idx::Int,
                 g::Real,
                 𝓂::Real, # 1/cos(θ)
                 fT::S,
                 fμ::U,
                 tol::Float64
                 ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters, copying the functions to avoid race conditions
    param = (𝒜, idx, g, 𝓂, copyprofile(fT), copyprofile(fμ))
    #integrate with the ODE solver (appears to be faster than quadrature)
    radau!(τ, x, dτdx, zero(eltype(τ)), x₁, x₂, param, atol=tol, rtol=tol)
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere, no storage

function 𝓇stream(dIdx::Q, #version of schwarzschild equation
                 I₀::Real, #initial irradiance
                 x₁::Real, #initial pressure coordinate
                 x₂::Real, #final pressure coordinate
                 𝒜::R,
                 idx::Int,
                 g::Real, #gravity [m/s^2]
                 𝓂::Real, #1/cos(θ), where θ is the stream angle
                 fT::S, #temperature profile fT(P)
                 fμ::U, #mean molar mass μ(T,P)
                 tol::Real #integrator error tolerance
                 ) where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (𝒜, idx, g, 𝓂, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords and return
    radau(dIdx, I₀, x₁, x₂, param, atol=tol, rtol=tol)
end

function 𝓇streams(dIdx::Q, #version of schwarzschild equation
                  I₀::Real, #initial irradiance
                  x₁::Real, #initial pressure coordinate
                  x₂::Real, #final pressure coordinate
                  𝒜::R,
                  idx::Int, #index of wavenumber
                  g::Real, #gravity [m/s^2]
                  fT::S, #temperature profile fT(P)
                  fμ::U, #mean molar mass μ(T,P)
                  nstream::Int,
                  tol::Real #integrator error tolerance
                  ) where {Q,R<:AbstractAbsorber,S,U}
    #setup gaussian quadrature nodes
    𝓂, 𝒲 = streamnodes(nstream)
    #copy profiles to avoid interpolator race conditions
    fT, fμ = copyprofiles(fT, fμ)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    M = zero(I₀)
    for i ∈ 1:nstream
        I = 𝓇stream(dIdx, I₀, x₁, x₂, 𝒜, idx, g, 𝓂[i], fT, fμ, tol)
        # integral over hemisphere: ∫∫ I cos(θ) sin(θ) dθ dϕ, where θ∈[0,π/2], ϕ∈[0,2π]
        M += 𝒲[i]*I #𝒲 = 2π*w*cos(θ)*sin(θ), precomputed
    end
    return M
end

#-------------------------------------------------------------------------------
# functions for streams and fluxes up/down the atmosphere with in-place storage

function 𝓇stream!(I, #output/solution vector
                  x, #output/solution coordinates
                  dIdx::Q, #version of schwarzschild equation
                  I₀::Real, #initial irradiance
                  𝒜::R,
                  idx::Int,
                  g::Real, #gravity [m/s^2]
                  𝓂::Real, #1/cos(θ), where θ is the stream angle
                  fT::S, #temperature profile fT(P)
                  fμ::U, #mean molar mass μ(T,P)
                  tol::Real #integrator error tolerance
                  )::Nothing where {Q,R<:AbstractAbsorber,S,U}
    #pack parameters
    param = (𝒜, idx, g, 𝓂, fT, fμ)
    #integrate the Schwarzschild equation in log pressure coords, in-place
    radau!(I, x, dIdx, I₀, x[1], x[end], param, atol=tol, rtol=tol)
    return nothing
end

function 𝓇streams!(M, #output/solution vector
                   x, #output/solution coordinates
                   dIdx::Q, #version of schwarzschild equation
                   I₀::R, #initial irradiance
                   𝒜::S,
                   idx::Int,
                   g::Real, #gravity [m/s^2]
                   fT::U, #temperature profile fT(P)
                   fμ::V, #mean molar mass μ(T,P)
                   nstream::Int,
                   tol::Real #integrator error tolerance
                   )::Nothing where {Q,R<:Real,S<:AbstractAbsorber,U,V}
    @assert length(M) == length(x)
    #setup gaussian quadrature nodes
    𝓂, 𝒲 = streamnodes(nstream)
    #copy profiles to avoid interpolator race conditions
    fT, fμ = copyprofiles(fT, fμ)
    #solve schwarzschild w multiple streams, integrating over hemisphere
    for i ∈ 1:nstream-1
        𝓇stream!(M, x, dIdx, I₀, 𝒜, idx, g, 𝓂[i], fT, fμ, tol)
        M .*= 𝒲[i]/𝒲[i+1] #fancy way to do accumulating dot product
    end
    𝓇stream!(M, x, dIdx, I₀, 𝒜, idx, g, 𝓂[end], fT, fμ, tol)
    M .*= 𝒲[end]
    return nothing
end

#-------------------------------------------------------------------------------
# core function for whole atmosphere upward and downward monochromatic fluxes

function 𝓇monoflux!(M⁺, #downward monochromatic fluxes [W/m^2/cm^-1]
                    M⁻, #upward monochromatic fluxes [W/m^2/cm^-1]
                    P, #pressure coordinates of output
                    ω, #transformed pressure coords
                    ι, #transformed pressure coords
                    𝒜::Q,
                    idx::Int,
                    g::Real, #gravity [m/s^2]
                    fT::R, #temperature profile fT(P)
                    fμ::S, #mean molar mass μ(T,P)
                    fS::U, #incoming stellar radiation fS(ν) [W/m^2]
                    fa::V, #surface albedo fa(ν)
                    θₛ::Real, #stellar radiation angle, corresponds to cos(θ) = 2/3
                    nstream::Int, #number of streams to integrate in both directions
                    tol::Real) where {Q<:AbstractAbsorber,R,S,U,V}
    #setup
    @assert length(M⁻) == length(M⁺) == length(P) == length(ω) == length(ι)
    #surface pressure, assuming ascending pressures
    Pₛ = P[end]
    #surface temperature
    Tₛ = fT(Pₛ)
    #wavenumber
    ν = 𝒜.ν[idx]
    #angle factor for incoming stellar radiation
    𝓂 = 1/cos(θₛ)
    #downward stellar irradiance at ν
    Iₜ⁻ = fS(ν)
    #copy profiles to avoid possible interpolator race conditions
    fT, fμ, fS, fa = copyprofiles(fT, fμ, fS, fa)
    #wipe any previous values
    M⁺[:] .= zero(eltype(M⁺))
    M⁻[:] .= zero(eltype(M⁻))

    #===================================
    downgoing flux throughout atmosphere
    ===================================#
    #cosine of the stellar zenith angle
    c = cos(θₛ)
    #atmospheric contribution to downward flux
    𝓇streams!(M⁻, ι, dIdι, zero(Iₜ⁻), 𝒜, idx, g, fT, fμ, nstream, tol)
    #divide by c before adding the stellar irradiance in-place
    M⁻ ./= c
    #add downward stellar irradiance, **assuming no emission**
    𝓇stream!(M⁻, ι, dJdι, Iₜ⁻, 𝒜, idx, g, 𝓂, fT, fμ, tol)
    #multiply everything by c to get the true flux
    M⁻ .*= c

    #=================================
    upgoing flux throughout atmosphere
    =================================#
    #some of the downward stellar flux is reflected
    Iₛ⁺ = M⁻[end]*fa(ν)/π #Lambertian
    #and the surface emits some radiation
    Iₛ⁺ += planck(ν, Tₛ)
    #upward radiation streams
    𝓇streams!(M⁺, ω, dIdω, Iₛ⁺, 𝒜, idx, g, fT, fμ, nstream, tol)
    #reverse the upward flux to match the coordinate ordering of P and ι
    reverse!(M⁺)

    return nothing
end