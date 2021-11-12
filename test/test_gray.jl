using Cubature
using BasicInterpolators: CubicSplineInterpolator, NoBoundaries
using Base.Threads: @threads
using Test

using ClearSky
using ClearSky: 𝓇stream, dIdω, P2ω

const R = 8.31446262

##

totalopticaldepth(σ, g, μ, Pₛ, Pₜ=0) = dτdP(σ, g, μ)*(Pₛ - Pₜ)

function analyticaloutgoing(σ, g, μ, cₚ, Pₛ, Tₛ)
    #whole atmosphere optical depth, to a pressure of 1 Pa
    τ∞ = totalopticaldepth(σ, g, μ, Pₛ)
    #adiabat exponent
    γ = R/(μ*cₚ)
    #blackbody power
    P = stefanboltzmann(Tₛ)
    #analytical OLR solution (Principles of Planetary Climate equation 4.32)
    P*(exp(-τ∞) + τ∞^(-4γ)*hquadrature(τ′->exp(-τ′)*τ′^(4γ), 0, τ∞)[1])
end

function numericaloutgoing(σ, g, μ, cₚ, Pₛ, Tₛ, Pₜ=1e-3, tol=1e-9)
    #nice big, wide wavenumber sample
    ν = [logrange(1e-6, 1e5, 10000, 4); 1e6]
    #absorber wrapper
    U = UnifiedAbsorber(GrayGas(σ, ν))
    #dry adiabat, whole atmosphere
    Γ = DryAdiabat(Tₛ, Pₛ, cₚ, μ)
    # 1/cos(θ)
    m = 1.0
    #mean molar mass
    fμ(T,P) = μ
    #transformed coords
    ω₁, ω₂ = P2ω(Pₛ), P2ω(Pₜ)
    #compute OLR
    olr = zeros(length(ν))
    @threads for i ∈ eachindex(ν)
        I₀ = planck(ν[i], Tₛ)
        olr[i] = π*stream(dIdω, I₀, ω₁, ω₂, U, i, g, m, Γ, fμ, tol)
    end
    #make an interpolator for smooth integration
    ϕ = CubicSplineInterpolator([0.0; ν], [0.0; olr], NoBoundaries())
    #integrate it
    hquadrature(ν->ϕ(ν), 0.0, maximum(ν))[1]
end

##

#define atmospheric parameters
g = 10 #gravity [m/s^2]
μ = 0.01 #molar mass [kg/molecule]
cₚ = 1e3 #heat capacity [J/kg/K]
Pₛ = 1e5 #surface pressure [Pa]
Tₛ = 300 #surface temperature [K]
σ = 10 .^ range(-29, -23, length=10) #gray gas absorption coefs [cm^2/molecule]

#test accuracy with different optical depths
abserr = similar(σ)
relerr = similar(σ)
τ = similar(σ)
for (i,σᵢ) ∈ enumerate(σ)
    τ[i] = totalopticaldepth(σᵢ, g, μ, Pₛ)
    OLRₐ = analyticaloutgoing(σᵢ, g, μ, cₚ, Pₛ, Tₛ)
    OLRₙ = numericaloutgoing(σᵢ, g, μ, cₚ, Pₛ, Tₛ, 1e-6)
    err = OLRₙ - OLRₐ
    abserr[i] = abs(err)
    relerr[i] = abserr[i] / OLRₐ
    @test relerr[i] < 0.01
end

##

"""using Plots

p = plot(τ, relerr,
    xaxis=:log,
    yaxis=:log,
    xlabel="Whole Atmosphere Optical Depth",
    ylabel="OLR [W/m^2] Relative Error",
    legend=false)"""
    
