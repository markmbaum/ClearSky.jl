using Cubature
using BasicInterpolators: CubicSplineInterpolator, NoBoundaries
using Base.Threads: @threads

using ClearSky
using ClearSky: stream, dIdω, P2ω

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
    #gray gas object with uniform cross-section at all ν
    gray = GrayGas(σ, ν)
    #absorber wrapper
    G = GroupedAbsorber(gray)
    #dry adiabat, whole atmosphere
    Γ = DryAdiabat(Tₛ, Pₛ, cₚ, μ)
    #compute OLR
    olr = zeros(length(ν))
    @threads for i ∈ eachindex(ν)
        olr[i] = π*stream(dIdω, planck(ν[i], Tₛ), P2ω(Pₛ), P2ω(Pₜ), G, i, g, 1, Γ, (T,P)->μ, tol)
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
σ = 10 .^ range(-28, -23, length=50) #gray gas absorption coefs [cm^2/molecule]

τ = totalopticaldepth.(σ, g, μ, Pₛ)
OLRₐ = analyticaloutgoing.(σ, g, μ, cₚ, Pₛ, Tₛ)
OLRₙ = numericaloutgoing.(σ, g, μ, cₚ, Pₛ, Tₛ, 1e-6)
err = OLRₙ .- OLRₐ
abserr = abs.(err)
relerr = abserr ./ OLRₐ
@test all(relerr .< 0.1)

##

"""using Plots

p = plot(τ, relerr,
    xaxis=:log,
    yaxis=:log,
    xlabel="Whole Atmosphere Optical Depth",
    ylabel="OLR [W/m^2] Relative Error",
    legend=false)"""
    
