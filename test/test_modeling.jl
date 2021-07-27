Pₛ = 1e5
Pₜ = 1
g = 9.8
fμ(T,P) = 0.029
θ = acos(2/3)
fstellar(ν) = 1336*normplanck(ν, 5500)
falbedo(ν) = 0.2

ν = [1:2500; 2510:10:75000]
Ω = AtmosphericDomain((100,350), 8, (0.9,2e5), 16)

co2 = WellMixedGas("HITRAN/CO2.par", 400e-6, ν, Ω, PHCO2!)
ch4 = WellMixedGas("HITRAN/CH4.par", 1e-5, ν, Ω, lorentz!)
h2o = VariableGas("HITRAN/H2O.par", (T,P)->0.01, ν, Ω, voigt!)

co2co2 = CIATables("HITRAN/CO2-CO2_2018.cia", verbose=false)
co2ch4 = CIATables("HITRAN/CO2-CH4_2018.cia", verbose=false)

Γ = MoistAdiabat(280, Pₛ, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Tstrat=190)
h2o = condensibleprofile(h2o, Γ, T->0.8*psatH2O(T))
P = reverse(logrange(1e5, 1, 250))
T = Γ.(P)
A = AcceleratedAbsorber(P, T, co2, ch4, h2o, co2co2, co2ch4, (ν,T,P)->1e-30)

issane(x) = all(@. !isnan(x)) & all(@. !isinf(x))

τ = opticaldepth(Pₛ, Pₜ, g, Γ, fμ, θ, A)
@test issane(τ)

t = transmittance(Pₛ, Pₜ, g, Γ, fμ, θ, A)
@test issane(t)

olr = outgoing(Pₛ, g, Γ, fμ, A)
@test issane(olr)

Fₜ⁻, Fₜ⁺ = topfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
@test issane(Fₜ⁻)
@test issane(Fₜ⁺)

Fₛ⁻, Fₛ⁺ = bottomfluxes(Pₛ, g, Γ, fμ, fstellar, falbedo, A)
@test issane(Fₛ⁻)
@test issane(Fₛ⁺)