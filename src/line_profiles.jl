#precompute a few numbers
const sqπ = √π
const osqπln2 = 1/sqrt(π/log(2.0))
const sqln2 = sqrt(log(2.0))
const c2 = 100.0*𝐡*𝐜/𝐤

#-------------------------------------------------------------------------------
# wavenumber truncation of line shapes

cutline(ν, νl, Δνcut)::Bool = abs(ν - νl) > Δνcut ? true : false

function includedlines(ν::Real,
                       νl::AbstractVector{<:Real},
                       Δνcut::Real)::Vector{Int64}
    findall(x -> !cutline(ν, x, Δνcut), νl)
end

function includedlines(ν::AbstractVector{<:Real},
                       νl::AbstractVector{<:Real},
                       Δνcut::Real)::Vector{Int64}
    findall((νl .> minimum(ν) - Δνcut) .& (νl .< maximum(ν) + Δνcut))
end

#-------------------------------------------------------------------------------
# Chebyshev polynomial fit for Qref/Q

function chebyQrefQ(T::Real, n::Int64, a::Vector{Float64})::Float64
    #check the temperature range
    @assert TMIN <= T <= TMAX "temperature outside of Qref/Q interpolation range"
    #map T to [-1,1]
    τ = 2*(T - TMIN)/(TMAX - TMIN) - 1
    #values of first two chebys at τ
    c₁ = 1.0
    c₂ = τ
    #value of expansion after first two terms
    y = a[1] + a[2]*c₂
    for k = 3:n
        #next cheby value
        c₃ = 2*τ*c₂ - c₁
        #contribute to expansion
        y += a[k]*c₃
        #swap values
        c₁ = c₂
        c₂ = c₃
    end
    #return the inverse, Qref/Q
    return 1.0/y
end

#-------------------------------------------------------------------------------
# special strategy for line profiles from sorted vectors of wavenumbers

function surf!(σ::AbstractVector,
               f::F,
               ν::AbstractVector,
               νl::AbstractVector,
               Δνcut::Real,
               A::Vararg{Vector{Float64},N}) where {F<:Function, N}
    @assert all(diff(ν) .> 0) "wavenumber vectors must be sorted in ascending order"
    L = length(νl)
    jstart = 1 #tracking index to avoid searching from beginning every time
    for i = 1:length(ν)
        #find the first line that isn't cut off
        j = jstart
        while (j <= L) && cutline(ν[i], νl[j], Δνcut)
            j += 1
        end
        #only proceed if there is a line to include
        if j <= L
            #update the starting index for the search
            jstart = j
            #evaluate line profiles until one gets cut off, then move on
            while (j <= L) && !cutline(ν[i], νl[j], Δνcut)
                #let block is required for good performance
                args = let k = j
                    ntuple(n->A[n][k], N)
                end
                σ[i] += f(ν[i], νl[j], args...)
                j += 1
            end
        end
    end
end

#-------------------------------------------------------------------------------
# temperature scaling of line intensity

export scaleintensity

"""
    scaleintensity(S, νl, Epp, M, I, T)

Compute the [temperature scaling for line intensity](https://hitran.org/docs/definitions-and-units/#mjx-eqn-eqn-intensity-temperature-dependence).

# Arguments
* `S`: spectal line intensity at 296 K [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)]
* `νl`: wavenumber of line [cm``^{-1}``]
* `Epp`: lower-state energy of transition [cm``^{-1}``]
* `M`: [HITRAN molecular identification number](https://hitran.org/docs/molec-meta)
* `I`: [HITRAN local isotopologue number](https://hitran.org/docs/iso-meta/)
* `T`: temperature [K]
"""
function scaleintensity(S, νl, Epp, M::Int16, I::Int16, T)::Float64
    #arguments to exp
    a = -c2*Epp
    b = -c2*νl
    #numerator and denominator
    n = exp(a/T)*(1 - exp(b/T))
    d = exp(a/Tref)*(1 - exp(b/Tref))
    #check if there is an approximating function
    if MOLPARAM[M][10][I]
        QrefQ = chebyQrefQ(T, MOLPARAM[M][11][I], MOLPARAM[M][13][I])
    else
        throw("no interpolating polynomial available to compute Qref/Q for isotopologue $I of $(MOLPARAM[M][3]) ($(MOLPARAM[M][2]))")
        #QrefQ = (Tref/T)^1.5
    end
    #shifted line intensity
    S*QrefQ*(n/d)
end

function scaleintensity(sl::SpectralLines, i::Vector{Int64}, T)::Vector{Float64}
    S = view(sl.S, i)
    ν = view(sl.ν, i)
    Epp = view(sl.Epp, i)
    I = view(sl.I, i)
    scaleintensity.(S, ν, Epp, sl.M, I, T)
end

#-------------------------------------------------------------------------------
# doppler broadening

export αdoppler, fdoppler, doppler, doppler!

"""
    αdoppler(νl, μ, T)

Compute doppler (gaussian) broadening coefficient from line wavenumber `νl` [cm``^{-1}``], gas molar mass `μ` [kg/mole], and temperature `T` [K].
"""
αdoppler(νl, μ, T)::Float64 = (νl/𝐜)*sqrt(2.0*𝐑*T/μ)

function αdoppler(sl::SpectralLines, i::Vector{Int64}, T)::Vector{Float64}
    αdoppler.(view(sl.ν,i), view(sl.μ,i), T)
end

"""
    fdoppler(ν, νl, α)

Evaluate doppler (gaussian) profile

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `α`: doppler (gaussian) broadening coefficient
"""
fdoppler(ν, νl, α)::Float64 = exp(-(ν - νl)^2/α^2)/(α*sqπ)

"""
    doppler(ν, νl, S, α)

Evaluate doppler (gaussian) absoption cross-section [cm``^2``/molecule]

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `S`: line absoption intensity [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)]
* `α`: doppler (gaussian) broadening coefficient
"""
doppler(ν, νl, S, α)::Float64 = S*fdoppler(ν, νl, α)

"""
    doppler(ν, sl, T, P, Pₚ, Δνcut=25)

Evaluate a single doppler (gaussian) absoption cross-section [cm``^2``/molecule]. Temperature scaling and doppler profiles are evaluated along the way.

# Arguments
* `ν`: wavenumber indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function doppler(ν, sl::SpectralLines, T, P, Pₚ, Δνcut=25.0)::Float64
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    sum(doppler.(ν, view(sl.ν,i), S, α))
end

"""
    doppler!(σ, ν, sl, T, P, Pₚ, Δνcut=25)

Identical to [`doppler`](@ref), but fills the vector of cross-sections (`σ`) in-place.
"""
function doppler!(σ::AbstractVector,
                  ν::AbstractVector,
                  sl::SpectralLines,
                  T,
                  P,
                  Pₚ,
                  Δνcut=25.0)
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    surf!(σ, doppler, ν, view(sl.ν,i), Δνcut, S, α)
end

"""
    doppler(ν, sl, T, P, Pₚ, Δνcut=25)

Compute a vector of doppler (gaussian) absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and doppler profiles are evaluated along the way.

# Arguments
* `ν`: vector of wavenumbers indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function doppler(ν::AbstractVector,
                 sl::SpectralLines,
                 T,
                 P,
                 Pₚ,
                 Δνcut=25.0)::Vector{Float64}
    σ = zeros(Float64, length(ν))
    doppler!(σ, ν, sl, T, P, Pₚ, Δνcut)
    return σ
end

#-------------------------------------------------------------------------------
# pressure broadening

export γlorentz, florentz, lorentz, lorentz!

"""
    γlorentz(γa, γs, na, T, P, Pₚ)

Compute lorentzian broadening coefficient

# Arguments
* `γa`: air-broadened half width at half maximum (HWHM) [cm``^{-1}``/atm] at 296 K and 1 atm
* `γs`: self-broadened half width at half maximum (HWHM) [cm``^{-1}``/atm] at 296 K and 1 atm
* `na`: coefficient of temperature dependence of air-broadened half width
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
"""
function γlorentz(γa, γs, na, T, P, Pₚ)::Float64
    ((Tref/T)^na)*(γa*(P - Pₚ) + γs*Pₚ)/𝐀
end

function γlorentz(sl::SpectralLines, i::Vector{Int64}, T, P, Pₚ)::Vector{Float64}
    γlorentz.(view(sl.γa,i), view(sl.γs,i), view(sl.na,i), T, P, Pₚ)
end

"""
    florentz(ν, νl, γ)

Evaluate lorentz profile

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `γ`: lorentzian broadening coefficient
"""
florentz(ν, νl, γ)::Float64 = γ/(π*((ν - νl)*(ν - νl) + γ*γ))

"""
    lorentz(ν, νl, S, γ)

Evaluate lorentzian absoption cross-section [cm``^2``/molecule]

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `S`: line absoption intensity [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)]
* `γ`: lorentzian broadening coefficient
"""
lorentz(ν, νl, S, γ)::Float64 = S*florentz(ν, νl, γ)

"""
    lorentz(ν, sl, T, P, Pₚ, Δνcut=25)

Compute a single lorentzian absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and lorentzian profiles are evaluated along the way.

# Arguments
* `ν`: single wavenumber indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function lorentz(ν, sl::SpectralLines, T, P, Pₚ, Δνcut=25.0)::Float64
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    sum(lorentz.(ν, view(sl.ν,i), S, γ))
end

"""
    lorentz!(σ, ν, sl, T, P, Pₚ, Δνcut=25)

Identical to [`lorentz`](@ref), fills the vector of cross-sections (`σ`) in-place.
"""
function lorentz!(σ::AbstractVector,
                  ν::AbstractVector,
                  sl::SpectralLines,
                  T,
                  P,
                  Pₚ,
                  Δνcut=25.0)
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    surf!(σ, lorentz, ν, view(sl.ν,i), Δνcut, S, γ)
end

"""
    lorentz(ν, sl, T, P, Pₚ, Δνcut=25)

Compute a vector of lorentzian absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and lorentzian profiles are evaluated along the way.

# Arguments
* `ν`: vector of wavenumbers indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function lorentz(ν::AbstractVector,
                 sl::SpectralLines,
                 T,
                 P,
                 Pₚ,
                 Δνcut=25.0)::Vector{Float64}
    σ = zeros(Float64, length(ν))
    lorentz!(σ, ν, sl, T, P, Pₚ, Δνcut)
    return σ
end

#-------------------------------------------------------------------------------
# voigt profile

export fvoigt, voigt, voigt!

"""
    fvoigt(ν, νl, α, γ)

Evaluate Voigt profile

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `α`: doppler (gaussian) broadening coefficient
* `γ`: lorentzian broadening coefficient
"""
function fvoigt(ν, νl, α, γ)::Float64
    #inverse of the doppler parameter
    β = 1/α
    #factor for real and complex parts of Faddeeva args, avoiding β division
    d = sqln2*β
    #arguments to Faddeeva function
    x = (ν - νl)*d
    y = γ*d
    #evaluate real part of Faddeeva function
    f = faddeyeva(x,y)
    #final calculation, avoiding α division by using β again
    osqπln2*β*f
end

"""
    voigt(ν, νl, S, α, γ)

Evaluate Voigt absoption cross-section [cm``^2``/molecule]

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `S`: line absoption intensity [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)]
* `α`: doppler (gaussian) broadening coefficient
* `γ`: lorentzian broadening coefficient
"""
voigt(ν, νl, S, α, γ)::Float64 = S*fvoigt(ν, νl, α, γ)

"""
    voigt(ν, sl::SpectralLines, T, P, Pₚ, Δνcut=25)

Evaluate Voigt absorption cross-section at a single wavenumber.
"""
function voigt(ν, sl::SpectralLines, T, P, Pₚ, Δνcut=25.0)::Float64
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    sum(voigt.(ν, view(sl.ν,i), S, α, γ))
end

"""
    voigt!(σ, ν, sl, T, P, Pₚ, Δνcut=25)

Identical to [`voigt`](@ref), but fills the vector of cross-sections (`σ`) in-place.
"""
function voigt!(σ::AbstractVector,
                ν::AbstractVector,
                sl::SpectralLines,
                T,
                P,
                Pₚ,
                Δνcut=25.0)
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    surf!(σ, voigt, ν, view(sl.ν,i), Δνcut, S, α, γ)
end

"""
    voigt(ν, sl, T, P, Pₚ, Δνcut=25)

Compute a vector of Voigt absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and Voigt profiles are evaluated along the way.

# Arguments
* `ν`: vector of wavenumbers indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function voigt(ν::AbstractVector,
               sl::SpectralLines,
               T,
               P,
               Pₚ,
               Δνcut=25.0)::Vector{Float64}
    σ = zeros(Float64, length(ν))
    voigt!(σ, ν, sl, T, P, Pₚ, Δνcut)
    return σ
end

#-------------------------------------------------------------------------------
# sublorentzian profile for CO2
# Perrin, M. Y., and J. M. Hartman. “Temperature-Dependent Measurements and Modeling of Absorption by CO2-N2 Mixtures in the Far Line-Wings of the 4.3 Μm CO2 Band.” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 42, no. 4, 1989, pp. 311–17.

export ΧPHCO2, PHCO2, PHCO2!

"""
    ΧPHCO2(ν, νl, T)

Compute the `Χ` (Chi) factor for sub-lorentzian CO2 line profiles, as in
* [Perrin, M. Y., and J. M. Hartman. “Temperature-Dependent Measurements and Modeling of Absorption by CO2-N2 Mixtures in the Far Line-Wings of the 4.3 Μm CO2 Band.” Journal of Quantitative Spectroscopy and Radiative Transfer, vol. 42, no. 4, 1989, pp. 311–17.](https://www.sciencedirect.com/science/article/abs/pii/0022407389900770)

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `T`: temperature [K]
"""
function ΧPHCO2(ν, νl, T)::Float64
    Δν = abs(ν - νl)
    if Δν < 3.0
        return 1.0
    end
    B1 = 0.0888 - 0.16*exp(-0.0041*T)
    if Δν < 30.0
        return exp(-B1*(Δν - 3.0))
    end
    B2 = 0.0526*exp(-0.00152*T)
    if Δν < 120.0
        return exp(-B1*27.0 - B2*(Δν - 30.0))
    end
    return exp(-B1*27.0 - B2*90.0 - 0.0232*(Δν - 120.0))
end

"""
    PHCO2(ν, νl, S, α)

Evaluate Perrin & Hartman sub-lorentzian absoption cross-section [cm``^2``/molecule] for CO2

# Arguments
* `ν`: profile evaluation wavenumber [cm``^{-1}``]
* `νl`: wavenumber of absorption line [cm``^{-1}``]
* `T`: temperature [K]
* `S`: line absoption intensity [cm``^{-1}``/(molecule``\\cdot``cm``^{-2}``)]
* `α`: doppler (gaussian) broadening coefficient
* `γ`: lorentzian broadening coefficient
"""
function PHCO2(ν, νl, T, S, α, γ)::Float64
    Χ = ΧPHCO2(ν, νl, T)
    voigt(ν, νl, S, α, Χ*γ)
end

"""
    PHCO2(ν, sl, T, P, Pₚ, Δνcut=500)

Compute a single Perrin & Hartman sub-lorentzian CO2 absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and profiles are evaluated along the way.

# Arguments
* `ν`: single wavenumber indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function PHCO2(ν, sl::SpectralLines, T, P, Pₚ, Δνcut=500.0)::Float64
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    sum(PHCO2.(ν, view(sl.ν,i), T, S, α, γ))
end

"""
    PHCO2!(σ, ν, sl, T, P, Pₚ, Δνcut=500)

Identical to [`PHCO2`](@ref), but fills the vector of cross-sections (`σ`) in-place.
"""
function PHCO2!(σ::AbstractVector,
                ν::AbstractVector,
                sl::SpectralLines,
                T,
                P,
                Pₚ,
                Δνcut=500.0)
    i = includedlines(ν, sl.ν, Δνcut)
    S = scaleintensity(sl, i, T)
    α = αdoppler(sl, i, T)
    γ = γlorentz(sl, i, T, P, Pₚ)
    f(ν, νl, S, α, γ) = PHCO2(ν, νl, T, S, α, γ) #shove T into the function
    surf!(σ, f, ν, view(sl.ν,i), Δνcut, S, α, γ)
end

"""
    PHCO2(ν, sl, T, P, Pₚ, Δνcut=500)

Compute a vector of Perrin & Hartman sub-lorentzian CO2 absorption cross-sections [cm``^2``/molecule] from a [`SpectralLines`](@ref) object. Temperature scaling and profiles are evaluated along the way.

# Arguments
* `ν`: vector of wavenumbers indicating where to evaluate [cm``^{-1}``]
* `sl`: [`SpectralLines`](@ref)
* `T`: temperature [K]
* `P`: air pressure [Pa]
* `Pₚ`: partial pressure [Pa]
* `Δνcut`: profile truncation distance [cm``^{-1}``]
"""
function PHCO2(ν::AbstractVector,
               sl::SpectralLines,
               T,
               P,
               Pₚ,
               Δνcut=500.0)::Vector{Float64}
    σ = zeros(Float64, length(ν))
    PHCO2!(σ, ν, sl, T, P, Pₚ, Δνcut)
    return σ
end
