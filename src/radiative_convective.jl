export RCM

#------------------------------------------------------------------------------
# definition & construction

struct RCM{A<:Real,𝒻₁,𝒻₂,𝒻₃,𝒻₄,C<:AbstractNumericalCore}
    Pₑ::Vector{A}
    Tₑ::Vector{A} #fossil initial condition
    Pᵣ::Vector{A}
    P::Vector{A}
    T::Vector{A}
    g::A #gravity [m/s^2]
    𝒻μ::𝒻₁ #mean molar mass function (T,P)->μ [kg/mole]
    𝒻S::𝒻₂ #stellar insolation at TOA ν->I [W/m^2/cm^-1]
    𝒻a::𝒻₃ #surface albedo ν->a [unitless]
    𝒻cₚ::𝒻₄ #heat capacity function (T,P)->cₚ [J/kg/K]
    cₛ::A #surface heat capacity [J/m^2/K] 
    𝒜::AcceleratedAbsorber
    ν::Vector{Float64}
    core::C
    F::FluxPack{A}
    np::Int64 #number of atmospheric cells (including surface)
    nν::Int64
    nᵤ::Vector{Int64} #update counters for rows, just curious
    R::Vector{A} #net radiative fluxes at cell edges
    H::Vector{A} #cell heating rate [K/s]
    J::Matrix{A} #Jacobian
end

function Base.show(io::IO, ℛ::RCM{A}) where {A}
    @unpack Pₑ, ν, T, core, g, np, nν = ℛ
    print(io, "RCM{$A,$(typeof(core))}\n")
    print(io, "  gravity is $g m/s^2\n")
    Pmin, Pmax = round(minimum(Pₑ), sigdigits=4), round(maximum(Pₑ), sigdigits=4)
    print(io, "  $np pressure layers ∈ [$Pmin, $Pmax] Pa\n")
    νmin, νmax = round(minimum(ν), sigdigits=4), round(maximum(ν), sigdigits=4)
    print(io, "  $nν wavenumber coordinates ∈ [$νmin, $νmax] cm^-1\n")
    Tmin, Tmax = round(minimum(T), sigdigits=4), round(maximum(T), sigdigits=4)
    print(io, "  temperature currently ∈ [$Tmin, $Tmax] K")
end

function RCM(Pₑ::AbstractVector{<:Real},
             Tₑ::AbstractVector{<:Real},
             g::Real,
             𝒻μ,
             𝒻S,
             𝒻a,
             𝒻cₚ,
             cₛ::Real,
             absorbers...;
             core::AbstractNumericalCore=Discretized(),
             radmul::Int=2)
    #demand pressure coordiantes are ascending
    W = typeof(float(Tₑ[1]))
    idx = sortperm(Pₑ)
    Pₑ = collect(W, Pₑ[idx])
    Tₑ = collect(W, Tₑ[idx])
    #check the numbers
    np = length(Pₑ)
    @assert length(Tₑ) == np "must have same number of initial temperature and pressure values"
    #values at midpoints between edges (cell centers) and surface
    P = similar(Pₑ)
    T = similar(Tₑ)
    for i ∈ 1:np-1
        P[i] = (Pₑ[i] + Pₑ[i+1])/2
        T[i] = (Tₑ[i] + Tₑ[i+1])/2
    end
    P[end] = Pₑ[end]
    T[end] = Tₑ[end]
    #build up extra radiative nodes with weighted averaging
    @assert ((radmul % 2 == 0) | (radmul == 1)) "radmul must be an even integer or 1"
    nrad = radmul*(np - 1) + 1
    Pᵣ = similar(Pₑ, nrad)
    P₁ = @view Pₑ[1:end-1]
    Pᵣ[1:np-1] .= P₁
    P₂ = @view Pₑ[2:end]
    i = np
    for j ∈ 2:radmul
        w₁ = j - 1
        w₂ = radmul - w₁
        @. Pᵣ[i:i+np-2] = (w₁*P₁ + w₂*P₂)/radmul
        i += np - 1
    end
    Pᵣ[end] = Pₑ[end]
    sort!(Pᵣ)
    #parse absorbers
    𝒰, _, nν = unifyabsorbers(absorbers)
    #prepare accelerated absorption calculations
    𝒜 = AcceleratedAbsorber(Tₑ, Pₑ, 𝒰)
    #handle types
    cₛ = convert(W, cₛ)
    g = convert(W, g)
    #make room for radiative calculations
    F = FluxPack(nrad, nν, W)
    #cross-section update tracker
    nᵤ = zeros(Int64, np)
    #make room for heating calculations
    R = zeros(W, np)
    H = zeros(W, np)
    J = zeros(W, np, np)
    #construct
    RCM(Pₑ, Tₑ, Pᵣ, P, T, g, 𝒻μ, 𝒻S, 𝒻a, 𝒻cₚ, cₛ, 𝒜, 𝒜.ν, core, F, np, nν, nᵤ, R, H, J)
end

#------------------------------------------------------------------------------
# time derivatives and such

export heating!
function heating!(ℛ::RCM, δᵤ::Real=1)::Nothing
    @unpack Pₑ, Tₑ, Pᵣ, P, T, g, 𝒻μ, 𝒻S, 𝒻a, 𝒻cₚ, cₛ, 𝒜, core, F, np, nν, nᵤ, R, H = ℛ
    
    𝒻T = AtmosphericProfile(P, T)
    radiate!(F, core, Pᵣ, g, 𝒻T, 𝒻μ, 𝒻S, 𝒻a, 𝒜)

    #====================================================================
    Coordinates are a little tricky here. The "upward" flux Fup is going
    out of the atmosphere, away from the surface. However, the positive
    direction is the other way because pressure is increasing toward the
    surface and the model discretization is along pressure coordinates.
    So, although the net flux is typically F⁺ - F⁻, it must be reversed
    here to respect the true positive direction.
    ====================================================================#
    𝒻Fnet = AtmosphericProfile(Pᵣ, F.Fnet)
    @. R = -𝒻Fnet(Pₑ) #F.F⁻ - F.F⁺

    #==========================
    heating rates
    ==========================#
    for i ∈ 1:np-1
        #heat capacity [J/kg/K]
        cₚ = 𝒻cₚ(T[i],P[i])
        #pressure layer thickness [Pa]
        ΔP = Pₑ[i+1] - Pₑ[i]
        #radiation out of the cell minus radiation into it
        ΔR = R[i] - R[i+1]
        #cell heating rate [K/s]
        H[i] = (g/cₚ)*ΔR/ΔP
    end
    #surface heating
    H[end] = R[end]/cₛ
    #H[end] = (F.F⁺[1] - F.F⁻[end])/cₛ

    return nothing
end

export step!
function step!(ℛ::RCM, Δt::Real)::Nothing
    heating!(ℛ)
    ℛ.T .+= Δt*ℛ.H
    return nothing
end

export jacobian!
function jacobian!(ℛ::RCM, ϵ::Real=1)::Nothing
    heating!(ℛ)
    H = ℛ.H[:]
    T = ℛ.T
    for i ∈ 1:ℛ.np
        #perturb
        T[i] += ϵ
        #evaluate
        heating!(ℛ)
        #finite diff
        for j ∈ 1:length(H)
            ℛ.J[j,i] = (ℛ.H[j] - H[j])/ϵ
        end
        #restore
        T[i] -= ϵ
    end
    return nothing
end