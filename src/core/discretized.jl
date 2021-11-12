#shifts and memoizes the Gauss-Lobatto nodes for optical depth in a layer
@memoize function lobattonodes(n::Int64)::NTuple{2,Vector{Float64}}
    #gauss-lobatto nodes for integrating optical depth in pressure layer
    x, w = gausslobatto(n)
    #simply shift to [0,1]
    𝓍 = @. (x + 1)/2
    𝓌 = @. w/2
    return 𝓍, 𝓌
end

function lobattoevaluations(P, 𝒻T::Q, 𝒻μ::R, nlobatto::Int) where {Q,R}
    #setup new arrays
    np = length(P)
    U = typeof(𝒻T(P[1]))
    T = Array{U,2}(undef, nlobatto, np-1)
    μ = Array{U,2}(undef, nlobatto, np-1)
    #Gauss-Lobatto quadrature nodes, prescaled to [0,1]
    𝓍, _ = lobattonodes(nlobatto)
    for i ∈ 1:nlobatto
        @inbounds for j ∈ 1:np-1
            ΔP = P[j+1] - P[j]
            Pᵢ = P[j] + ΔP*𝓍[i]
            Tᵢ = 𝒻T(Pᵢ)
            μᵢ = 𝒻μ(Tᵢ,Pᵢ)
            T[i,j] = Tᵢ
            μ[i,j] = μᵢ
        end
    end
    return T, μ
end

function temperaturederivs(P, 𝒻T::Q, ϵ::Real=0.01) where {Q}
    n = length(P)
    ∂T = similar(P)
    δP = (P[2] - P[1])*ϵ
    ∂T[1] = (𝒻T(P[1] + δP) - 𝒻T(P[1]))/δP
    for i = 2:n-1
        δP = min(P[i+1] - P[i], P[i] - P[i-1])*ϵ
        ∂T[i] = (𝒻T(P[i] + δP) - 𝒻T(P[i] - δP))/(2*δP)
    end
    δP = (P[n] - P[n-1])*ϵ
    ∂T[n] = (𝒻T(P[n]) - 𝒻T(P[n] - δP))/δP
    return ∂T
end

function planckevaluations(P, ν, 𝒻T::Q) where {Q}
    U = typeof(planck(ν[1], 𝒻T(P[1])))
    B = Array{U,2}(undef, length(P), length(ν))
    #B′ = similar(B)
    @threads for i ∈ eachindex(P)
        Tᵢ = 𝒻T(P[i])
        @inbounds for j ∈ eachindex(ν)
            B[i,j] = planck(ν[j], Tᵢ)
            #B′[i,j] = dplanck(ν[j], Tᵢ)
        end
    end
    return B#, B′
end

function βevaluations(P, g, 𝒻T::Q, 𝒻μ::R, 𝒜::S) where {Q,R,S<:AbstractAbsorber}
    ν = 𝒜.ν
    C = 1e-4*𝐍𝐚/g
    β₀ = 𝜷(P[1], 𝒻T(P[1]), 𝒻μ(P[1], 𝒻T(P[1])), C, 1, 𝒜)
    β = Array{typeof(β₀),2}(undef, length(P), length(ν))
    @threads for j ∈ eachindex(ν)
        @inbounds for i ∈ eachindex(P)
            Pᵢ = P[i]
            Tᵢ = 𝒻T(Pᵢ)
            μᵢ = 𝒻μ(Tᵢ, Pᵢ)
            β[i,j] = 𝜷(Pᵢ, Tᵢ, μᵢ, C, j, 𝒜)
        end
    end
    return β
end

function 𝜷(P, T, μ, C, idx::Int, 𝒜::Q) where {Q<:AbstractAbsorber}
    #total cross-section [cm^2/molecule]
    σ = Σ(𝒜, idx, T, P)
    #absorption coefficient [m*s^2/kg], represents dτ/dP
    C*(σ/μ)
end

#the "linear in τ" approximation for planck emission of a discrete layer, see:
#  Clough, S. A., Iacono, M. J. & Moncet, J.-L. Line-by-line calculations of atmospheric fluxes and cooling rates: Application to water vapor. J. Geophys. Res. 97, 15761 (1992).
layerplanck(B₁, B₂, τ, t) = B₂*(1.0 - t) - (B₁ - B₂)*t + (1.0 - t)*(B₁ - B₂)/τ

layerplanck(B₁, B₂, τ) = layerplanck(B₁, B₂, τ, exp(-τ))

#-------------------------------------------------------------------------------
# optical depth paths

function 𝒹depth(P::AbstractVector{<:Real}, #ascending
                T::Matrix{<:Real}, #evaluated at Lobatto nodes
                μ::Matrix{<:Real}, #evaluated at Lobatto nodes
                𝒜::Q,
                idx::Int,
                C::Real,
                𝓂::Real, #1/cos(θ)
                nlobatto::Int) where {Q<:AbstractAbsorber}
    #setup
    np = length(P)
    𝓍, 𝓌 = lobattonodes(nlobatto)

    #first quadrature point is right at the surface
    @inbounds β₁ = 𝜷(P[1], T[1,1], μ[1,1], C, idx, 𝒜)
    #accumulating value
    τ = zero(β₁)
    #loop over pressure layers
    @inbounds for i ∈ 1:np-1
        #beginning and end of pressure layer
        P₁, P₂ = P[i], P[i+1]
        ΔP = P₂ - P₁
        #integrate optical depth over layer with Gauss-Lobotto quadrature
        τᵢ = zero(eltype(τ))
        #first β is reused from end of previous layer
        τᵢ += (ΔP*𝓌[1])*β₁
        for n ∈ 2:nlobatto-1
            #pressure coordinate of quadrature node, tricky because decreasing!
            Pₙ = P₁ + ΔP*𝓍[n]
            #absorption coefficient
            β = 𝜷(Pₙ, T[n,i], μ[n,i], C, idx, 𝒜)
            #contribute to optical depth integral
            τᵢ += (ΔP*𝓌[n])*β
        end
        #final point at end of pressure layer
        βₙ = 𝜷(P₂, T[end,i], μ[end,i], C, idx, 𝒜)
        τᵢ += (ΔP*𝓌[end])*βₙ
        #the end is now the beginning
        β₁ = βₙ
        #store τ
        τ += τᵢ*𝓂
    end
    return τ
end

function 𝒹depth!(τ::AbstractVector{<:Real}, #to be filled
                 P::AbstractVector{<:Real}, #ascending
                 T::Matrix{<:Real}, #evaluated at Lobatto nodes
                 μ::Matrix{<:Real}, #evaluated at Lobatto nodes
                 𝒜::Q,
                 idx::Int,
                 C::Real,
                 nlobatto::Int)::Nothing where {Q<:AbstractAbsorber}
    #setup
    np = length(P)
    𝓍, 𝓌 = lobattonodes(nlobatto)
    τmin = 1e-6*one(eltype(τ))

    #first quadrature point
    β₁ = 𝜷(P[1], T[1,1], μ[1,1], C, idx, 𝒜)
    #loop over pressure layers
    @inbounds for i ∈ 1:np-1
        #beginning and end of pressure layer
        P₁, P₂ = P[i], P[i+1]
        ΔP = P₂ - P₁
        #integrate optical depth over layer with Gauss-Lobotto quadrature
        τᵢ = zero(eltype(τ))
        #first β is reused from end of previous layer
        τᵢ += (ΔP*𝓌[1])*β₁
        for n ∈ 2:nlobatto-1
            #pressure coordinate of quadrature node, tricky because decreasing!
            Pₙ = P₁ + ΔP*𝓍[n]
            #intra-layer absorption coefficient isn't pre-calculated
            βₙ = 𝜷(Pₙ, T[n,i], μ[n,i], C, idx, 𝒜)
            #contribute to optical depth integral
            τᵢ += (ΔP*𝓌[n])*βₙ
        end
        #final point at end of pressure layer
        βₙ = 𝜷(P₂, T[end,i], μ[end,i], C, idx, 𝒜)
        τᵢ += (ΔP*𝓌[end])*βₙ
        #the end is now the beginning for the next layer
        β₁ = βₙ
        #put a floor on τ and store it
        τ[i] = max(τᵢ, τmin)
    end
    nothing
end

#-------------------------------------------------------------------------------
# streams/fluxes up the atmosphere, no storage

function 𝒹streams(P::AbstractVector{<:Real},
                  T::Matrix{<:Real}, #evaluated at Lobatto nodes
                  μ::Matrix{<:Real}, #evaluated at Lobatto nodes
                  𝒜::Q,
                  idx::Int,
                  I₀::Real,
                  C::Real,
                  nstream::Int,
                  nlobatto::Int) where {Q<:AbstractAbsorber}
    #setup
    np = length(P)
    𝓂, 𝒲 = streamnodes(nstream)
    𝓍, 𝓌 = lobattonodes(nlobatto)
    ν = 𝒜.ν[idx]
    τmin = 1e-6*one(eltype(τ))

    #output flux [W/m^2/cm^-1]    
    F = zero(I₀)
    #𝜷 coefficient at the surface
    βₛ = 𝜷(P[1], T[1,1], μ[1,1], C, idx, 𝒜)
    #iterate through streams of different angles
    
    @inbounds for k ∈ 1:nstream
        #initial planck emission from surface
        I = convert(typeof(F), I₀)
        #first quadrature point is right at the surface
        β₁ = βₛ
        #loop over pressure layers
        for i ∈ 1:np-1
            #beginning and end of pressure layer
            P₁, P₂ = P[i], P[i+1]
            ΔP = P₁ - P₂
            #integrate optical depth over layer with Gauss-Lobotto quadrature
            τ = zero(I)
            #first β is reused from end of previous layer
            τ += (ΔP*𝓌[1])*β₁*𝓂[k]
            for n ∈ 2:nlobatto-1
                #pressure coordinate of quadrature node, tricky because decreasing!
                Pₙ = P₁ - ΔP*𝓍[n]
                #absorption coefficient
                β = 𝜷(Pₙ, T[n,i], μ[n,i], C, idx, 𝒜)
                #contribute to optical depth integral
                τ += (ΔP*𝓌[n])*β*𝓂[k]
            end
            #final point at end of pressure layer
            βₙ = 𝜷(P₂, T[end,i], μ[end,i], C, idx, 𝒜)
            τ += (ΔP*𝓌[end])*βₙ*𝓂[k]
            τ = max(τ, τmin)
            #the end is now the beginning
            β₁ = βₙ
            #transmittance through layer
            t = exp(-τ)
            #CIM planck emission approximation
            B₁ = planck(ν, T[1,i])
            B₂ = planck(ν, T[end,i])
            Bₑ = layerplanck(B₁, B₂, τ, t)
            #update irradiance
            I = I*t + Bₑ
        end
        F += 𝒲[k]*I
    end
    return F
end

#-------------------------------------------------------------------------------
# core function for whole atmosphere upward and downward monochromatic fluxes

function 𝒹monoflux!(M⁺, #downward monochromatic fluxes [W/m^2/cm^-1]
                    M⁻, #upward monochromatic fluxes [W/m^2/cm^-1]
                    τ, #optical depth of each layer [-] length(P) - 1
                    P::AbstractVector{<:Real}, #pressure coordinates of output
                    B::AbstractVector{<:Real}, #pre-evaluated Planck at P coordinates
                    ν::Real,
                    fS::Q, #incoming stellar radiation fS(ν) [W/m^2]
                    fa::R, #surface albedo fa(ν)
                    θₛ::Real, #stellar radiation angle, corresponds to cos(θ) = 2/3
                    nstream::Int #number of streams to integrate in both directions
                    )::Nothing where {Q,R}
    #setup
    np = length(P)
    #number of layers
    L = np - 1 
    #check sizes
    @assert np == length(M⁻) == length(M⁺)
    @assert L == length(τ) 
    #gauss nodes for integrating streams
    𝓂, 𝒲 = streamnodes(nstream)
    #cosine of stellar zenith angle
    c = cos(θₛ)
    #avoid possible interpolator race conditions
    fS, fa = copyprofiles(fS, fa)
    #wipe any previous values
    M⁺[:] .= zero(eltype(M⁺))
    M⁻[:] .= zero(eltype(M⁻))

    #=============================
    #downward atmospheric emission
    =============================#
    @inbounds for k ∈ 1:nstream
        #initial irradiance
        I = zero(eltype(M⁻))
        for i ∈ 1:L
            τᵢ = τ[i]*𝓂[k] #weight by 1/cos(θ)
            tᵢ = exp(-τᵢ)
            Bₑ = layerplanck(B[i], B[i+1], τᵢ, tᵢ)
            I = I*tᵢ + Bₑ
            M⁻[i+1] += 𝒲[k]*I
        end
    end

    #===========================================================
    downward stellar flux, assuming absorption only, no emission
    ===========================================================#
    M⁻[1] += c*fS(ν) #initial flux scaled by angle
    Mₛ = M⁻[1]
    @inbounds for i ∈ 1:L
        τᵢ = τ[i]/c
        tᵢ = exp(-τᵢ)
        Mₛ *= tᵢ
        M⁻[i+1] += Mₛ
    end

    #=================
    outgoing radiation
    =================#
    Iₛ⁺ = M⁻[end]*fa(ν)/π + B[end] #Lambertian reflection and planck emission
    M⁺[end] = Iₛ⁺*π #upward surface flux
    @inbounds for k ∈ 1:nstream
        I = Iₛ⁺
        for i ∈ L:-1:1
            τᵢ = τ[i]*𝓂[k] #weight by 1/cos(θ)
            tᵢ = exp(-τᵢ)
            Bₑ = layerplanck(B[i+1], B[i], τᵢ, tᵢ)
            I = I*tᵢ + Bₑ
            M⁺[i] += 𝒲[k]*I
        end
    end

    #stop at nothing!
    nothing
end