#-------------------------------------------------------------------------------
#caching function for gaussian quadrature of multiple streams over the azimuth

@memoize function streamnodes(n::Int64)::NTuple{2,Vector{Float64}}
    #gauss-legendre quadrature points and weights in [-1,1]
    x, w = gausslegendre(n)
    #map to θ ∈ [0,2π] and precompute trig factors
    𝓂 = Vector{Float64}(undef,n)
    𝒲 = Vector{Float64}(undef,n)
    for i ∈ 1:n
        #azimuth angle of stream
        θᵢ = (π/2)*(x[i] + 1)/2
        #mapped weight
        wᵢ = (π/2)*w[i]/2
        #factor of relevance to optical depth and such
        𝓂[i] = 1/cos(θᵢ)
        #weight with precomputed trig factors for hemispheric geometry
        𝒲[i] = 2π*wᵢ*cos(θᵢ)*sin(θᵢ)
    end
    return 𝓂, 𝒲
end

#-------------------------------------------------------------------------------

copyprofile(x::MoistAdiabat) = copy(x)
copyprofile(x::AtmosphericProfile) = copy(x)
copyprofile(x::OneDimensionalInterpolator) = copy(x)
copyprofile(x) = x
copyprofiles(x...) = map(copyprofile, x)

#-------------------------------------------------------------------------------
#little types to dispatch on numerical cores/methods

export Radau, Discretized

abstract type AbstractNumericalCore end

#--------------------------------
 
struct Radau <: AbstractNumericalCore
    nstream::Int64
    tol::Float64
end

function Radau(; nstream::Int=5, tol::AbstractFloat=1e-5)
    Radau(nstream, tol)
end

function Base.show(io::IO, r::Radau)
    print(io, "Radau(nstream=$(r.nstream), tol=$(r.tol))")
end

#--------------------------------

struct Discretized <: AbstractNumericalCore
    nstream::Int64
    nlobatto::Int64
end

function Discretized(; nstream::Int=5, nlobatto::Int=2)
    Discretized(nstream, nlobatto)
end

function Base.show(io::IO, d::Discretized)
    print(io, "Discretized(nstream=$(d.nstream), nlobatto=$(d.nlobatto))")
end

#-------------------------------------------------------------------------------
#a simple container of arrays needed for whole atmosphere characterization

export FluxPack

struct FluxPack{T<:Real}
    τ::Matrix{T}
    M⁺::Matrix{T}
    M⁻::Matrix{T}
    F⁺::Vector{T}
    F⁻::Vector{T}
    Fnet::Vector{T}
end

function Base.show(io::IO, M::FluxPack{T}) where {T}
    np, nν = size(M.M⁺)
    print(io, "size ($np x $nν) FluxPack{$T}\n")
    print("  TOA:\n")
    print(io, "    outgoing radiation = $(M.F⁺[1]) W/m^2\n")
    print(io, "    incoming radiation = $(M.F⁻[1]) W/m^2\n")
    print("  surface:\n")
    print(io, "    outgoing radiation = $(M.F⁺[end]) W/m^2\n")
    print(io, "    incoming radiation = $(M.F⁻[end]) W/m^2\n")
end

function FluxPack(np::Int, nν::Int, T::Type=Float64)
    FluxPack(
        zeros(T, np - 1, nν),
        zeros(T, np, nν),
        zeros(T, np, nν),
        zeros(T, np),
        zeros(T, np),
        zeros(T, np)
    )
end

function FluxPack(P::AbstractVector, ν::AbstractVector, T::Type=Float64)
    FluxPack(length(P), length(ν), T)
end

Base.size(F::FluxPack) = size(F.M⁺)

for op ∈ (:+, :-, :\, :*)
    @eval begin
        Base.$op(A::FluxPack, B::FluxPack) = FluxPack(
            broadcast($op, A.τ, B.τ),
            broadcast($op, A.M⁺, B.M⁺),
            broadcast($op, A.M⁻, B.M⁻),
            broadcast($op, A.F⁺, B.F⁺),
            broadcast($op, A.F⁻, B.F⁻),
            broadcast($op, A.Fnet, B.Fnet)
        )
    end
end

#-------------------------------------------------------------------------------

function ∫F!(F⁺, F⁻, M⁺, M⁻, ν)::Nothing
    n, m = size(M⁻)
    @assert size(M⁺) == size(M⁻)
    @assert length(ν) == m 
    @assert length(F⁻) == length(F⁺) == n
    @inbounds for i ∈ 1:n
        Mᵢ⁺ = @view M⁺[i,:]
        Mᵢ⁻ = @view M⁻[i,:]
        F⁺[i] = trapz(ν, Mᵢ⁺)
        F⁻[i] = trapz(ν, Mᵢ⁻)
    end
    return nothing
end