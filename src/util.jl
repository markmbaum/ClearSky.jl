export trapz, logrange, meshgrid, shellintegral

"""
    trapz(x, y)

Integrate a sorted group of coordinates using the composite [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).
"""
function trapz(x::AbstractVector, y::AbstractVector)::Float64
    @assert length(x) == length(y) "vectors must be equal length"
    s = 0.0
    for i = 1:length(x) - 1
        @inbounds s += (x[i+1] - x[i])*(y[i] + y[i+1])/2
    end
    return s
end

function meshgrid(x::AbstractVector, y::AbstractVector)
    X = x' .* ones(length(y))
    Y = y .* ones(length(x))'
    return X, Y
end

function logrange(a, b, N::Int=101, γ::Real=1)::Vector{Float64}
    ((10 .^ LinRange(0, γ, N)) .- 1)*(b - a)/(10^γ - 1) .+ a
end

function interp(q::Float64,
                x::AbstractVector{Float64},
                y::AbstractVector{Float64})::Float64
    #find the proper cell
    i = findcell(q, x)
    #interpolate
    (q - x[i])*(y[i+1] - y[i])/(x[i+1] - x[i]) + y[i]
end

function falseposition(F::T,
                       x₁::Real,
                       x₂::Real,
                       param=nothing;
                       tol::Float64=1e-6)::Float64 where {T}
    y₁ = F(x₁, param)
    if y₁ == 0
        return x₁
    end
    y₂ = F(x₂, param)
    if y₂ == 0
        return x₂
    end
    @assert sign(y₁) != sign(y₂) "false position non-bracketing"
    yₘ = Inf
    yₚ = NaN
    n = 0
    while abs(yₚ - yₘ) > (tol + tol*abs(yₘ)) || (n < 2)
        #store previous evaluation
        yₚ = yₘ
        #approximate zero
        xₘ = x₁ - y₁*(x₂ - x₁)/(y₂ - y₁)
        yₘ = F(xₘ, param)
        #reduce the interval
        if y₁*yₘ > 0
            x₁ = xₘ
        else
            x₂ = xₘ
        end
        #count
        n += 1
    end
    return (x₁ + x₂)/2
end

# θ: latitude [-π/2,π/2]
# ϕ: longitude [0,2π]
# f(θ, ϕ)
# approximates ∫∫ f cos(θ) dθ dϕ , with θ∈[-π/2,π/2] and ϕ∈[0,2π]
function shellintegral(f::T,
                       param=nothing;
                       nθ::Int=360,
                       nϕ::Int=720)::Float64 where {T}

    I = 0.0
    Δθ = π/nθ
    Δϕ = 2π/nϕ
    #integrate
    θ = -π/2 + Δθ/2
    for θ ∈ LinRange(-π/2 + Δθ/2, π/2 - Δθ/2, nθ)
        for ϕ ∈ LinRange(Δϕ/2, 2π - Δϕ/2, nϕ)
            I += f(θ, ϕ, param)*cos(θ)*Δθ*Δϕ
        end
    end
    return I
end

function derivative!(dydx::AbstractVector{T},
                     x::AbstractVector{T},
                     y::AbstractVector{T})::Nothing where {T}
    @assert length(dydx) == length(x) == length(y)
    n = length(x)
    dydx[1] = (y[2] - y[1])/(x[2] - x[1])
    for i ∈ 2:n-1
        dydx[i] = ((y[i+1] - y[i])/(x[i+1] - x[i]) + (y[i] - y[i-1])/(x[i] - x[i-1]))/2
    end
    dydx[n] = (y[n] - y[n-1])/(x[n] - x[n-1])
    return nothing
end

#-------------------------------------------------------------------------------
#log coordinates are handy

#upward calculations
P2ω(P)::Float64 = -log(P)
ω2P(ω)::Float64 = exp(-ω)
P2ω(Pₛ, Pₜ)::NTuple{2,Float64} = P2ω(Pₛ), P2ω(Pₜ)

#downward calculations
P2ι(P)::Float64 = log(P)
ι2P(ι)::Float64 = exp(ι)
P2ι(Pₜ, Pₛ)::NTuple{2,Float64} = P2ι(Pₜ), P2ι(Pₛ)
