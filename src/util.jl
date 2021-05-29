export trapz

"""
    trapz(x, y)

Integrate a sorted group of coordinates using the composite [trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule).
"""
function trapz(x::AbstractVector, y::AbstractVector)::Float64
    @assert length(x) == length(y) "vectors must be equal length"
    s = 0.0
    for i = 1:length(x)-1
        s += (x[i+1] - x[i])*(y[i] + y[i+1])/2
    end
    return s
end

function meshgrid(x::AbstractVector, y::AbstractVector)
    X = x' .* ones(length(y))
    Y = y .* ones(length(x))'
    return X, Y
end

function logrange(a, b, N::Int=101, γ=1)::Vector{Float64}
    ((10 .^ LinRange(0, γ, N)) .- 1)*(b - a)/(10^γ - 1) .+ a
end

function falseposition(F::T, x₁, x₂)::Float64 where {T}
    y₁ = F(x₁)
    if y₁ == 0
        return x₁
    end
    y₂ = F(x₂)
    if y₂ == 0
        return x₂
    end
    @assert sign(y₁) != sign(y₂) "false position non-bracketing"
    yₘ = Inf
    yₚ = NaN
    n = 0
    while !(yₚ ≈ yₘ) || (n < 3)
        #store previous evaluation
        yₚ = yₘ
        #approximate zero
        xₘ = x₁ - y₁*(x₂ - x₁)/(y₂ - y₁)
        yₘ = F(xₘ)
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

#-------------------------------------------------------------------------------
#log coordinates are handy

#upward calculations
P2ω(P)::Float64 = -log(P)
ω2P(ω)::Float64 = exp(-ω)
P2ω(Pₛ, Pₜ)::Tuple{Float64,Float64} = P2ω(Pₛ), P2ω(Pₜ)

#downward coordinates
P2ι(P)::Float64 = log(P)
ι2P(ι)::Float64 = exp(ι)
P2ι(Pₜ, Pₛ)::Tuple{Float64,Float64} = P2ι(Pₜ), P2ι(Pₛ)
