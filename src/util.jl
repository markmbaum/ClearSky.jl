#-------------------------------------------------------------------------------
#log coordinates are handy

#upward calculations
P2ω(P) = -log(P)
ω2P(ω) = exp(-ω)
P2ω(Pₛ, Pₜ) = P2ω(Pₛ), P2ω(Pₜ)

#downward calculations
P2ι(P) = log(P)
ι2P(ι) = exp(ι)
P2ι(Pₜ, Pₛ) = P2ι(Pₜ), P2ι(Pₛ)

#-------------------------------------------------------------------------------

export trapz
function trapz(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "vectors must be equal length"
    s = zero(eltype(y))
    for i = 1:length(x) - 1
        @inbounds s += (x[i+1] - x[i])*(y[i] + y[i+1])/2
    end
    return s
end

export meshgrid
function meshgrid(x::AbstractVector, y::AbstractVector)
    X = x' .* ones(length(y))
    Y = y .* ones(length(x))'
    return X, Y
end

export logrange
function logrange(a, b, N::Int=101, γ::Real=1)
    ((10 .^ LinRange(0, γ, N)) .- 1)*(b - a)/(10^γ - 1) .+ a
end

# θ: latitude [-π/2,π/2]
# ϕ: longitude [0,2π]
# f(θ, ϕ)
# approximates ∫∫ f cos(θ) dθ dϕ , with θ∈[-π/2,π/2] and ϕ∈[0,2π]
export shellintegral
function shellintegral(f::Function, param=nothing; nθ::Int=360, nϕ::Int=720)
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

#function epsdiff(x, )

export insertdiff
function insertdiff(x::AbstractVector, ϵ::Real=1e-3)
    n = length(x)
    #output version of x with additional coordinates
    ξ = zeros(eltype(x), 3n)
    #spacing for each original point in x
    δ = zeros(eltype(x), n)
    #first element
    δ[1] = ϵ*(x[2] - x[1])
    ξ[1] = x[1]
    ξ[2] = x[1] + δ[1]
    ξ[3] = x[1] + 2δ[1]
    #interior elements
    @inbounds for i ∈ 2:n-1
        δ[i] = ϵ*min(x[i] - x[i-1], x[i+1] - x[i])
        ξ[3i-2] = x[i] - δ[i]
        ξ[3i-1] = x[i]
        ξ[3i]   = x[i] + δ[i]
    end
    #last element
    δ[n] = ϵ*(x[n] - x[n-1])
    ξ[3n-2] = x[n] - 2δ[n]
    ξ[3n-1] = x[n] - δ[n]
    ξ[3n]   = x[n]
    return ξ, δ
end

export evaldiff!
function evaldiff!(∂, y, δ)::Nothing
    #check the dimensions line up
    @assert mod(length(y), 3) == 0
    n = length(y) ÷ 3
    @assert length(∂) == length(δ) == n
    #2nd order forward diff at first element
    ∂[1] = (-3y[1] + 4y[2] - y[3])/(2δ[1])
    #2nd order central diff at interior points
    for i ∈ 2:n-1
        @inbounds ∂[i] = (-y[3i-2] + y[3i])/(2δ[i]) 
    end
    #2nd order backward diff at last element
    ∂[n] = (y[3n-2] - 4y[3n-1] + 3y[3n])/(2δ[n]) 
    return nothing
end

export evaldiff
function evaldiff(y, δ)
    ∂ = similar(y, length(y) ÷ 3)
    evaldiff!(∂, y, δ)
    return ∂
end

#performing a dot product as if you only know one element of one of the vectors at a time
function dots(x, y)
    d = zero(y)
    for i = 1:length(x)-1
        d += x[i]
        d *= y[i]/y[i+1]
    end
    d += x[end]
    d *= y[end]
    return d
end

#------------------------------------------------------------------------------
# a couple of root finding methods

terminate(a, b, tol)::Bool = abs(a - b) < (tol + tol*abs(b)) ? true : false

function terminate(x₁, x₂, y₁, y₂, tol)::Bool
    terminate(x₁, x₂, tol) && terminate(y₁, y₂, tol) && return true
    return false
end

export regulafalsi
function regulafalsi(F, x₁, x₂, p=nothing; tol=1e-6)
    @assert x₁ != x₂ "starting points must not be identical"
    y₁ = F(x₁, p)
    y₁ == 0 && return x₁
    y₂ = F(x₂, p)
    y₂ == 0 && return x₂
    @assert sign(y₁) != sign(y₂) "regula falsi non-bracketing"
    yₘ = floatmax(y₂)
    yₚ = zero(y₂)
    n = 0
    while !terminate(x₁, x₂, yₚ, yₘ, tol) || (n < 2)
        #store previous evaluation
        yₚ = yₘ
        #approximate zero
        xₘ = x₁ - y₁*(x₂ - x₁)/(y₂ - y₁)
        yₘ = F(xₘ, p)
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

export secant
function secant(F, x₁, x₂, p=nothing; tol=1e-6)
    @assert x₁ != x₂ "starting points must not be identical"
    y₁ = F(x₁, p)
    y₁ == 0 && return x₁
    y₂ = F(x₂, p)
    y₂ == 0 && return x₂
    y₃ = floatmax(y₂)
    x₃ = zero(x₂)
    n = 0
    while !terminate(x₁, x₂, y₁, y₂, tol) || (n < 2)
        #approximate zero
        x₃ = x₁ - y₁*(x₂ - x₁)/(y₂ - y₁)
        y₃ = F(x₃, p)
        #swap values
        x₁, x₂ = x₂, x₃
        y₁, y₂ = y₂, y₃
        #count
        n += 1
    end
    return x₃
end