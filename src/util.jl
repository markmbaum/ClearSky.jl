#-------------------------------------------------------------------------------
#transformed coordinates are handy

#upward calculations
P2ω(P) = -sqrt(P)
ω2P(ω) = ω*ω
dωfac(P) = 2*sqrt(P)
P2ω(Pₛ, Pₜ) = P2ω(Pₛ), P2ω(Pₜ)

#downward calculations
P2ι(P) = sqrt(P)
ι2P(ι) = ι*ι
dιfac(P) = 2*sqrt(P)
P2ι(Pₜ, Pₛ) = P2ι(Pₜ), P2ι(Pₛ)

#-------------------------------------------------------------------------------

export pressuregrid
function pressuregrid(Pₜ, Pₛ, n)
    @assert Pₛ > Pₜ
    @assert n >= 3
    exp.(chebygrid(log(Pₜ), log(Pₛ), n))
end

export trapz
function trapz(x::AbstractVector, y::AbstractVector)
    @assert length(x) == length(y) "vectors must be equal length"
    s = zero(eltype(y))
    @inbounds for i ∈ 1:length(x) - 1
        s += (x[i+1] - x[i])*(y[i] + y[i+1])/2
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

#evaluates the 0th, 1st, and 2nd derivative of a parabola through three points, anywhere
function quaddiff(x, x₁, x₂, x₃, y₁, y₂, y₃)
    z₁ = x₁^2 - x₂^2
    z₂ = x₂^2 - x₃^2
    w = z₁/z₂
    b = (y₁ - y₂ - (y₂ - y₃)*w)/(x₁ - x₂ - (x₂ - x₃)*w)
    a = (y₂ - y₃ - b*(x₂ - x₃))/z₂
    c = y₁ - (a*x₁^2 + b*x₁)
    f = a*x^2 + b*x + c
    f′ = 2*a*x + b
    f′′ = 2*a
    return f, f′, f′′
end

function quaddiff(x, xₚ, yₚ)
    @assert length(xₚ) == length(yₚ) == 3
    @inbounds quaddiff(x, xₚ[1], xₚ[2], xₚ[3], yₚ[1], yₚ[2], yₚ[3])[2]
end

export deriv
function deriv(x, y, T)
    d = similar(x)
    d[1] = (y[2] - y[1])/(x[2] - x[1])
    for i in 2:length(d)-1
        #d₁ = (y[i] - y[i-1])/(x[i] - x[i-1])
        #d₂ = (y[i+1] - y[i])/(x[i+1] - x[i])
        #if i == 25
        #    d[i] = (d₁ + d₂)/2
        #elseif i > 25
        #    d[i] = d₁
        #else
        #    d[i] = d₂
        #end
        #if d₁*d₂ > 0
        #    d[i] = (d₁ + d₂)/2
        #else
        #    if T[i-1] < T[i+1]
        #        d[i] = d₁
        #    else
        #        d[i] = d₂
        #    end
        #end
        d[i] = quaddiff(x[i], view(x,i-1:i+1), view(y,i-1:i+1))
    end
    d[end] = (y[end] - y[end-1])/(x[end] - x[end-1])
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