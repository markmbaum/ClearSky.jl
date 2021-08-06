"""
This module implements an efficient approximation for the Faddeyeva function
(sometimes called Faddeeva). Two methods are provided, one that is refined for
evaluation of the real part only and one for the whole complex result. The
approximation is due to Mofreh R. Zaghloul, as descriped in the paper:
* Mofreh R. Zaghloul. 2017. Algorithm 985: Simple, Efficient, and Relatively Accurate Approximation for the Evaluation of the Faddeyeva Function. ACM Trans. Math. Softw. 44, 2, Article 22 (August 2017), 9 pages. https://doi.org/10.1145/3119904
"""

module Faddeyeva

export faddeyeva

const θ = 1/√π

const α0 = 122.60793
const α1 = 214.38239
const α2 = 181.92853
const α3 = 93.15558
const α4 = 30.180142
const α5 = 5.9126262
const α6 = 1/√π

const β0 = 122.60793
const β1 = 352.73063
const β2 = 457.33448
const β3 = 348.70392
const β4 = 170.35400
const β5 = 53.992907
const β6 = 10.479857

const γ0 = 36183.31
const γ1 = 3321.99
const γ2 = 1540.787
const γ3 = 219.031
const γ4 = 35.7668
const γ5 = 1.320522
const γ6 = 1/√π

const λ0 = 32066.6
const λ1 = 24322.84
const λ2 = 9022.228
const λ3 = 2186.181
const λ4 = 364.2191
const λ5 = 61.57037
const λ6 = 1.841439

const s0 = 38000.0
const s1 = 256.0
const s2 = 62.0
const s3 = 30.0
const t3 = 1.0e-13
const s4 = 2.5
const t4 = 5.0e-9
const t5 = 0.072

#region 4: Laplace continued fractions, 4 convergents
function regionIV(z::Complex, x::Real, y::Real)::Complex
    z² = z^2
    (θ*(-y + im*x))*(z² - 2.5)/(z²*(z² - 3.0) + 0.75)
end

#region 5: Humlicek's w4 (Region IV), part a
function regionVa(z::Complex, x²::Real)::Complex
    z² = z^2
    r = γ0 + z²*(γ1 + z²*(γ2 + z²*(γ3 + z²*(γ4 + z²*(γ5 + z²*γ6)))))
    t = λ0 + z²*(λ1 + z²*(λ2 + z²*(λ3 + z²*(λ4 + z²*(λ5 + z²*(λ6 + z²))))))
    exp(-x²) + (im*z*r/t)
end

#region 5: Humlicek's w4 (Region IV), part b
function regionVb(z::Complex)::Complex
    z² = z^2
    r = γ0 + z²*(γ1 + z²*(γ2 + z²*(γ3 + z²*(γ4 + z²*(γ5 + z²*γ6)))))
    t = λ0 + z²*(λ1 + z²*(λ2 + z²*(λ3 + z²*(λ4 + z²*(λ5 + z²*(λ6 + z²))))))
    exp(-z²) + (im*z*r/t)
end

#region 6: Hui's p-6 Approximation
function regionVI(x::Real, y::Real)::Complex
    q = y - im*x
    r = α0 + q*(α1 + q*(α2 + q*(α3 + q*(α4 + q*(α5 + q*α6)))))
    t = β0 + q*(β1 + q*(β2 + q*(β3 + q*(β4 + q*(β5 + q*(β6 + q))))))
    r/t
end

"""
    faddeyeva(z::Complex)

Evaluate the [Faddeyeva function](https://en.wikipedia.org/wiki/Faddeeva_function) with a complex argument, returning a complex number.
"""
function faddeyeva(z::Complex)::Complex

    x = real(z)
    y = imag(z)
    x² = x*x
    y² = y*y
    s = x² + y²

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s0
        return (y + im*x)*θ/s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s1
        a = y*(0.5 + s)
        b = x*(s - 0.5)
        d = s^2 + (y² - x²) + 0.25
        return (a + im*b)*(θ/d)
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s2
        q = y² - x² + 1.5
        r = 4.0*x²*y²
        a = y*((q - 0.5)*q + r + x²)
        b = x*((q - 0.5)*q + r - y²)
        d = s*(q*q + r)
        return θ*(a + im*b)/d
    end

    #region 4: Laplace continued fractions, 4 convergents
    if s >= s3 && y² >= t3
        return regionIV(z, x, y)
    end

    #region 5: Humlicek's w4 (Region IV)
    if s > s4 && y² < t4
        return regionVa(z, x²)
    elseif s > s4 && y² < t5
        return regionVb(z)
    end

    #region 6: Hui's p-6 Approximation
    return regionVI(x, y)

end

#real arguments representint z = x + im*y, returning only the real part
"""
    faddeyeva(x::Real, y::Real)

Evaluate only the real part of the [Faddeyeva function](https://en.wikipedia.org/wiki/Faddeeva_function) from the complex argument ``x + iy``
"""
function faddeyeva(x::T, y::T)::Float64 where {T<:Real}

    x²::T = x*x
    y²::T = y*y
    s::T = x² + y²

    #region 1: Laplace continued fractions, 1 convergent
    if s >= s0
        return y*θ/s
    end

    #region 2: Laplace continued fractions, 2 convergents
    if s >= s1
        return y*(0.5 + s)*(θ/((s^2 + (y² - x²)) + 0.25))
    end

    #region 3: Laplace continued fractions, 3 convergents
    if s >= s2
        q = y² - x² + 1.5
        r = 4.0*x²*y²
        return θ*(y*((q - 0.5)*q + r + x²))/(s*(q*q + r))
    end

    if s >= s3 && y² >= t3
        return real(regionIV(x + im*y, x, y))
    end

    if s > s4 && y² < t4
        return real(regionVa(x + im*y, x²))
    elseif s > s4 && y² < t5
        return real(regionVb(x + im*y))
    end

    return real(regionVI(x, y))

end

end
