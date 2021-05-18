#computes "exact" real part of faddeyeva
wofz(z::Complex)::Complex = erfcx(-im*z)

#grid characteristics for tests
xa = [-500, -200, -10, -500]
xb = [500, 200, 10, 500]
nx = [40001, 40001, 40001, 1000]
ya = [1e-5, 1e-20, 1e-5, 1e-20]
yb = [1e5, 1e4, 1e5, 1e5]
ny = [71, 71, 71, 25000]

for i = 1:4
    #make a grid
    x = collect(LinRange(xa[i], xb[i], nx[i]))
    y = 10 .^ collect(LinRange(log10(ya[i]), log10(yb[i]), ny[i]))
    X = ones(ny[i])*x'
    Y = y*ones(nx[i])'
    #complex argument to faddeyeva function
    Z = X .+ im*Y
    #whole complex result
    W = wofz.(Z)
    F = faddeyeva.(Z)
    relerr = @. abs(real(W) - real(F))/abs(real(W))
    @test maximum(relerr) < 1e-4
    relerr = @. abs(imag(W) - imag(F))/abs(imag(W))
    relerr = relerr[@. !isnan(relerr)]
    @test maximum(relerr) < 1e-4
    #real part only
    W = real.(W)
    F = faddeyeva.(X, Y)
    relerr = @. abs(W - F)/abs(W)
    @test maximum(relerr) < 1e-4
end
