@test scaleheight(9.8, 0.029, 288) ≈ 8425.634182125264

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Tstrat=200)
@test altitude(1e4, 1e5, 9.8, Γ, (T,P)->0.029) ≈ 15010.645624480181

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Ptropo=2e4)
H = Hydrostatic(1e5, 0.1, 9.8, Γ, (T,P)->0.029)
@test H(1e3) ≈ 88625.38177948173

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Ptropo=1e4)
H = Hydrostatic(1e5, 0.1, 9.8, Γ, (T,P)->0.029)
@test altitude(H, 1e3) ≈ 24226.80062521415

Γ = MoistAdiabat(220, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Tstrat=200)
@test tropopause(Γ)[2] ≈ 70721.89893432238

Γ = MoistAdiabat(220, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Ptropo=1e4)
@test tropopause(Γ)[1] ≈ 116.63675563727014
