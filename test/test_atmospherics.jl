@test 8425 <scaleheight(9.8, 0.029, 288) < 8426

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Tstrat=200)
@test 15010 < altitude(1e4, 1e5, 9.8, Γ, (T,P)->0.029) < 15011

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Ptropo=2e4)
H = Hydrostatic(1e5, 0.1, 9.8, Γ, (T,P)->0.029)
@test 88625 < H(1e3) < 88626

Γ = DryAdiabat(288, 1e5, 1e3, 0.029, Ptropo=1e4)
H = Hydrostatic(1e5, 0.1, 9.8, Γ, (T,P)->0.029)
@test 24226 < altitude(H, 1e3) < 24227

Γ = MoistAdiabat(220, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Tstrat=200)
@test 70719 < tropopause(Γ)[2] < 70720

Γ = MoistAdiabat(220, 1e5, 1040, 1996, 0.029, 0.018, 2.3e6, psatH2O, Ptropo=1e4)
@test 116 < tropopause(Γ)[1] < 117
