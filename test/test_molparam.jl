for i in 1:length(MOLPARAM)
    if length(MOLPARAM[i].I) > 1
        bcheb = MOLPARAM[i].hascheb
        ncheb = MOLPARAM[i].ncheb
        rlerr = MOLPARAM[i].maxrelerr
        coefs = MOLPARAM[i].cheb
        @test all(rlerr .<= 0.01)
        for j = 1:length(ncheb)
            @test ncheb[j] == length(coefs[j])
            if bcheb[j]
                @test !any(isnan.(coefs[j]))
            else
                @test length(coefs[j]) == 0
            end
        end
        @test sum(MOLPARAM[i].A) <= 1.001
    end
end
