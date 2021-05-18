for i in 1:length(MOLPARAM)
    if length(MOLPARAM[i]) > 1
        bcheb = MOLPARAM[i][10]
        ncheb = MOLPARAM[i][11]
        rlerr = MOLPARAM[i][12]
        coefs = MOLPARAM[i][13]
        @test all(rlerr .<= 0.01)
        for j = 1:length(ncheb)
            @test ncheb[j] == length(coefs[j])
            if bcheb[j]
                @test !any(isnan.(coefs[j]))
            else
                @test length(coefs[j]) == 0
            end
        end
        @test sum(MOLPARAM[i][7]) <= 1.001
    end
end
