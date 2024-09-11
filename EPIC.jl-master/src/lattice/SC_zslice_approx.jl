
#initialize transverse beam size, 1> first slice and so on
#long code and wasting time

σx1 = Vector{Float64}(undef, n)
σy1 = Vector{Float64}(undef, n)


for i in 1:n-1
    get_emittance!(bbSC)
    σz = bbSC.beamsize[5]
    num_macro = length(bbSC.dist.x)
    λz = Vector{Float64}(undef, num_macro)
    factor1 = Vector{Float64}(undef, num_macro)

    σz = bbSC.beamsize[5]/n
    nzslice = 3
    X = bbSC.beamsize[1]
    Y = bbSC.beamsize[3]

    for i in 1:nzslice
        if σz < 0.06/3
            σx1[i] = max(pairwise(dist, X))
            σy1[i] = max(pairwise(dist, Y))
        end
        if σz < 0.06/3 && σz < 0.06/3*2 
            σx2[i] = max(pairwise(dist, X))
            σy2[i] = max(pairwise(dist, Y))
        end
    end

end

@inbounds for j in 1:num_macro
    Bassetti_ErskineSC!(fieldvecSC, bbSC.dist.x[j], bbSC.dist.y[j], σx1[i], σy1[i]) #see line 90
    #BE for each slice, with different σx / σy
    
    # gaussian function for part. distribution, lambda_z
    λz[j] = 1 /(σz*sqrt(2*pi))*exp((-1/2) * (bbSC.dist.z[j] /σz)^2)
    
    # lambda_z*p0
    factor1[j] = 2*bbSC.particle.classrad0 * λz[j] /(bbSC.beta^2* bbSC.gamma^3)

    # delta p_/p_0 = gb.charge[j] charge for BunchedBeam -> particle (ParticleType)
    bbSC.dist.px[j] += factor1[j]* fieldvecSC[1]* ds
    bbSC.dist.py[j] += factor1[j]* fieldvecSC[2]* ds


end


