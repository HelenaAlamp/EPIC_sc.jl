struct SC_lens <: AbstractSpaceCharge
    optics::AbstractOptics4D
    ds::Float64   #ring arc length
    SC_lens(optics, ds) = new(optics, ds)
end

using DelimitedFiles

function track!(σx, σy, σz, temp1, bbSC::BunchedBeam, factorSC::Float64)
#function track!(dist::AbstractVector{ps6d{T}}, σx, σy, σz, temp1, factorSC::Float64) where T

    fieldvec_thread=[MVector{3}(0.0, 0.0, 0.0)  for j = 1:Threads.nthreads()]

    @inbounds Threads.@threads :static for j in eachindex(bbSC.dist.x)

        Bassetti_Erskine!(fieldvec_thread[Threads.threadid()], bbSC.dist.x[j], bbSC.dist.y[j], σx, σy)
        
        # gaussian function for particle distribution, lambda_z  (bbSC.dist.z[j]= temp1[j])
        temp1[j] = 1.0 /(σz*sqrt(2*π))*exp((-0.5) * bbSC.dist.z[j]^2 /σz^2)

        # delta p_/p_0 = ≈ xe-5/e-6
        
        bbSC.dist.px[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][1]
        bbSC.dist.py[j] += factorSC* temp1[j]* fieldvec_thread[Threads.threadid()][2]


    end

end


function SC_kick!(SC::SC_lens, bbSC::BunchedBeam)
    factorSC = 2*bbSC.particle.classrad0/bbSC.beta^2/bbSC.gamma^3*bbSC.num_particle*SC.ds  #if nSC=10, ≈3e-12
    
    get_emittance!(bbSC)
    betax = SC.optics.optics_x.beta
    betay = SC.optics.optics_y.beta
    
    σz = bbSC.beamsize[5]
    σx = sqrt(bbSC.emittance[1] * betax)  # beamsize x at SC_point
    σy = sqrt(bbSC.emittance[2] * betay) #beamsize y at SC_point
    
    track!(σx, σy, σz, bbSC.temp1, bbSC, factorSC)

end


#factor_tune = 2*bbSC.particle.classrad0*bbSC.num_particle/bbSC.beta^2/bbSC.gamma^3/(2*pi)^1.5/σz
#Δνx = factor_tune* betax * 3800 /(σx*(σy+σy))