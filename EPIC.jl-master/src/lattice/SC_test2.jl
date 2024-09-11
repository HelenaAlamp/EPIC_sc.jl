
struct SC_lens <: AbstractSpaceCharge
    optics::AbstractOptics4D  #useful for defining twiss parameters
    ds::Float64   #ring arc length, anche 3800/n ed n da definire
    SC_lens(optics, ds) = new(optics, ds)
end

#prova solo con abstract ps6d e non BunchedBeam
#dist::AbstractVector{ps6d{T}}
#I need of "beamsize", so BunchedBeam for get emittance and beamsize on z ??? 
function track!(temp1, temp2, temp3, temp4, temp5, bbSC::BunchedBeam, SC::SC_lens, factor1::Float64) where T 
    

    fieldvecSC = MVector{3,Float64}(0.0, 0.0, 0.0)
    get_emittance!(bbSC)
    σz = bbSC.beamsize[5]
    num_macro = bbSC.num_macro
    λz = Vector{Float64}(undef, num_macro)

    betax = SC.optics.optics_x.beta
    betay = SC.optics.optics_y.beta
    alphax = SC.optics.optics_x.alpha
    alphay = SC.optics.optics_y.alpha
    gammax = SC.optics.optics_x.gamma
    gammay = SC.optics.optics_y.gamma

    M = zeros(Float64, num_macro, 3)

    

    @inbounds for j in 1:num_macro

        temp1[j] = bbSC.dist.z[j]  # z location at SC_point
        temp4[j] = betax + gammax * temp1[j] * temp1[j] - 2.0 * alphax * temp1[j]  # beta x at SC_point
        temp5[j] = betay + gammay * temp1[j] * temp1[j] - 2.0 * alphay * temp1[j]   # beta y at SC_point
        temp2[j] = bbSC.beamsize[1] * sqrt(temp4[j] / betax)  # sqrt(bbSC.emittance.x[j] * betax)  # beamsize x at SC_point
        temp3[j] = bbSC.beamsize[2] * sqrt(temp5[j] / betay)    # sqrt(bbSC.emittance.y[j] * betay) #beamsize y at SC_point
        
        Bassetti_Erskine!(fieldvecSC, bbSC.dist.x[j], bbSC.dist.y[j], temp2[j], temp3[j]) #dovrebbe riprendere funzione in strogbb.jl

        # gaussian function for part. distribution, lambda_z  (bbSC.dist.z[j]= temp1[j])
        λz[j] = 1 /(σz*sqrt(2*pi))*exp((-1/2) * (bbSC.dist.z[j] /σz)^2)
        
        #factor1 = 2*bbSC.particle.classrad0/(bbSC.beta^2* bbSC.gamma^3)
        

        # delta p_/p_0
        bbSC.dist.px[j] += num_macro* factor1* λz[j]* fieldvecSC[1]* SC.ds
        bbSC.dist.py[j] += num_macro* factor1* λz[j]* fieldvecSC[2]* SC.ds
        M[j, 1] = temp2[j]
        M[j, 2] = temp3[j]
        M[j, 3] = temp2[j].-temp3[j]
    end

    println(M[:,3]) #-2.4156766915000354e-5 
    #println(fieldvecSC[1])  #303.1738534973828

end


function SC_kick!(SC::SC_lens, bbSC::BunchedBeam)
    factor1 = 2*bbSC.particle.classrad0/(bbSC.beta^2* bbSC.gamma^3)
    track!(bbSC.temp1, bbSC.temp2, bbSC.temp3, bbSC.temp4, bbSC.temp5, bbSC, SC,factor1)

end

#Analytical tune shift
#=emi_factorx = 1.0/bbSC.emittance[1]/(1.0+sqrt(bbSC.emittance[2]*betay/betax/bbSC.emittance[1]))
emi_factory = 1.0/bbSC.emittance[2]/(1.0+sqrt(bbSC.emittance[1]*betax/betay/bbSC.emittance[2]))

Δνx = []
Δνy = []
Δνx = push!(Δνx, factorSC/SC.ds*emi_factorx/(2*π))
Δνy = push!(Δνy, factorSC/SC.ds*emi_factory/(2*π))

#save in .txt
open("tunex_file.txt", "a") do io
    writedlm(io, Δνx, "\n")
end
open("tuney_file.txt", "a") do io
    writedlm(io, Δνy, "\n")
end=#
#println("Tune shift x,y:", Δνx, " " ,Δνy)

#Equations for tune shift ps
#Δνx = ro*N/beta^2/gamma^3/(2*pi)/Ex * 2/(1+sqrt(betay*Ey/betax/Ex))
#Δνy = ro*N/beta^2/gamma^3/(2*pi)/Ey * 2/(1+sqrt(betax*Ex/betay/Ey))


