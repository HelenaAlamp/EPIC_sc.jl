#electric field Ex

function Bassetti_Erskine_xgtySC!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64) # x size greater than y
    # Only positive y is valid for this function
    # for y<0, Ex = Ex, Ey = -Ey
    if y < 0.0
        Bassetti_Erskine_xgtySC!(res, x, -y, σx, σy)
        res[2] = -res[2]
        nothing
        return
    end
    termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
	res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/σx/σy
    nothing
end

#electric field Ex
function Bassetti_Erskine_ygtxSC!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64) # x size greater than y
    # Only negative x is valid for this function
    # for x>0, Ex = -Ex, Ey = Ey
    if x > 0.0
        Bassetti_Erskine_ygtxSC!(res, -x, y, σx, σy)
        res[1] = -res[1]
        nothing
        return
    end
    termexp=exp(-x*x/2/σx/σx-y*y/2/σy/σy)
	sqrtδσ2=sqrt(Complex(2*(σx*σx-σy*σy)))
	term1=erfcx(-1im*(x+1im*y)/sqrtδσ2)
	term2=erfcx(-1im*(x*σy/σx+1im*y*σx/σy)/sqrtδσ2)
	
	complex_e=-1im*2*sqrt(pi)/sqrtδσ2*(term1-termexp*term2)
    res[1]=real(complex_e)
    res[2]=-imag(complex_e)
    res[3]=termexp/2.0/π/σx/σy
	nothing
end

function Bassetti_ErskineSC!(res::AbstractVector, x::Float64, y::Float64, σx::Float64, σy::Float64)
    if σx > σy
        Bassetti_Erskine_xgtySC!(res, x, y, σx, σy)
        nothing
        return
    else
        Bassetti_Erskine_ygtxSC!(res, x, y, σx, σy)
        nothing
        return
    end
end

fieldvecSC = MVector{3,Float64}(0.0, 0.0, 0.0)
n=6 #ring divisions
ds = 3800/n

# initilize_sigma_x,y

σx = Vector{Float64}(undef, n)
σy = Vector{Float64}(undef, n)
σx_avg = Vector{Float64}(undef, n)
σy_avg = Vector{Float64}(undef, n)


function SC_point!(bbSC::BunchedBeam) where T


    for i in 1:n-1
        get_emittance!(bbSC)
        σz = bbSC.beamsize[5]
        num_macro = length(bbSC.dist.x)
        λz = Vector{Float64}(undef, num_macro)
        factor1 = Vector{Float64}(undef, num_macro)

        σx[i] = bbSC.beamsize[1]
        σx[i+1] = bbSC.beamsize[1]
        σy[i] = bbSC.beamsize[3]
        σy[i+1] = bbSC.beamsize[3]

        # σx_avg[i] = (σx[i]-σx[i+1])/2
        σx_avg[i] = middle(σx[i], σx[i+1])
        σy_avg[i] = middle(σy[i], σy[i+1])


        @inbounds for j in 1:num_macro
            Bassetti_ErskineSC!(fieldvecSC, bbSC.dist.x[j], bbSC.dist.y[j], σx_avg[i], σy_avg[i]) #see line 90
            
            # gaussian function for part. distribution, lambda_z
            λz[j] = 1 /(σz*sqrt(2*pi))*exp((-1/2) * (bbSC.dist.z[j] /σz)^2)
            
            # lambda_z*p0
            factor1[j] = 2*bbSC.particle.classrad0 * λz[j] /(bbSC.beta^2* bbSC.gamma^3)

            # delta p_/p_0
            bbSC.dist.px[j] += factor1[j]* fieldvecSC[1]* ds
            bbSC.dist.py[j] += factor1[j]* fieldvecSC[2]* ds

        
        end
        
    end
   
end

############# SC_ point with check on kx,ky and tune shift #############
#from line 68
#=

using Plots 
# initilize parameters for beta avg
emix = Vector{Float64}(undef, n)
emiy = Vector{Float64}(undef, n)
betax = Vector{Float64}(undef, n)
betay = Vector{Float64}(undef, n)
betax_avg = Vector{Float64}(undef, n)
betay_avg = Vector{Float64}(undef, n)

begin
    get_emittance!(bbSC)
    σz = bbSC.beamsize[5]
    num_macro = length(bbSC.dist.x)
    Δνx = Vector{Float64}(undef, num_macro)
    Δνy = Vector{Float64}(undef, num_macro)
end

function SC_point!(bbSC::BunchedBeam) where T

    for i in 1:n-1
        get_emittance!(bbSC)
        σz = bbSC.beamsize[5]
        num_macro = length(bbSC.dist.x)
        λz = Vector{Float64}(undef, num_macro)
        factor1 = Vector{Float64}(undef, num_macro)
        kx = Vector{Float64}(undef, num_macro)
        ky = Vector{Float64}(undef, num_macro)
        Δνx = Vector{Float64}(undef, num_macro)
        Δνy = Vector{Float64}(undef, num_macro)

        σx[i] = bbSC.beamsize[1]
        σx[i+1] = bbSC.beamsize[1]
        σy[i] = bbSC.beamsize[3]
        σy[i+1] = bbSC.beamsize[3]

        # σx_avg[i] = (σx[i]-σx[i+1])/2
        σx_avg[i] = middle(σx[i], σx[i+1])
        σy_avg[i] = middle(σy[i], σy[i+1])

        # emittance for beta and beta_avg

        emix[i] = bbSC.emittance[1]
        emix[i+1] = bbSC.emittance[1]

        emiy[i] = bbSC.emittance[2]
        emiy[i+1] = bbSC.emittance[2]

        betax[i] = σx[i]^2 / emix[i]
        betax[i+1] = σx[i+1]^2 / emix[i+1]
        
        betay[i] = σy[i]^2 / emiy[i]
        betay[i+1] = σy[i+1]^2 / emiy[i+1]
        
        betax_avg[i] = middle(betax[i],betax[i+1])
        betay_avg[i] = middle(betay[i],betay[i+1])


        @inbounds for j in 1:num_macro
            Bassetti_ErskineSC!(fieldvecSC, bbSC.dist.x[j], bbSC.dist.y[j], σx_avg[i], σy_avg[i]) #see line 90
            
            # gaussian function for part. distribution, lambda_z
            λz[j] = 1 /(σz*sqrt(2*pi))*exp((-1/2) * (bbSC.dist.z[j] /σz)^2)
            
            # lambda_z*p0
            factor1[j] = 2*bbSC.particle.classrad0 * λz[j] /(bbSC.beta^2* bbSC.gamma^3)

            # delta p_/p_0 = gb.charge[j] charge for BunchedBeam -> particle (ParticleType)
            bbSC.dist.px[j] += factor1[j]* fieldvecSC[1]* ds
            bbSC.dist.py[j] += factor1[j]* fieldvecSC[2]* ds


            kx[j] =  (2*bbSC.particle.classrad0* λz[j]) / (bbSC.beta^2* bbSC.gamma^3* σx_avg[i]* (σx_avg[i]+σy_avg[i]))
            #Δνx[j] = 1/(2*pi) * kx[j] * betax * ds

            push!(Δνx, 1/(2*pi)* kx[j]* betax* ds)

            ky[j] =  (2*bbSC.particle.classrad0* λz[j]) / (bbSC.beta^2* bbSC.gamma^3*σy_avg[i]*(σx_avg[i]+σy_avg[i]))
            #Δνy[j] = 1/(2*pi) * ky[j] * betay * ds

            push!(Δνy, 1/(2*pi)* ky[j]* betay* ds)
        end
        
    end
   
    plot(Δνx,Δνy,seriestype=:scatter, markersize=1)
end
=#

