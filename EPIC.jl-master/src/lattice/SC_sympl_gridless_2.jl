
#for testing
#include("/home/alampres/EPIC_new.jl/EPIC.jl-master/src/abstypes.jl")

#calculate generalized perveance / current #total charge for PROTON: 1.1023e-8
function calculate_K(bb::BunchedBeam, σz)

    K = 2*bb.particle.classrad0*bb.num_particle/(bb.beta^3*bb.gamma^3*σz)
    
    return K
end


###### test calculate_K ######
#=
num_particles=1000
beam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-3])
K = calculate_K(beam)
=#

###### for test ######
#generate beam distribution
#=
opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0) 
mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
αc=1.5e-3
lmap = LongitudinalRFMap(αc, mainRF)
initilize_6DGaussiandist!(beam, opIPp, lmap) #vector 6D phase space
=#


###### for test ######
#=
Nl = 15
Nm = 15
#dx = 3e-5
#dy = 3e-5
a = 10e-3 #beam pipe 10mm
b = 10e-3
Np = num_particles
dt = 3800/5/3e8 #ring division - time length
=#


function space_charge_gl(bb::BunchedBeam, K, Nl, Nm, a, b, Np, dt)

    gamma2lm = zeros(Nl, Nm)

    term1 = zeros(Np)
    #term2 = zeros(Np)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]
    #rinpx = @view bb.dist.px[:,1]
    #rinpy = @view bb.dist.py[:,1]

    factor = 16*pi*K*dt/Np/a/b #Np is number of macroparticles!!!
    
    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            for k in 1:Np
                for l in 1:Np
                    gamma2lm[i,j] = al^2 + bm^2
                    term1[l] += (al * sin(al * rinx[k]) * sin(bm * riny[k]) * cos(al * rinx[l]) * cos(bm * riny[l]))/gamma2lm[i,j]
                    #term2[l] += (bm * sin(al * rinx[k]) * sin(bm * riny[k]) * cos(al * rinx[l]) * cos(bm * riny[l]))/gamma2lm[i,j]
                    bb.dist.px[l] -= (factor) * term1[l]
                end
            end
        end
    end


    #bb.dist.px .-= (factor) * term1[l]
    #bb.dist.py .-= (factor) * term2
    
    
    #println(rinpx)

end

###### for test ######
#space_charge_gl(beam, Nl, Nm, a, b, Np, dt)


#add lost_flags before use
#=
function space_charge_gl_P(bb::BunchedBeam, Nl, Nm, a, b, Np, dt, lost_flags) #with paralellization, same calculate_philm_lf
    
    gamma2lm = zeros(Nl, Nm)

    term1 = zeros(Np)
    term2 = zeros(Np)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]
    rinpx = @view bb.dist.px[:,1]
    rinpy = @view bb.dist.py[:,1]

    factor = 16*pi*K*dt/Np/a/b

    nthreads = Threads.nthreads()
    term1_thread = [zeros(Np) for _ in 1:nthreads]  
    term2_thread = [zeros(Np) for _ in 1:nthreads]

    Threads.@threads for i in 1:Nl
        tid = Threads.threadid()
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i,j] = al^2 + bm^2
            term1_thread[tid] .+= (al .*sin.(al .* rinx[k]) .* sin.(bm .* riny[k]) .* cos.(al .* rinx[k]) .* cos.(bm .* riny[k]))./gamma2lm[i,j]
            term2_thread[tid] .+= (bm .*sin.(al .* rinx[k]) .* sin.(bm .* riny[k]) .* cos.(al .* rinx[k]) .* cos.(bm .* riny[k]))./gamma2lm[i,j]
        end
    end

    # Combine results from all threads
    for t in 1:nthreads
        term1 .+= term1_thread[t]
        term2 .+= term2_thread[t]
    end

    rinpx .-= factor * term1
    rinpy .-= factor * term2
    

    #println(rinpx)
end=#
 


mutable struct SPACECHARGE <: AbstractSpaceCharge #problem in abstypes?   # spectral space charge
    # this element is treated as an integrated effect of space charge over a length of effective_len
    len::Float64
    effective_len::Float64
    Nl::Int64
    Nm::Int64
    a::Float64
    b::Float64
    eletype::String

    function SPACECHARGE(len::Float64, effective_len::Float64, Nl::Int64, 
                        Nm::Int64, a::Float64, b::Float64)
        new(len, effective_len, Nl, Nm, a, b)
    end
end


function SC_gl_track!(ele::SPACECHARGE, bb::BunchedBeam, num_particles::Int64)
   
    get_emittance!(bb)
    σz = bb.beamsize[5]

    K = calculate_K(bb, σz)
    dt = ele.effective_len
    space_charge_gl(bb, K, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, dt)
    #space_charge_gl_P(bb, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, dt, lost_flags) #with paralellization
    
    return nothing
end


###### for test ######
#dt=3800/5/3e8
#sc = SPACECHARGE(dt, dt, Nl, Nm, a, b)
#SC_gl_track!(sc, beam, Np)