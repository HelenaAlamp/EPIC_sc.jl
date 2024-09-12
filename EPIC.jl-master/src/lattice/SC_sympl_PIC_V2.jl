
#include("/home/alampres/EPIC_new.jl/EPIC.jl-master/src/abstypes.jl")

#calculate generalized perveance / current #total charge for PROTON:1.1023e-8
function calculate_K(bb::BunchedBeam, I::Float64)
    m0 = bb.particle.mass * 1.782662e-36 # kg
    charge_e = 1.60217663e-19
    charge = abs(bb.particle.charge * charge_e) #bb.particle.charge => -+1
    K = charge * I / (2.0 * pi *  8.854187817e-12 * 4e-7*pi * 2.99792458^3 * bb.beta^3 * bb.gamma^3)
    return K
end

#I = qe N/t = qe N betac/l = 8.7023e-4 total current along 3800m ?
#test calculate_K
num_particles=100
beam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-2])
K = calculate_K(beam, 5e-3)



function shape_function(x, deltax) 
    if abs(x) <= deltax / 2.0
        return 3.0 / 4.0 - (x / deltax)^2
    elseif deltax / 2.0 < abs(x) <= 3.0 * deltax / 2.0
        return 0.5 * (1.5 - abs(x / deltax))^2
    else
        return 0.0
    end
end

function d_shape_function(x, deltax)
    if abs(x) <= deltax / 2.0
        return -2.0 * x / deltax^2
    elseif deltax / 2.0 < abs(x) <= 3.0 * deltax / 2.0 && x > 0
        return (-1.5 + x / deltax) / deltax
    elseif deltax / 2.0 < abs(x) <= 3.0 * deltax / 2.0 && x <= 0
        return (1.5 + x / deltax) / deltax
    else
        return 0.0
    end
end


# generate beam distribution
opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0)
mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
αc=1.5e-3
lmap = LongitudinalRFMap(αc, mainRF)
initilize_6DGaussiandist!(beam, opIPp, lmap) #vector 6D phase space


#define other parameters for symplect. meth. with grid
Nl = 15
Nm = 15
dx = 3e-5
dy = 3e-5
a = 10e-3 #beam pipe 10mm
b = 10e-3
Np = num_particles


function calculate_philm(bb::BunchedBeam, Nl, Nm, dx, dy, a, b, Np)
    rholm = zeros(Nl, Nm)
    gamma2lm = zeros(Nl, Nm)
    philm = zeros(Nl, Nm)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]

    for i in 1:Nl
        for j in 1:Nm
            for k in 1:Np
                grid_x = (i - 1) * dx
                grid_y = (j - 1) * dy
                rholm .+= shape_function(rinx[k] .- grid_x, dx) * shape_function(riny[k] .- grid_y, dy) #eq 27? #aggiungo io, ho messo il punto
            end
        end
    end
    rholm *= 1.0 / (dx * dy * Np) #constant factor in eq.27

    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i, j] = al^2 + bm^2
            philm[i, j] += 4.0 * pi * rholm[i, j] / gamma2lm[i, j] #eq 26
        end
    end

    return philm, gamma2lm, rholm

end

calculate_philm(beam, Nl, Nm, dx, dy, a, b, Np)
##########fine aggiungo io########




function calculate_philm_lf(bb::BunchedBeam, Nl, Nm, a, b, Np, lost_flags) #if lost_flags[k] == 1, ignore + go to  the next one
    philm = zeros(Nl, Nm)
    gamma2lm = zeros(Nl, Nm)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]

    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i, j] = al^2 + bm^2
            for k in 1:Np
                if lost_flags[k] == 1
                    continue
                end
                #philm[i, j] .+= sin(al * rin[(k-1)*6 + 1]) * sin(bm * rin[(k-1)*6 + 3]) / gamma2lm[i, j] #eq 28, integrals solved?
                philm[i, j] += sin(al * rinx[k]) * sin(bm * riny[k]) / gamma2lm[i, j] #aggiungo io, ho tolto il punto
            
            end
        end
    end
    philm .*= 4.0 * pi  * 4.0 / (a * b * Np) #constant in eq 28
    return philm, gamma2lm
end


#lost_flags is an array of Np elements of 0/1
lost_flags = zeros(Int, num_particles) #solo 0, così considera tutti gli elementi ???
calculate_philm_lf(beam, Nl, Nm, a, b, Np, lost_flags)



lost_flags = zeros(Float64, Np)

function calculate_philm_P(bb::BunchedBeam, Nl, Nm, dx, dy, a, b, Np, lost_flags) #with paralellization, same calculate_philm_lf
    philm = zeros(Nl, Nm)
    gamma2lm = zeros(Nl, Nm)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]


    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            gamma2lm[i, j] = al^2 + bm^2

            local_sum = zeros(Threads.nthreads())  
        
            Threads.@threads for k in 1:Np
                tid = Threads.threadid()
                if lost_flags[k] == 1
                    continue
                end
                #local_sum[tid] += sin(al * rin[(k-1)*6 + 1]) * sin(bm * rin[(k-1)*6 + 3]) / gamma2lm[i, j]
                local_sum[tid] += sin(al * rinx[k]) * sin(bm * riny[k]) / gamma2lm[i, j] #aggiunto io 
            
            end
            
            # Sum up all partial results from each thread
            philm[i, j] = sum(local_sum)
        end
    end
    philm .*= 4.0 * pi  * 4.0 / (a * b * Np)
    return philm, gamma2lm
end

calculate_philm_P(beam, Nl, Nm, dx, dy, a, b, Np, lost_flags)


#px_b = rin[:,3] #mi serve dopo 

function space_charge!(bb::BunchedBeam, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)
    philm, gamma2lm = calculate_philm_P(bb, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    term1 = zeros(Np)
    term2 = zeros(Np)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]
    rinpx = @view bb.dist.px[:,1]
    rinpy = @view bb.dist.py[:,1]

    for i in 1:Nl
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            for k in 1:Np #aggiungo io
                #term1 .+= (philm[i, j] * al .* cos.(al .* rin[1:6:end]) .* sin.(bm .* rin[3:6:end]))
                #term2 .+= (philm[i, j] * bm .* sin.(al .* rin[1:6:end]) .* cos.(bm .* rin[3:6:end]))
                term1 .+= (philm[i, j] * al .* cos.(al .* rinx[k]) .* sin.(bm .* riny[k]))
                term2 .+= (philm[i, j] * bm .* sin.(al .* rinx[k]) .* cos.(bm .* riny[k]))
            end
        end
    end

    #update px, py eq.35, ####### no dovrebbe essere K/2 ??
    #r_in[2:6:end] .-= (dt * K / 2.0) * term1
    #r_in[4:6:end] .-= (dt * K / 2.0) * term2
    rinpx .-= (dt * K / 2.0) .* term1 #aggiungo io
    rinpx .-= (dt * K / 2.0) .* term2 #aggiungo io
    
    println(rinpx)
end

dt=3800/10*3e8 #tau should be the path of the beam?
space_charge!(beam, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)

#delta_px = []
#delta_px = px_b .- rin[:,3]  #to see if there's the kick ?!



function space_charge_P!(bb::BunchedBeam, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)  #SC paralellization
    philm, gamma2lm = calculate_philm_P(bb, Nl, Nm, dx, dy, a, b, Np, lost_flags)
    term1 = zeros(Np)
    term2 = zeros(Np)

    rinx = @view bb.dist.x[:,1]
    riny = @view bb.dist.y[:,1]
    rinpx = @view bb.dist.px[:,1]
    rinpy = @view bb.dist.py[:,1]

    nthreads = Threads.nthreads()
    term1_thread = [zeros(Np) for _ in 1:nthreads]  
    term2_thread = [zeros(Np) for _ in 1:nthreads]

    Threads.@threads for i in 1:Nl   #non ho bisogno di j l ???
        tid = Threads.threadid()
        for j in 1:Nm
            al = i * pi / a
            bm = j * pi / b
            term1_thread[tid] .+= (philm[i, j] * al .* cos.(al .* rinx[1:1:end]) .* sin.(bm .* riny[1:1:end]))
            term2_thread[tid] .+= (philm[i, j] * bm .* sin.(al .* rinx[1:1:end]) .* cos.(bm .* riny[1:1:end]))
        end
    end

    # Combine results from all threads
    for t in 1:nthreads
        term1 .+= term1_thread[t]
        term2 .+= term2_thread[t]
    end
    rinpx[1:1:end] .-= (dt * K / 2.0) * term1
    rinpy[1:1:end] .-= (dt * K / 2.0) * term2

    println(rinpx)
end

space_charge_P!(beam, K, Nl, Nm, dx, dy, a, b, Np, dt, lost_flags)



mutable struct SPACECHARGE <: AbstractSpaceCharge # spectral space charge
    # this element is treated as an integrated effect of space charge over a length of effective_len
    name::String
    len::Float64
    effective_len::Float64
    Nl::Int64
    Nm::Int64
    a::Float64
    b::Float64
    eletype::String

    function SPACECHARGE(name::String = "SPACECHARGE", len::Float64 = 0.0, effective_len::Float64 = 0.0, Nl::Int64 = 15, 
                        Nm::Int64 = 15, a::Float64 = 10e-3, b::Float64 = 10e-3)
        new(name, len, effective_len, Nl, Nm, a, b, "SPACECHARGE")
    end
end



#function SC_track!(ele::SPACECHARGE, bb::BunchedBeam, num_particles::Int64, particles::Beam)
function SC_track!(ele::SPACECHARGE, bb::BunchedBeam, num_particles::Int64)
    # ele: SPACECHARGE
    # rin: 6-by-num_particles array
    # num_particles: number of particles

    lost_flags = zeros(Int, num_particles) #particles.lost_flag         #da aggiungere in BunchedBeam?
    I = 0.1                                #particles.current           #da aggiungere in BunchedBeam?

    dx = ele.a / ele.Nl
    dy = ele.b / ele.Nm
    # v = speed_of_light * sqrt(1.0 - 1.0 / particles.gamma^2)
    dt = ele.effective_len 
    K = calculate_K(bb, I)
    #m0 = particles.mass * 1.782662e-36 # kg
    # p0 = particles.gamma * particles.beta * m0 * speed_of_light

    space_charge!(bb, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end


sc = SPACECHARGE("SPACECHARGE", 1.0, 0.98, 15, 15, 0.01, 0.01)

#SC_track!(sc, beam, num_particles, particles::Beam)
SC_track!(sc, beam, num_particles)




#SC with paralellization
function pass_TPSA!(ele::SPACECHARGE, rin::Vector{CTPS{T, TPS_Dim, Max_TPS_Degree}}) where {T, TPS_Dim, Max_TPS_Degree}
    println("TPSA is not implemented for space charge yet.")
    return nothing
end

function pass_P!(ele::SPACECHARGE, rin::Array{Float64,1}, num_particles::Int64, particles::Beam)
    # println("Parallel computing is not implemented for space charge yet.")
    lost_flags = particles.lost_flag
    I = particles.current
    dx = ele.a / ele.Nl
    dy = ele.b / ele.Nm
    # v = speed_of_light * sqrt(1.0 - 1.0 / particles.gamma^2)
    dt = ele.effective_len
    K = calculate_K(particles, I)
    # convert the mass from eV to kg
    m0 = particles.mass * 1.782662e-36 # kg
    # p0 = particles.gamma * particles.beta * m0 * speed_of_light
    space_charge_P!(rin, K, ele.Nl, ele.Nm, dx, dy, ele.a, ele.b, num_particles, dt, lost_flags)
    return nothing
end
