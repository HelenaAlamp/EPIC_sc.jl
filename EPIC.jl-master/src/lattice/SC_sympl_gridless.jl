
#for test
#calculate generalized perveance / current #total charge for PROTON: 1.1023e-8
function calculate_K(bb::BunchedBeam, σz)

    K = 2*bb.particle.classrad0*bb.num_particle/(bb.beta^3*bb.gamma^3*σz)
    
    return K
end

############# 3 + 1 for cycle #############

function space_charge_gl(bb::BunchedBeam, K, Nl, Nm, a, b, Np, dt)

    gamma2lm = zeros(Nl, Nm)
    philm = zeros(Nl, Nm)
    temp = 0.0
    factor = 16*pi*K*dt/Np/a/b

    #x = zeros(Np)
    #y = zeros(Np)

    for i in 1:Nl
        al = i * pi / a
        for j in 1:Nm
            bm = j * pi / b

            @inbounds for k in 1:Np
                bb.dist.x[k] -= minimum(bb.dist.x)
                bb.dist.y[k] -= minimum(bb.dist.y)
                temp += sin(al*bb.dist.x[k])*sin(bm*bb.dist.y[k])

                #x[k] = bb.dist.x[k] - minimum(bb.dist.x)
                #y[k] = bb.dist.y[k] - minimum(bb.dist.y)
                #temp += sin(al*x[k])*sin(bm*y[k])

            end
            gamma2lm[i,j] = al^2 + bm^2
            philm[i,j] = factor*temp/gamma2lm[i,j]

        end
        
    end

    for i in 1:Nl
        al = i * pi / a
        for j in 1:Nm
            bm = j * pi / b

            @inbounds for l in 1:Np
                bb.dist.x[l] -= minimum(bb.dist.x)
                bb.dist.y[l] -= minimum(bb.dist.y)

                bb.dist.px[l] -= philm[i,j]* al*cos(al * bb.dist.x[l])*sin(bm * bb.dist.y[l])
                bb.dist.py[l] -= philm[i,j]* bm*sin(al * bb.dist.x[l])*cos(bm * bb.dist.y[l])

                #bb.dist.px[l] -= philm[i,j]* al*cos(al * x[l])*sin(bm * y[l])
                #bb.dist.py[l] -= philm[i,j]* bm*sin(al * x[l])*cos(bm * y[l])

            end

        end
    end

end



############# 3 + 1 for cycle #############


mutable struct SPACECHARGE <: AbstractSpaceCharge #problem in abstypes?   # spectral space charge
    # this element is treated as an integrated effect of space charge over a length of effective_len
    effective_len::Float64
    Nl::Int64
    Nm::Int64
    a::Float64
    b::Float64
    eletype::String

    function SPACECHARGE(effective_len::Float64, Nl::Int64, 
                        Nm::Int64, a::Float64, b::Float64)
        new(effective_len, Nl, Nm, a, b)
    end
end


function SC_gl_track!(ele::SPACECHARGE, bb::BunchedBeam, num_particles::Int64)
   
    get_emittance!(bb)
    σz = bb.beamsize[5]

    K = calculate_K(bb, σz)
    
    space_charge_gl(bb, K, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, ele.effective_len)
    #space_charge_gl_P(bb, ele.Nl, ele.Nm, ele.a, ele.b, num_particles, dt, lost_flags) #paralellization
    
    #return nothing
end


########################### test ###########################

###### test calculate_K ######
#=
num_particles=1000
beam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-3])
K = calculate_K(beam,0.006)
beam.dist.x

###### for test ######
#generate beam distribution
opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0) 
mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
αc=1.5e-3
lmap = LongitudinalRFMap(αc, mainRF)
initilize_6DGaussiandist!(beam, opIPp, lmap) #vector 6D phase space

Nl = 15
Nm = 15
#dx = 3e-5
#dy = 3e-5
a = 10e-3 #beam pipe 10mm
b = 10e-3
Np = num_particles
dt = 0.5 #3800/5/3e8 #ring division - time length

sc = SPACECHARGE(dt, Nl, Nm, a, b)

space_charge_gl(beam, K, Nl, Nm, a, b, Np, dt)

SC_gl_track!(sc, beam, Np)

using Plots
plot(beam.dist.x, beam.dist.px, seriestype=:scatter, markersize=2)
=#
