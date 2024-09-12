
using Revise, EPIC
using Test
using BenchmarkTools
using Plots
using DelimitedFiles
using Serialization
using StructArrays
using FFTW


theme(:ggplot2)

struct record
    cx::Float64
    cy::Float64
    cz::Float64
    sx::Float64
    sy::Float64
    sz::Float64
    ex::Float64
    ey::Float64
    ez::Float64
    lumi::Float64
end

turns = 3
num_particles = 1000 #num_macro
nSC = 10

# Create the proton beam at IP

begin
    pbeam = BunchedBeam(PROTON, 0.688e11, 275e9, num_particles, [11.3e-9, 1e-9, 3.7e-4])
    opIPp = optics4DUC(0.8, 0.0, 0.072, 0.0) 
    mainRF = AccelCavity(591e6, 15.8e6, 7560.0, π)
    αc=1.5e-3
    lmap = LongitudinalRFMap(αc, mainRF)
    initilize_6DGaussiandist!(pbeam, opIPp, lmap)

    # define crab cavity
    overcrab=1.0
    crab_ratio=0.333
    pcrabu = easyCrabCavity(197.0e6, overcrab*12.5e-3*(1+crab_ratio))
    pcrabu2nd = easyCrabCavity(197.0e6*2.0, -overcrab*12.5e-3*crab_ratio)
    
    pcrabd = easyCrabCavity(197.0e6, -overcrab*12.5e-3*(1+crab_ratio))
    pcrabd2nd = easyCrabCavity(197.0e6*2.0, overcrab*12.5e-3*crab_ratio)

    # define Lorentz boost
    lb = LorentzBoost(12.5e-3)
    invlb = InvLorentzBoost(12.5e-3)

    # define oneturnmap
    tunex=0.228
    tuney=0.21
    oneturn = TransferMap4DChrom(opIPp, tunex, tuney, 1.0, 1.0)
    opIPe = optics4DUC(0.45, 0.0, 0.056, 0.0)
    estrong = StrongGaussianBeam(ELECTRON, 1.72e11, 10e9,  opIPe, [95e-6, 8.5e-6, 0.007], 5)
    #distz=open(deserialize, "ebeam_dist_z_wake3x.jls")
    
    #initilize_zslice!(estrong, distz, :evennpar, 11)
    #crab_crossing_setup!(estrong, 12.5e-3)
    
    initilize_zslice!(estrong, :gaussian, :evennpar, 5.0)
    ilb = InvLorentzBoost(12.5e-3)
    #track!(ebeam, ilb)
    ezcrab2 = easyCrabCavity(394.0e6, 12.5e-3, π) #wrong, crabcrossing setup in new EPIC update

    #define SC_lens
    ds=3833.8451/(nSC+1)
    #opSC = optics4DUC(0.5, 0.0, 0.5, 0.0) #random number
    opSC = optics4DUC(1.0, 0.0, 1.0, 0.0)
    phi_advx = LinRange(0, 29.228, nSC+1)
    phi_advy = LinRange(0, 30.210, nSC+1)
    SC = SC_lens(opSC, ds)
    
    # useful for SC_kick
    TM1 = Array{TransferMap4D, 1}(undef, nSC)
    TM1_inv = Array{Inverse_TransferMap4D, 1}(undef, nSC)

    particles_turns = zeros(Float64, num_particles, 6, turns)
    particles_turns_SC = zeros(Float64, num_particles, 6, turns)
    delta_px = zeros(Float64, num_particles)

    records = StructArray{record}(undef, turns)

    dt = 3800/5 #0.5
    sc = SPACECHARGE(dt, dt, 15, 15, 0.01, 0.01) #or 3800/5/299792458=2.5e-6

    lumi=0.0
end


@benchmark begin

    for i in 1:turns
        track!(pbeam, oneturn)
        track!(pbeam, lmap)
        track!(pbeam, pcrabu)
        track!(pbeam, pcrabu2nd)
        track!(pbeam, lb)
        #lumi=track!(pbeam, estrong)
        track!(pbeam, invlb)
        track!(pbeam, pcrabd)
        track!(pbeam, pcrabd2nd)
        
    
        particles_turns[:, 1, i] .= pbeam.dist.x
        particles_turns[:, 2, i] .= pbeam.dist.px
        particles_turns[:, 3, i] .= pbeam.dist.y
        particles_turns[:, 4, i] .= pbeam.dist.py
        particles_turns[:, 5, i] .= pbeam.dist.z
        particles_turns[:, 6, i] .= pbeam.dist.dp    
        
        #=
        for j in 1:(nSC)
            
            TM1[j] = TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
            track!(pbeam, TM1[j])
            
            SC_kick!(SC, pbeam)
            
            TM1_inv[j] = Inverse_TransferMap4D(opIPp, opSC, phi_advx[j], phi_advy[j])
            track!(pbeam, TM1_inv[j])

        end =#

        
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)
        SC_gl_track!(sc, pbeam, num_particles)

    
        particles_turns_SC[:, 1, i] .= pbeam.dist.x
        particles_turns_SC[:, 2, i] .= pbeam.dist.px
        particles_turns_SC[:, 3, i] .= pbeam.dist.y
        particles_turns_SC[:, 4, i] .= pbeam.dist.py
        particles_turns_SC[:, 5, i] .= pbeam.dist.z
        particles_turns_SC[:, 6, i] .= pbeam.dist.dp

        delta_px =  particles_turns[:,2,i] .- particles_turns_SC[:,2,i]
    
        get_emittance!(pbeam)
        records[i]=record(pbeam.centroid[1], pbeam.centroid[3], pbeam.centroid[5],
                            pbeam.beamsize[1], pbeam.beamsize[3], pbeam.beamsize[5],
                            pbeam.emittance[1], pbeam.emittance[2], pbeam.emittance[3], lumi)


    end
    

end

particles_turns
particles_turns_SC
delta_px


plot(particles_turns[:, 1, turns], delta_px, seriestype=:scatter, markersize=1)    

######### tune diagram #########
using PyCall
@pyimport nafflib

half=turns ÷ 2

tunex_1=nafflib.multiparticle_tunes(particles_turns_SC[:, 1, 1:half ])
tunex_2=nafflib.multiparticle_tunes(particles_turns_SC[:, 1, half:end])
tuney_1=nafflib.multiparticle_tunes(particles_turns_SC[:, 3, 1:half])
tuney_2=nafflib.multiparticle_tunes(particles_turns_SC[:, 3, half:end])

diff_tunex = sqrt.((tunex_1 .- tunex_2) .^2  .+ (tuney_1 .- tuney_2) .^2 )

maximum(diff_tunex)


scatter(tunex_1, tuney_1, marker_z = log10.(diff_tunex .+ 1e-15), markersize = 0.5,  color = :jet, clim=(-10,-2), 
aspect_ratio=:equal, legend=:topleft, xlabel="Horizontal Tune", ylabel= "Vertical Tune", label="nSC = 10", dpi=300)

    #,size=(600,500)#, xlim=(0.228, 0.243), ylim=(0.21, 0.225), xlabel="Horizontal Tune", ylabel="Vertical Tune", label="something", dpi=300)


##### FFT #####
using LaTeXStrings

plot(records.cz, label = L"$σ_z$ centroid")
plot!(records.sz, label = L"$σ_z$")
plot(records.ex, label = L"$ε_x$", xlabel="# Turns")
plot(records.sx, label = L"$σ_x$", xlabel="# Turns")


fftx=fft(@view records.ex[251:end]) #N/2 + 1 to N/2
ffty=fft(@view records.cy[251:end])
fftz=fft(@view records.cz[251:end])
fftf=fftfreq(500, 1.0)

plot(@view(fftf[1:turns÷2]), abs.(@view fftx[1:turns÷2]), xlim=(0,0.25))
xlabel!(L"$σ_x$")
plot(@view(fftf[1:turns÷2]), abs.(@view ffty[1:turns÷2]), xlim=(0,0.25))

plot(abs.(@view(ffty[1:turns÷2])), abs.(@view fftx[1:turns÷2]), seriestype=:scatter, markersize=1, xlim=(0,2.5e-6))


###### plot #######

plot(pbeam.dist.z, pbeam.dist.dp, seriestype=:scatter, markersize=1)
plot(pbeam.dist.x, pbeam.dist.px, seriestype=:scatter, markersize=1)

#save and plot
open("pbeam_records_tracking.jls", "w") do f
    serialize(f, records)
end

records_1x = open(deserialize, "pbeam_records_tracking.jls")
plot(records_1x.sx, label="eps_x")
records_1y = open(deserialize, "pbeam_records_tracking.jls")
plot(records_1x.sy, label="eps_y")
records_1x = open(deserialize, "pbeam_records_tracking.jls")
plot(records_1x.sz, label="eps_x")


#calculate average:
#avg_emix = mean(records_1x.ex)
#avg_emiy = mean(records_1y.ey)
avg_sx = mean(records_1x.sx)
avg_sy = mean(records_1y.sy)
avg_sz = mean(records_1y.sz)

factor_tune = 2*pbeam.particle.classrad0*pbeam.num_particle/pbeam.beta^2/pbeam.gamma^3/(2*pi)^1.5/avg_sz
Δνx = factor_tune*3800.0 /(avg_sx*(avg_sx+avg_sy))*1.0
Δνy = factor_tune*3800.0 /(avg_sy*(avg_sx+avg_sy))*1.0



# Save the results into file
#=
open("pbeam_wake3x_tracking.jls", "w") do f
    serialize(f, records)
end

open("pbeam_wake3x_status.jls", "w") do f
    serialize(f, pbeam)
end
=#

# Read file and process it later
#=
records_1x=open(deserialize, "pbeam_wake1x_tracking.jls")
plot(records_1x.sy, label="cz")
records_3x=open(deserialize, "pbeam_wake3x_tracking.jls")
plot!(records_3x.sy, label="cz")
=#
