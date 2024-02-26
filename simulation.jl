include("active_sphere.jl")

N = 300 #number of particles
n_steps = 10000 #number of time steps
sigma = 0.1   # radius of the particle
R =  8.7* sigma # radius of the sphere
r_cutoff = 2.4 * sigma # range of interaction (XY)


const k = 1.0 #elastic constant
const gamma = 1.0 #friction coefficient
const mu = 1.0/gamma
const tau = 1.0/(mu*k)  # some time scale

const J = 1*tau # J: XY type model
dt = 0.001 * tau # time step size
v_0 = 1.0*sigma/tau #
v_r = 0.0/tau #noise strength or rotational diffusion constant

rt = Vector{Vector{SVector{3,Float64}}}()
nt = Vector{Vector{SVector{3,Float64}}}()

system, n = init_system(N=N)
#print(n)
r = system.positions

for step in 1:n_steps
    #print(step)
    #print('\n')

    dr, dn = calulate_changes(system, r, n)

    for i in eachindex(dn)
        n[i] += dn[i]*dt .+ sqrt(v_r*dt)*randn()
        n[i] = normalize(n[i])
        r[i] += dr[i]*dt
        system.positions[i] = r[i]
    end
    
    push!(rt, copy(r))
    push!(nt, copy(n))

end