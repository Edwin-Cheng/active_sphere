include("active_sphere.jl")
using Plots

spf = 100 #steps per frame
frames = 100

function snapshot(rt, nt, t)
    x,y,z = get_3dvector_coms(rt[t])
    nx,ny,nz = get_3dvector_coms(nt[t])    
    arrow3d!(x,y,z,nx,ny,nz)
end

#get trajectory of particle i
function get_traj_i(rt, nt, i)
    rt_i = Vector{SVector{3,Float64}}()
    nt_i = Vector{SVector{3,Float64}}()
    for t in eachindex(rt)
        push!(rt_i, rt[t][i])
        push!(nt_i, nt[t][i])
    end
    return rt_i, nt_i
end

function get_ps(r,v) #order parameter 
    sum = [0.0, 0.0, 0.0]
    for i in eachindex(r)
        sum += cross(r[i], v[i])
    end
    return norm(sum)/N/R/v_0
end

print(get_ps(r,n))

plt = plot3d(
1,
xlim = (-R-1, R+1),
ylim = (-R-1, R+1),
zlim = (-R-1, R+1),
legend = false,
)
snapshot(rt,nt,n_steps)
display(plt)

gif=true
if gif == true 
    @gif for i in 1:frames
        plt = plot3d(
        1,
        xlim = (-R-1, R+1),
        ylim = (-R-1, R+1),
        zlim = (-R-1, R+1),
        legend = false,
        )
        snapshot(rt,nt,i*spf)
    end
end

