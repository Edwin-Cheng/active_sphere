using Random
using LinearAlgebra
using CellListMap.PeriodicSystems, StaticArrays

mutable struct FAndNcross
    f::Vector{SVector{3,Float64}}
    n_ij::Vector{SVector{3,Float64}}
end
 
function ran_unit(N) #N random unit vectors
    v = rand(SVector{3,Float64},N)
    for i in eachindex(v)
        v[i] = normalize(v[i].-0.5)
    end
    return v
end

function get_3dvector_coms(vs)
    M = stack(vs,dims=1)
    return M[:,1],M[:,2],M[:,3]
end

function arrow3d!(x, y, z,  u, v, w; as=0.1, lc=:black, la=1, lw=0.4, scale=:identity)
    (as < 0) && (nv0 = -maximum(norm.(eachrow([u v w]))))
    for (x,y,z, u,v,w) in zip(x,y,z, u,v,w)
        nv = sqrt(u^2 + v^2 + w^2)
        v1, v2 = -[u,v,w]/nv, nullspace(adjoint([u,v,w]))[:,1]
        v4 = (3*v1 + v2)/3.1623  # sqrt(10) to get unit vector
        v5 = v4 - 2*(v4'*v2)*v2
        (as < 0) && (nv = nv0) 
        v4, v5 = -as*nv*v4, -as*nv*v5
        plot!([x,x+u], [y,y+v], [z,z+w], lc=lc, la=la, lw=lw, scale=scale, label=false)
        plot!([x+u,x+u-v5[1]], [y+v,y+v-v5[2]], [z+w,z+w-v5[3]], lc=lc, la=la, lw=lw, label=false)
        plot!([x+u,x+u-v4[1]], [y+v,y+v-v4[2]], [z+w,z+w-v4[3]], lc=lc, la=la, lw=lw, label=false)
    end
end

P_T(r, a) = a-dot(normalize(r),a)*normalize(r)

P_N(r, a) = dot(a,normalize(r))

#interface

import CellListMap.PeriodicSystems: copy_output
copy_output(x::FAndNcross) = FAndNcross(copy(x.f), copy(x.n_ij))

import CellListMap.PeriodicSystems: reset_output!
function reset_output!(output::FAndNcross)
    for i in eachindex(output.f)
        output.f[i] = SVector(0.0, 0.0, 0.0)
    end
    for i in eachindex(output.n_ij)
        output.n_ij[i] = SVector(0.0, 0.0, 0.0)
    end
    return output
end

import CellListMap.PeriodicSystems: reducer
function reducer(x::FAndNcross, y::FAndNcross)
    x.f .+= y.f
    x.n_ij .+= y.n_ij
    return FAndNcross(x.f, x.n_ij)
end


#initialization
function init_system(;N::Int=200)
    positions = ran_unit(N)*R
    n = ran_unit(N)
    for i in eachindex(n)
        n[i] = normalize(P_T(positions[i], n[i]))
    end

    system = PeriodicSystem(
        xpositions = positions,
        unitcell=[3*R,3*R,3*R], 
        cutoff = r_cutoff, 
        output = FAndNcross(similar(positions), similar(positions)),
        output_name = :f_and_n_cross
    );
    return system, n
end

#update equations
function update_interaction!(x,y,i,j,d2,output::FAndNcross,n)
    #F_ij
    r = y - x
    d = sqrt(d2)
    if d < r_cutoff
        f_ij = -k*(2*sigma-d)*r/d
    else
        f_ij =0
    end
    output.f[i] += f_ij
    output.f[j] -= f_ij
    #n_i cross n_j
    n_cross = cross(n[i],n[j])
    output.n_ij[i] += n_cross
    output.n_ij[j] -= n_cross
    return output
end

function calulate_changes(system, r, n)
    interactions = map_pairwise((x,y,i,j,d2,output) -> update_interaction!(x,y,i,j,d2,output,n), system)
    Fs = interactions.f
    ns = interactions.n_ij

    dn = zeros(SVector{3,Float64},N)
    dr = zeros(SVector{3,Float64},N)
    for i in eachindex(dn)
        dn[i] = (P_N(r[i],-J*ns[i])) * cross(normalize(r[i]),n[i])
        dr[i] = P_T(r[i], v_0*n[i]+mu*Fs[i])
    end
    if false
        print(dn[1])
        print(dr[1])
        print('\n')
    end

    return dr, dn
end


;


