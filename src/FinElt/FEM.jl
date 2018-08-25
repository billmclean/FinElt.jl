"""
Data structure describing the FEM degrees of freedom.

isfree[nd]     is true if nd is free, false is nd is fixed.
freenode[k]    is the gmsh node number of the kth free node.
fixednode[k]   is the gmsh node number of the kth fixed node.
node2free[nd]  = k if nd = freenode[k], -1 otherwise.
node2fixed[nd] = k if nd = fixednode[k], -1 otherwise.
"""
struct DoF
    isfree     :: Vector{Bool}    
    freenode   :: Vector{Integer}
    fixednode  :: Vector{Integer}
    node2free  :: Vector{Integer}
    node2fixed :: Vector{Integer}
end

struct VariationalProblem
    mesh         :: Mesh
    dof          :: DoF
    essential_bc :: Vector{String}
    bilin_form   :: Vector{Any}
    lin_functnl  :: Vector{Any}
    ufixed       :: Vector{Float64}
end

struct EigenProblem
    mesh           :: Mesh
    dof            :: DoF
    essential_bc   :: Vector{String}
    LHS_bilin_form :: Vector{Any}
    RHS_bilin_form :: Vector{Any}
end

function VariationalProblem(mesh::Mesh, 
                            essential_bc::Array{String})
    for name in essential_bc
        if !(name in keys(mesh.nodes_of))
            error("$name: unknown physical name")
        end
    end
    dof = degrees_of_freedom(mesh, essential_bc)
    bilin_form   = Any[]
    lin_functnl  = Any[]
    nofixed = length(dof.fixednode)
    # Zero boundary condition by default.
    ufixed = zeros(nofixed)
    return VariationalProblem(mesh, dof, essential_bc, 
                              bilin_form, lin_functnl, ufixed)
end

# Use in case of pure Neumann bc.
VariationalProblem(mesh::Mesh) = VariationalProblem(mesh, String[])

function EigenProblem(mesh::Mesh, essential_bc::Array{String})
    for name in essential_bc
        if !(name in keys(mesh.nodes_of))
            error("$name: unknown physical name")
        end
    end
    dof = degrees_of_freedom(mesh, essential_bc)
    LHS_bilin_form   = Any[]
    RHS_bilin_form   = Any[]
    nofixed = length(dof.fixednode)
    return EigenProblem(mesh, dof, essential_bc,  
                        LHS_bilin_form, RHS_bilin_form)
end

function assign_bdry_vals!(vp::VariationalProblem, 
                           name::String, g::Float64)
    if !(name in vp.essential_bc)
        error("$name: not listed in essential_bc")
    end
    for nd in vp.mesh.nodes_of[name]
        i = vp.dof.node2fixed[nd]
        vp.ufixed[i] = g
    end
end

function assign_bdry_vals!(vp::VariationalProblem, 
                           name::String, g::Function)
    if !(name in vp.essential_bc)
        error("$name: not listed in essential_bc")
    end
    for nd in vp.mesh.nodes_of[name]
        i = vp.dof.node2fixed[nd]
        x = vp.mesh.coord[:,nd]
        vp.ufixed[i] = g(x)
    end
end

"""
Returns the DoF object for the given mesh and essential boundary
conditions.
"""
function degrees_of_freedom(mesh::Mesh, essential_bc::Array{String})
    nonodes = size(mesh.coord, 2)
    isfree = Array{Bool}(undef, nonodes)
    fill!(isfree, true)
    for name in essential_bc
        if !(name in values(mesh.physname))
            error("Unknown physical region $name")
        end
        elm = mesh.elms_of[name]
        for k = 1:size(elm, 2)
            for p = 1:size(elm, 1)
                isfree[elm[p,k]] = false
            end                
        end
    end
    nofree = 0
    nofixed = 0
    for nd = 1:nonodes
        if isfree[nd]
            nofree += 1
        else
            nofixed += 1
        end
    end
    kfree = 0
    kfixed = 0
    freenode  = zeros(Integer, nofree)
    fixednode = zeros(Integer, nofixed)
    for nd = 1:nonodes
        if isfree[nd]
            kfree += 1
            freenode[kfree] = nd
        else
            kfixed += 1
            fixednode[kfixed] = nd
        end
    end
    node2free  = zeros(Integer, nonodes)
    fill!(node2free, -1)
    for k = 1:nofree
        nd = freenode[k]   
        node2free[nd] = k
    end
    node2fixed = zeros(Integer, nonodes)
    fill!(node2fixed, -1)
    for k = 1:nofixed
        nd = fixednode[k]
        node2fixed[nd] = k
    end
    return DoF(isfree, freenode, fixednode, node2free, node2fixed)
end

function assembled_linear_system(vp::VariationalProblem)
    mesh = vp.mesh
    dof = vp.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    A = assembled_matrix(vp)
    b = zeros(nofree)
    assembled_vector!(b, vp)
    if nofixed > 0
        b .= b - A[:,nofree+1:end] * vp.ufixed
    end
    return A[:,1:nofree], b
end

function assembled_eigenproblem_matrices(ep::EigenProblem)
    mesh = ep.mesh
    dof = ep.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    A = sparse(Int64[], Int64[], Float64[], nofree, nofree)
    B = sparse(Int64[], Int64[], Float64[], nofree, nofree)
    for (name, elm_mat!, coef) in ep.LHS_bilin_form
        next = assembled_matrix(name, elm_mat!, coef, mesh, dof)
        A += next[:,1:nofree]
    end
    for (name, elm_mat!, coef) in ep.RHS_bilin_form
        next = assembled_matrix(name, elm_mat!, coef, mesh, dof)
        B += next[:,1:nofree]
    end
    return A, B
end

function element_matrix(elm_A::Matrix{Float64}, z::Matrix{Float64}, 
                        coef::Union{Float64,Function}, k::Int64, 
                        elm_mat!::Function)
    elm_mat!(elm_A, z, coef)
end

function element_matrix(elm_A::Matrix{Float64}, z::Matrix{Float64}, 
                        coef::Vector{Float64}, k::Int64, elm_mat!::Function)
    elm_mat!(elm_A, z, coef[k])
end

function assembled_matrix(vp::VariationalProblem)
    mesh = vp.mesh
    dof = vp.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    A = sparse(Int64[], Int64[], Float64[], nofree, nofree+nofixed)
    for (name, elm_mat!, coef) in vp.bilin_form
        next = assembled_matrix(name, elm_mat!, coef, mesh, dof)
        A += next
    end
    return A
end

function assembled_matrix(name::String, elm_mat!::Function, 
                          coef::Union{Float64, Function, Vector{Float64}},
                          mesh::Mesh, dof::DoF)
    # name = physical name of sub-domain
    # coef    = constant coefficient
    # coef(x) = value of coef at x
    # coef[j] = value of coef at centroid of jth element
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    isfree = dof.isfree
    I = Int64[]
    J = Int64[]
    V = Float64[]
    elm = mesh.elms_of[name] 
    nonodes_per_elm = size(elm, 1)
    elm_A = zeros(nonodes_per_elm, nonodes_per_elm)
    z = zeros(size(mesh.coord,1), nonodes_per_elm)
    for k = 1:size(elm, 2)
        for i = 1:size(z,1), p = 1:nonodes_per_elm
            z[i,p] = mesh.coord[i,elm[p,k]]
        end
        element_matrix(elm_A, z, coef, k, elm_mat!)
        for p = 1:nonodes_per_elm
            if !isfree[elm[p,k]]
                continue
            end
            for q = 1:nonodes_per_elm
                push!(I, dof.node2free[elm[p,k]])
                if isfree[elm[q,k]]
                    push!(J, dof.node2free[elm[q,k]])
                else
                    push!(J, nofree + dof.node2fixed[elm[q,k]])
                end
                push!(V, elm_A[p,q])
            end
        end
    end
    A = sparse(I, J, V, nofree, nofree+nofixed)
    return A
end

function assembled_vector!(loadvec::Vector{Float64}, vp::VariationalProblem)
    mesh = vp.mesh
    dof = vp.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    @assert length(loadvec) == nofree
    fill!(loadvec, 0.0)
    for (name, elm_vec!, f) in vp.lin_functnl
        assembled_vector!(loadvec, name, elm_vec!, f, mesh, dof)
    end
end

function assembled_vector(name::String, elm_vec!::Function, 
                          f::Function, mesh::Mesh, dof::DoF)
    nofree = length(dof.freenode)
    loadvec = zeros(nofree)
    assembled_vector!(loadvec, name, elm_vec!, f, mesh, dof)
    return loadvec
end

function assembled_vector!(loadvec::Vector{Float64}, 
                           name::String, elm_vec!::Function, 
                           f::Function, mesh::Mesh, dof::DoF)
    nofree = length(dof.freenode)
    @assert length(loadvec) == nofree
    nofixed = length(dof.fixednode)
    isfree = dof.isfree
    elm = mesh.elms_of[name]
    nonodes_per_elm = size(elm, 1)
    elm_v = zeros(nonodes_per_elm)
    z = zeros(size(mesh.coord,1), nonodes_per_elm)
    for k = 1:size(elm, 2)
        for p = 1:nonodes_per_elm
            for i = 1:size(z,1)
                z[i,p] = mesh.coord[i,elm[p,k]]
            end
        end
        elm_vec!(elm_v, z, f)
        for p = 1:nonodes_per_elm
            if isfree[elm[p,k]]
                i = dof.node2free[elm[p,k]]
                loadvec[i] += elm_v[p]
            end
        end
    end
end

"""
Returns the vector of nodal values, including both the free
and fixed nodes.
"""
function complete_soln(ufree::Vector{Float64}, vp::VariationalProblem)
    dof = vp.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    u = zeros(nofree+nofixed)
    for k = 1:nofree
        nd = dof.freenode[k]
        u[nd] = ufree[k]
    end    
    for k = 1:nofixed
        nd = dof.fixednode[k]
        u[nd] = vp.ufixed[k]
    end
    return u
end

function complete_soln!(u::Vector{Float64}, ufree::Vector{Float64},
                        ufixed::Vector{Float64}, dof::DoF)
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    @assert length(ufree) == nofree
    @assert length(ufixed) == nofixed
    for k = 1:nofree
        nd = dof.freenode[k]
        u[nd] = ufree[k]
    end    
    for k = 1:nofixed
        nd = dof.fixednode[k]
        u[nd] = ufixed[k]
    end
end

function complete_soln(vfree::Matrix{Float64}, ep::EigenProblem)
    dof = ep.dof
    nofree = length(dof.freenode)
    nofixed = length(dof.fixednode)
    nev = size(vfree, 2)
    v = zeros(nofree+nofixed, nev)
    for j = 1:nev, k = 1:nofree
        nd = dof.freenode[k]
        v[nd,j] = vfree[k,j]
    end    
    return v
end
