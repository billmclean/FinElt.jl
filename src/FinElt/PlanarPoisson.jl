module PlanarPoisson

importall FinElt

export add_bilin_form!, add_lin_functnl!
export barycentric
export grad_dot_grad!, func_times_func!, source_times_func!
export bdry_func_times_func!, bdry_source_times_func!
export shape_params

function add_bilin_form!(vp::VariationalProblem, name::ASCIIString, 
                         elm_mat!::Function, coef=1.0)
    if !(name in keys(vp.mesh.elmtype))
        error("$name: unknown physical name")
    end
    if !(elm_mat! in keys(BILIN2GEOMTYPE))
        error("$name: unknown element matrix routine")
    end
    if BILIN2GEOMTYPE[elm_mat!] != vp.mesh.elmtype[name]
        error("Element type does not match bilinear form")
    end
    push!(vp.bilin_form, (name, elm_mat!,coef))
    return
end

function add_lin_functnl!(vp::VariationalProblem, name::ASCIIString, 
                          elm_vec!::Function, f::Function)
    if !(name in keys(vp.mesh.elmtype))
        error("$name: unknown physical name")
    end
    if !(elm_vec! in keys(LIN2GEOMTYPE))
        error("$name: unknown element matrix routine")
    end
    if LIN2GEOMTYPE[elm_vec!] != vp.mesh.elmtype[name]
        error("Element type does not match bilinear form")
    end
    push!(vp.lin_functnl, (name, elm_vec!, f))
    return
end

function add_bilin_form!(ep::EigenProblem, name::ASCIIString, 
                         elm_mat!::Function, side::ASCIIString, 
                         coef=1.0)
    if !(name in keys(ep.mesh.elmtype))
        error("$name: unknown physical name")
    end
    if !(elm_mat! in keys(BILIN2GEOMTYPE))
        error("$name: unknown element matrix routine")
    end
    if BILIN2GEOMTYPE[elm_mat!] != ep.mesh.elmtype[name]
        error("Element type does not match bilinear form")
    end
    if side == "LHS"
        push!(ep.LHS_bilin_form, (name, elm_mat!,coef))
    elseif side == "RHS"
        push!(ep.RHS_bilin_form, (name, elm_mat!,coef))
    else
        error("$side: illegal argument")
    end
    return
end

"""
b, centroid, area = barycentric(z)

For a planar triangle with vertices z_1, z_2, z_3, 
compute vectors b1, b2, b3 such that if

x = lambda_1 z_1 + lambda_2 z_2 + lambda_3 z_3

then the barycentric coordinates are given by

lambda_p = 1/3 + b_p dot ( x - centroid )

Also compute the centroid and area of the triangle.
"""
function barycentric(z::Matrix)
    b = zeros(2, 3)
    centroid = zeros(2)
    b[1,1] = z[1,1] - z[1,3]
    b[2,1] = z[2,1] - z[2,3]
    b[1,2] = z[1,2] - z[1,3]
    b[2,2] = z[2,2] - z[2,3]
    # overwrite B[1:2,1:2] with its inverse transpose
    area = b[1,1] * b[2,2] - b[1,2] * b[2,1]
    b[1,1], b[2,2] =  b[2,2]/area,  b[1,1]/area
    b[1,2], b[2,1] = -b[2,1]/area, -b[1,2]/area

    b[1,3] = -b[1,1] - b[1,2]
    b[2,3] = -b[2,1] - b[2,2]
    area /= 2
    centroid[1] = ( z[1,1] + z[1,2] + z[1,3] ) / 3
    centroid[2] = ( z[2,1] + z[2,2] + z[2,3] ) / 3
    return b, centroid, abs(area)
end

function grad_dot_grad!(A::Matrix,         # output 3x3
                        z::Matrix)         # input
    # A = element stiffness matrix
    b, centroid, area = barycentric(z)
    for q=1:3
        for p = q:3
            A[p,q] = area * ( b[1,p] * b[1,q] + b[2,p] * b[2,q] )
        end
    end
    for q = 2:3
        for p = 1:q-1
            A[p,q] = A[q,p]
        end
    end
    return
end                    

function func_times_func!(A::Matrix,       # output 3x3
                          z::Matrix)       # inputs 2x3
    # A = element mass matrix
    b, centroid, area = barycentric(z)
    fill!(A, area/12)
    for p=1:3
        A[p,p] *= 2
    end
    return
end

function bdry_func_times_func!(A::Matrix,   # output 2x2
                               z::Matrix)   # input  2x2
    # A = element mass matrix for the line segment from z[:,1]
    # to z[:,2]
    h = hypot(z[1,2]-z[1,1], z[2,2]-z[2,1])
    A[1,1] = h / 3
    A[2,1] = h / 6
    A[1,2] = A[2,1]
    A[2,2] = A[1,1]
    return
end

function source_times_func!(v::Vector,      # output 3
                            z::Matrix,      # input  2x3
                            f::Function)
    # v = element load vector, using piecewise-linear interpolation
    # of f
    A = zeros(3, 3)
    func_times_func!(A, z)
    fz = zeros(3)
    for p = 1:3
        fz[p] = f(z[:,p])
    end
    v[:] = A * fz    
    return
end

function bdry_source_times_func!(v::Vector,     # output 3
                                z::Matrix,      # input  2x2
                                g::Function)
    # v = element load vector, using piecewise-linear interpolation
    # of g
    A = zeros(2, 2)
    bdry_func_times_func!(A, z)
    gz = zeros(2)
    for p = 1:2
        gz[p] = g(z[:,p])
    end
    v[:] = A * gz    
    return
end

"""
Returns the diameter h and shape regularity measure h^2/area 
for the triangle with vertices z[:,1], z[:,2], z[:,3] 
(assumed to lie in the plane, i.e., only the first two 
components of each column are used).
"""
function shape_params(z::Matrix)
    b, centroid, area = barycentric(z)
    side = zeros(2,3)
    for p = 1:2
        side[p,1] = z[p,3] - z[p,2]
        side[p,2] = z[p,1] - z[p,3]
        side[p,3] = z[p,2] - z[p,1]
    end
    h = min(norm(side[:,1]), norm(side[:,2]), norm(side[:,3]))
    sr = h^2 / area
    return h, sr
end

function shape_params(mesh::Mesh)
    h  = Dict{ASCIIString, Vector{Float64}}()
    sr = Dict{ASCIIString, Vector{Float64}}()
    z = zeros(size(mesh.coord,1), 3)
    for name in keys(mesh.elms_of)
        if mesh.elmtype[name] != TRIANGLE
            continue
        end
        elms = mesh.elms_of[name]
        noelms = size(elms, 2)
        h[name]  = zeros(noelms)
        sr[name] = zeros(noelms)
        for k = 1:noelms
            for i = 1:size(z,1), p = 1:3
                z[i,p] = mesh.coord[i,elms[p,k]]
            end
            h[name][k], sr[name][k] = shape_params(z)
        end
    end
    return h, sr
end

const BILIN2GEOMTYPE = Dict( 
    grad_dot_grad!          => TRIANGLE,
    func_times_func!        => TRIANGLE,
    source_times_func!      => TRIANGLE,
    bdry_func_times_func!   => LINE,
    bdry_source_times_func! => LINE)

const LIN2GEOMTYPE = Dict( 
    source_times_func!      => TRIANGLE,
    bdry_source_times_func! => LINE)

end # module
