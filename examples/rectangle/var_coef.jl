#
#Simple example of a Poisson problem with a variable coefficient:
#
#           y
#             |
#             |           u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |   -(au_x)_x - (au_y)_y = f     | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0               Lx
#
#Here, 
#
#      u = sin(ωx x) * cos(ωy y) 
#      a = exp(x-y)
#      f = ω exp(x-y) ( 2ω sin(ω x) cos(ω y) - cos ω(x-y) )
#      ωx = mπ / Lx
#      ωy = nπ / Ly
#

using FinElt
using FinElt.PlanarPoisson

include("params.jl")
const ωx = pi / Lx
const ωy = pi / Ly

function exact_u(x)
    return sin(ωx*x[1]) * cos(ωy*x[2])
end

function a(x)
    return exp(x[1]-x[2])
end

function f(x)
    return exp(x[1]-x[2]) * ( 
             ωx * ( ωx * sin(ωx*x[1]) - cos(ωx*x[1]) ) * cos(ωy*x[2])
           + ωy * ( ωy * cos(ωy*x[2]) - sin(ωy*x[2]) ) * sin(ωx*x[1])
           )
end

ell_f(v, z) = P1SourceTimesFunc!(v, z, f)
lin_functionals = [ ("Omega", ell_f) ]

maxerr = zeros(refinements+1)
@printf("%10s  %12s  %8s  %8s\n\n", 
        "N", "max error", "rate", "seconds")
for k = 0:refinements
    start = time()
    mesh = read_msh_file("rect$k.msh")    
    vp = VariationalProblem(mesh, ["Left", "Right"])
    add_bilin_form!(vp, "Omega", grad_dot_grad!, a)
    add_lin_functnl!(vp, "Omega", source_times_func!, f)
    # Implicit zero boundary conditions.
    A, b = assembled_linear_system(vp)
    ufree = A \ b
    uh = complete_soln(ufree, vp)
    u = get_nodal_vals(exact_u, mesh)
    finish = time()
    maxerr[k+1] = maximum(abs.(uh-u))
    N = length(ufree)
    if k == 0
        @printf("%10d  %12.4e\n", N, maxerr[k+1])
    else
        rate = log2(maxerr[k]/maxerr[k+1])
        elapsed = finish - start
        @printf("%10d  %12.4e  %8.4f  %8.4f\n", 
                N, maxerr[k+1], rate, elapsed)
    end
end

