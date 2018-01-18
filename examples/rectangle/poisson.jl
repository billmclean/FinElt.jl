#
#Simple example of a Poisson problem:
#
#           y
#             |
#             |           u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |    u - u_xx - u_yy = f         | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0               Lx
#
#Here, 
#
#      u = sin(omegax x) * cos(omegay y) 
#      f = ( 1 + omegax^2 + omegay^2 ) * sin(omegax x) * cos(omegay y)
#      omegax = nx pi / Lx
#      omegay = ny pi / Ly
#

using FinElt
using FinElt.PlanarPoisson

include("params.jl")

const omegax = 3 * pi / Lx
const omegay = pi / Ly

function exact_u(x)
    return sin(omegax*x[1]) * cos(omegay*x[2]) 
end

function f(x)
    c = 1 + omegax^2 + omegay^2
    return c * exact_u(x)
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
    add_bilin_form!(vp, "Omega", grad_dot_grad!)
    add_bilin_form!(vp, "Omega", func_times_func!)
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

