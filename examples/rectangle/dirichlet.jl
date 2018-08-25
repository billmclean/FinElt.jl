#Simple example of a Dirichlet problem:
#
#           y
#             |
#             |     u = sin(omega x)
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |        u_xx + u_yy = 0         | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                          u = 0              Lx
#
#Here, omega = 2 n pi / Lx, and the exact solution is
#
#      u = sin(omega x) * sinh(omega y) / sinh(omega Ly)
#

using FinElt
using FinElt.PlanarPoisson
using Printf

include("params.jl")

const omega = 2 * pi / Lx

function exact_u(x)
    return sin(omega*x[1]) * sinh(omega*x[2]) / sinh(omega*Ly)
end

essential_bc = [ "Top", "Bottom", "Left", "Right" ]

maxerr = zeros(refinements+1)
@printf("%10s  %12s  %8s  %8s\n\n", 
        "N", "max error", "rate", "seconds")
for k = 0:refinements
    start = time()
    mesh = read_msh_file("rect$k.msh")    
    vp = VariationalProblem(mesh, essential_bc)
    assign_bdry_vals!(vp, "Top", exact_u)
    add_bilin_form!(vp, "Omega", grad_dot_grad!)
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

