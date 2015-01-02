# Simple example of a Neumann problem:
#
#           y
#             |
#             |     u_y = cos(omega x)
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#     u_x = 0 |    -u_xx - u_yy + u = 0        | u_x = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0              Lx
#
# Here, omega = 2 n pi / Lx, and the exact solution is
#
#      u = cos(omega x) * cosh(alpha y) / (alpha sinh(alpha Ly)
#
# where alpha^2 = 1 + omega^2.

using FinElt
using FinElt.PlanarPoisson

include("params.jl")

const omega = 2 * pi / Lx
const alpha = sqrt(1.0 + omega^2)

function exact_u(x)
    return cos(omega*x[1]) * cosh(alpha*x[2]) / (alpha*sinh(alpha*Ly))
end

function gN(x)
    return cos(omega*x[1])
end

maxerr = zeros(refinements+1)
@printf("%10s  %12s  %8s  %8s\n\n", 
        "N", "max error", "rate", "seconds")
for k = 0:refinements
    start = time()
    mesh = read_msh_file("rect$k.msh")    
    vp = VariationalProblem(mesh)
    add_bilin_form!(vp, "Omega", grad_dot_grad!)
    add_bilin_form!(vp, "Omega", func_times_func!)
    add_lin_functnl!(vp, "Top", bdry_source_times_func!, gN)
    A, b = assembled_linear_system(vp)
    ufree = A \ b
    uh = complete_soln(ufree, vp)
    u = get_nodal_vals(exact_u, mesh)
    finish = time()
    maxerr[k+1] = maximum(abs(uh-u))
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

