# Simple example of an inhomogeneous heat equation:
#
#           y
#             |
#             |     u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#     u_x = 0 |    u_t - a(u_xx + u_yy) = f    | u_x = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                         u_y = 0              Lx
#
# Here, a>0, the initial data is u0 = 0 and the source term is
#
#      f(x, y, t) = t cos(pi x / Lx) cos(pi y / Ly)
#
# The exact solution is
#
#                    lambda t - 1 - exp(-lambda t) 
#      u(x, y, t) =  ----------------------------- phi(x,y)
#                          lambda^2
#
# where
#
#      phi(x,y) = cos(pi x / Lx ) cos(pi y / Ly)
#      lambda = a pi^2 ( 1 / Lx^2 + 1 / Ly^2 )
#
using FinElt
using FinElt.PlanarPoisson
include("ode23s.jl")

const a = 0.1
const Lx = 1.0
const Ly = 1.0
const lambda = a * pi^2 * ( 1/Lx^2 + 1/Ly^2 )
const T = 2.0

function phi(x)
    return cos(pi*x[1]/Lx)* cos(pi*x[2]/Ly)
end

function exact_u(x, t)
    v = ( lambda*t - 1 + exp(-lambda*t) ) / lambda^2
    return v * phi(x)
end

function f(x, t)
    return t * phi(x)
end

function RHS(t, u, S, mesh, dof)
    F = assembled_vector("Omega", source_times_func!, 
                         x -> f(x,t), mesh, dof)
    return F - S * u
end

maxnorm(x) = norm(x, Inf)

start = time()
mesh = read_msh_file("../rectangle/rect2.msh")    
dof = degrees_of_freedom(mesh, String[])
M = assembled_matrix("Omega", func_times_func!, 1.0, mesh, dof)
S = assembled_matrix("Omega", grad_dot_grad!, a, mesh, dof)
u0 = zeros(length(dof.freenode))
t, u = ode23s((t,u) -> RHS(t, u, S, mesh, dof), u0, [0,T]; 
             reltol=1.0e-4, abstol=1.0e-4, mass=M, norm=maxnorm)

N = length(t)
err = zeros(N)
for n = 1:N
    un = get_nodal_vals(x->exact_u(x,t[n]), mesh)
    err[n] = maxnorm(un-u[n])
end

finish = time()
elapsed = finish - start

figure(1)
plot(t, err)
xlabel(L"$t$")
title("Error solving the inhomogeneous heat equation")
grid(true)
