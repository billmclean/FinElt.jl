# Simple example of an eigenproblem:
#
#           y
#             |
#             |            u_y = 0
#          Ly ----------------------------------
#             |                                |
#             |                                |
#             |                                |
#       u = 0 |     u_xx + u_yy + λ u = 0      | u = 0
#             |                                |
#             |                                |
#             |                                |
#           -----------------------------------|-----  x
#                          u_y = 0              Lx
#
# Here, exact solutions are
#
#      u = sin( n π x / Lx ) cos( k π y / Ly )
#
#                      2               2
#      λ = ( n π / Lx )  + ( k π / Ly )
#
# for n = 1, 2, 3, ... and k = 0, 1, 2, ....

using FinElt
using FinElt.PlanarPoisson
using Printf
using Arpack

include("params.jl")
const nev = 4

λ = Float64[]
for n = 1:nev, k = 0:nev-1
    push!(λ, (n*pi/Lx)^2+(k*pi/Ly)^2)
end

sort!(λ)

essential_bc = [ "Left", "Right" ]

err = zeros(nev, refinements+1)
@printf("Eigenvalue errors\n\n")
@printf("%7s|", "N")
for j = 1:nev
    @printf("%16s|", "λ_$j    ")
end 
@printf(" elapsed\n")
@printf("%84s\n", "-"^84)
for k = 0:refinements
    start = time()
    mesh = read_msh_file("rect$k.msh")    
    ep = EigenProblem(mesh, essential_bc)
    add_bilin_form!(ep, "Omega", grad_dot_grad!,   :LHS)
    add_bilin_form!(ep, "Omega", func_times_func!, :RHS)
    A, B = assembled_eigenproblem_matrices(ep)
    d, nconv, niter, nmult, resid = eigs(A, B, nev=nev, which=:SM, 
                                         ritzvec=false)
    finish = time()
    err[:,k+1] = abs.(d-λ[1:nev])
    N = size(A, 1)
    if k == 0
        elapsed = finish - start
        @printf("%7d|", N) 
        for j = 1:nev
            @printf("%9.2e %6s|", err[j,k+1], "")
        end
        @printf(" %7.4f\n", elapsed)
    else
        elapsed = finish - start
        @printf("%7d|", N)
        for j = 1:nev
            rate = log2(err[j,k]/err[j,k+1])
            @printf("%9.2e %6.4f|", err[j,k+1], rate)
        end
        @printf(" %7.4f\n", elapsed)
    end
end

