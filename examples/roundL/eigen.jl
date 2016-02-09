# An eigenproblem described in the LaTeX file eigen.tex.

using FinElt
using FinElt.PlanarPoisson

const h = 0.05
run(`gmsh -2 -clmax $h -o roundL.msh roundL.geo`)

mesh = read_msh_file("roundL.msh")
essential_bc = [ "Dirichlet" ]
ep = EigenProblem(mesh, essential_bc)
add_bilin_form!(ep, "Omega", grad_dot_grad!,   "LHS")
add_bilin_form!(ep, "Omega", func_times_func!, "RHS")

A, B = assembled_eigenproblem_matrices(ep)

d, vfree, nconv, niter, nmult, resid = eigs(A, B, which=:SM)

v = complete_soln(vfree, ep)

write_pos_file("roundL.pos") do fid
    for k = 1:size(v,2)
        save_warp_nodal_scalar_field(v[:,k], "v$k", fid)
    end
end

#run(`gmsh eigen.script`)
