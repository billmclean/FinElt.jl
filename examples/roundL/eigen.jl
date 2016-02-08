# An eigenproblem described in the LaTeX file eigen.tex.

using FinElt
using FinElt.PlanarPoisson

h = 0.05
run(`gmsh -2 -clmax $h -o roundL.msh roundL.geo`)

mesh = read_msh_file("roundL.msh")
essential_bc = [ "Dirichlet" ]
ep = EigenProblem(mesh, essential_bc)
add_bilin_form!(ep,  "Omega", grad_dot_grad!,   "LHS")
add_bilin_form!(ep,  "Omega", func_times_func!, "RHS")

A, B = assembled_eigenproblem_matrices(ep)

#u = complete_soln(ufree, vp)
#write_pos_file("roundL.pos") do fid
#    save_warp_nodal_scalar_field(u, "u", fid)
#end

#run(`gmsh eign.script`)
