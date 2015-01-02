# A mixed boundary-value problem.

using FinElt
using FinElt.PlanarPoisson

h = 0.05
run(`gmsh -2 -clmax $h -o roundL.msh roundL.geo`)

mesh = read_msh_file("roundL.msh")
essential_bc = [ "Dirichlet" ]
vp = VariationalProblem(mesh, essential_bc)
add_bilin_form!(vp,  "Omega",     grad_dot_grad!)
add_lin_functnl!(vp, "Omega",     source_times_func!,     x->2.0)
add_lin_functnl!(vp, "NeumannSW", bdry_source_times_func!, x->-0.75)

A, b = assembled_linear_system(vp)

ufree = A \ b

u = complete_soln(ufree, vp)
write_pos_file("roundL.pos") do fid
    save_warp_nodal_scalar_field(u, "u", mesh, fid)
end

run(`gmsh mixed.script`)
