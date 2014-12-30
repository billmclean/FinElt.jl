using FinElt
using FinElt.PlanarPoisson

mesh = read_msh_file("keyhole.msh")
#h, sr = shape_params(mesh)
essential_bc = [ "Gamma" ]
f(x) = 4.0
vp = VariationalProblem(mesh, essential_bc)
add_bilin_form!(vp,  "Omega", grad_dot_grad!)
add_lin_functnl!(vp, "Omega", source_times_func!, f)

A, b = assembled_linear_system(vp)

ufree = A \ b

u = complete_soln(ufree, vp)
write_pos_file("keyhole.pos") do fid
    save_warp_nodal_scalar_field(u, "u", mesh, fid)
end
