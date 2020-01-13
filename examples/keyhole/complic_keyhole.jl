using FinElt
using FinElt.PlanarPoisson

run(`gmsh -2 -o complic_keyhole.msh -clmax 0.05 -format msh22
complic_keyhole.geo`)
mesh = read_msh_file("complic_keyhole.msh")
essential_bc = [ "DarkBlue", "Black" ]
vp = VariationalProblem(mesh, essential_bc)
g(x) = -hypot(x[1],x[2])/2
assign_bdry_vals!(vp, "DarkBlue", g)
add_bilin_form!(vp,  "LightBlue", grad_dot_grad!, 1.0)
add_bilin_form!(vp,  "LightGreen", grad_dot_grad!, 10.0)
add_lin_functnl!(vp, "LightBlue", source_times_func!, x->1.0)
add_lin_functnl!(vp, "LightGreen", source_times_func!, x->4.0)
add_lin_functnl!(vp, "Violet",  bdry_source_times_func!, x->-1.0)

A, b = assembled_linear_system(vp)

ufree = A \ b

u = complete_soln(ufree, vp)
write_pos_file("complic_keyhole.pos") do fid
    save_warp_nodal_scalar_field(u, "u", fid)
end    

run(`gmsh complic_keyhole.script`)
