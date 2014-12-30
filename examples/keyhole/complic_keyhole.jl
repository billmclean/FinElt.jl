using FinElt
using FinElt.PlanarPoisson

run(`gmsh -2 -o complic_keyhole.msh -clmax 0.05 complic_keyhole.geo`)
mesh = read_msh_file("complic_keyhole.msh")
essential_bc = [ "North", "South" ]
vp = VariationalProblem(mesh, essential_bc)
g(x) = -hypot(x[1],x[2])/2
assign_bdry_vals!(vp, "North", g)
add_bilin_form!(vp,  "Major", grad_dot_grad!, 1.0)
add_lin_functnl!(vp, "Major", source_times_func!, x->1.0)
add_bilin_form!(vp,  "Minor", grad_dot_grad!, 10.0)
add_lin_functnl!(vp, "Minor", source_times_func!, x->4.0)
add_lin_functnl!(vp, "East",  bdry_source_times_func!, x->-1.0)
add_lin_functnl!(vp, "West",  bdry_source_times_func!, x->-1.0)

A, b = assembled_linear_system(vp)

ufree = A \ b

u = complete_soln(ufree, vp)
write_pos_file("complic_keyhole.pos") do fid
    save_warp_nodal_scalar_field(u, "u", mesh, fid)
end    

run(`gmsh complic_keyhole.script`)
