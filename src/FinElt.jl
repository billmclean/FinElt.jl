module FinElt

export Mesh, read_msh_file, write_pos_file
export GeomType, LINE, TRIANGLE, TETRAHEDRON
export save_nodal_scalar_field, save_nodal_vector_field
export save_warp_nodal_scalar_field, write_format_version
export get_node_coords, get_nodal_vals

export VariationalProblem, assign_bdry_vals!
export DoF, degrees_of_freedom
export assembled_linear_system
export assembled_matrix, assembled_vector, complete_soln

include("FinElt/Gmsh.jl")
include("FinElt/FEM.jl")
include("FinElt/PlanarPoisson.jl")

end # module
