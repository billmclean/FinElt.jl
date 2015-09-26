const FILE_FORMAT = ("2.2", "0", "8")

immutable GeomType
    gmsh_code :: Int
    dimen     :: Int
    nonodes   :: Int
end

const LINE        = GeomType(1, 1, 2)
const TRIANGLE    = GeomType(2, 2, 3)
const TETRAHEDRON = GeomType(4, 3, 4)

const GETGEOMTYPE = Dict(1 => LINE,  
                         2 => TRIANGLE, 
                         4 => TETRAHEDRON)

# Mesh data type with fields
#
# coord[:,n]         = x, y, z coordinates of node n
# physname[m]        = string label of physical group with
#                      the numerical label 'm'
# physnum[s]         = numerical label of physical group
#                      the string label 's'
# elmtype[physname]  = element type for elements in physname
# elms_of[physname]  = array whose columns give the node numbers
#                      of the elements in physname
# nodes_of[physname] = nodes in physname
#

immutable Mesh
    coord    :: Array{Float64, 2}
    physdim  :: Dict{ASCIIString, Int}
    physnum  :: Dict{ASCIIString, Int}
    physname :: Dict{Int, ASCIIString}
    elmtype  :: Dict{ASCIIString, GeomType}
    elms_of  :: Dict{ASCIIString, Matrix{Int}}
    nodes_of :: Dict{ASCIIString, Vector{Int}}
end

function verify_next_line(s, f)
    line = readline(f)
    @assert strip(line) == s
end

function read_msh_file(fname::ASCIIString)
    # Read a .msh file created by Gmsh and return the
    # corresponding Mesh object.
    coord = Float64[]
    nonodes = 0
    line = ""
    elms = ASCIIString[]
    elmcount = Dict{Int,Int}()
    physdim  = Dict{ASCIIString,Int}()
    physnum  = Dict{ASCIIString,Int}()
    physname = Dict{Int,ASCIIString}()
    f = open(fname,"r") 
    while !eof(f)
        header = strip( readline(f) )
        if header == "\$MeshFormat"
            fmt = split( readline(f) )
            for j = 1:3
                @assert fmt[j] == FILE_FORMAT[j]
            end
            verify_next_line("\$EndMeshFormat", f)
        elseif header == "\$Nodes"
            line = readline(f)
            nonodes = parse(Int, line)
            for n = 1:nonodes
                line = readline(f)
                s = split(line)
                @assert n == parse(Int, s[1])
                for j = 1:3
                   push!(coord, parse(Float64, s[j+1]))
                end
            end
            verify_next_line("\$EndNodes", f)
        elseif header == "\$Elements"
            line = readline(f)
            noelms = parse(Int, line)
            for n = 1:noelms
                line = readline(f)
                push!(elms, line)
                row = split(line)
                physno = parse(Int, row[4])
                elmcount[physno] = get(elmcount, physno, 0) + 1
            end
            verify_next_line("\$EndElements", f)
        elseif header == "\$PhysicalNames"
            line = readline(f)
            nonames = parse(Int, line)
            for n = 1:nonames
                line = readline(f)
                s = split(line)
                dimen = parse(Int, s[1])
                num   = parse(Int, s[2])
                name  = strip(s[3], '"')
                physdim[name] = dimen
                physnum[name] = num
                physname[num] = name
            end
            verify_next_line("\$EndPhysicalNames", f)
        elseif header == "\$Comment"
            while true
                line = readline(f)
                if strip(line) == "\$EndComment"
                    break
                end
            end
        else
            msg = string("Unknown section header: ", line)
            error(msg)
        end
    end
    close(f)
    coord = reshape(coord, (3,nonodes))

    elmtype = Dict{ASCIIString, GeomType}()
    elms_of = Dict{ASCIIString,Matrix{Int}}()
    nodes_of = Dict{ASCIIString,Vector{Int}}()
    k = Dict{ASCIIString,Int}()
    for elm in elms
        s = split(elm)
        name = physname[parse(Int, s[4])]
        gmsh_code = parse(Int, s[2])
        if haskey(elmtype, name)
            @assert elmtype[name] == GETGEOMTYPE[gmsh_code]
        else
            elmtype[name] = GETGEOMTYPE[gmsh_code]
        end
    end
    for (name, element_type) in elmtype
        elms_of[name] = zeros(Float64, element_type.nonodes, 
                              elmcount[physnum[name]])
        k[name] = 0
    end
    for elm in elms
        line = [parse(Int, s) for s in split(elm)]
        # line uses the format
        # elm-number elm-type number-of-tags <tags> node-number-list
        # where the first tag is the numerical label of the physical
        # name.
        notags = line[3]
        name = physname[line[4]]
        k[name] += 1 
        elms_of[name][:,k[name]] = line[notags+4:end]
    end
    for (name, elm) in elms_of
        nodes_of[name] = noduplicates(elm[:])
    end
    return Mesh(coord, physdim, physnum, physname, elmtype, 
                elms_of, nodes_of)
end

function noduplicates(a::Vector{Int})
    # Return a vector containg the elements of a with
    # no duplicates, and sorted in increasing order.
    sort!(a)
    b = Array(Int, length(a))
    j = 1
    b[j] = a[1]
    for k = 2:length(a)
        if a[k] != a[k-1]
            j += 1
            b[j] = a[k]
        end
    end
    return b[1:j]
end

function save_nodal_scalar_field(u::Vector{Float64}, name::ASCIIString, 
                                 fid::IOStream, 
                                 time=0.0, timeidx=0)
    write(fid, "\$NodeData\n")
    write(fid, "1\n")             # number of string tags
    write(fid, "\"$name\"\n")     # label for scalar field u
    write(fid, "1\n")             # number of real tags
    s = @sprintf("%.15e\n", time)
    write(fid, s)
    write(fid, "3\n")             # number of integer tags
    s = @sprintf("%d\n", timeidx)
    write(fid, s)
    write(fid, "1\n")             # number of field components
    nonodes = length(u)
    s = @sprintf("%d\n", nonodes)
    write(fid, s)
    @assert length(u) == nonodes
    for nd = 1:nonodes
        s = @sprintf("%d  %.15e\n", nd, u[nd])
        write(fid, s)
    end
    write(fid, "\$EndNodeData\n")
end 

function save_warp_nodal_scalar_field(u::Array{Float64}, name::ASCIIString, 
                                      fid::IOStream, 
                                      time=0.0, timeidx=0)
    # Gmsh kludge to facilitate surface plotting of 2D scalar field via
    # the warp plugin.
    n = length(u)
    uwarp = zeros(3, length(u))
    uwarp[3,:] = u
    save_nodal_scalar_field(u, name, fid, time, timeidx)
    save_nodal_vector_field(uwarp, "$(name)_warp", fid, time, timeidx)   
end

function write_format_version(fid::IOStream)
    write(fid, "\$MeshFormat\n")
    write(fid, "$(FILE_FORMAT[1]) $(FILE_FORMAT[2]) $(FILE_FORMAT[3])\n")
    write(fid, "\$EndMeshFormat\n")
end

function save_nodal_vector_field(u::Array{Float64,2}, name::ASCIIString, 
                                 fid::IOStream, 
                                 time=0.0, timeidx=0)
    write(fid, "\$NodeData\n")
    write(fid, "1\n")             # number of string tags
    write(fid, "\"$name\"\n")     # label for scalar field u
    write(fid, "1\n")             # number of real tags
    s = @sprintf("%.15e\n", time)
    write(fid, s)
    write(fid, "3\n")             # number of integer tags
    s = @sprintf("%d\n", timeidx)
    write(fid, s)
    write(fid, "3\n")             # number of field components
    nonodes = size(u, 2)
    s = @sprintf("%d\n", nonodes)
    write(fid, s)
    @assert size(u) == (3, nonodes)
    for nd = 1:nonodes
        s = @sprintf("%d  %.15e  %.15e  %.15e\n", 
                     nd, u[1,nd], u[2,nd], u[3,nd])
        write(fid, s)
    end
    write(fid, "\$EndNodeData\n")
end 

function write_pos_file(f::Function, fname::ASCIIString)
    # Convenience function for use with do-syntax.
    fid = open(fname, "w")
    try
        write_format_version(fid)
        f(fid)
    finally
        close(fid)
    end
end

function write_pos_file(f::Function, posfname::ASCIIString, 
                        meshfname::ASCIIString)
    # Convenience function for use with do-syntax.
    # Include the mesh data.
    cp(meshfname, posfname)
    fid = open(posfname, "a")
    try
        f(fid)
    finally
        close(fid)
    end
end

function get_node_coords(mesh::Mesh)
   nonodes = size(mesh.coord, 2)
   x = zeros(nonodes)
   y = zeros(nonodes)
   z = zeros(nonodes)
   for nd = 1:nonodes
       x[nd] = mesh.coord[1,nd]
       y[nd] = mesh.coord[2,nd]
       z[nd] = mesh.coord[3,nd]
   end
   return x, y, z
end

function get_nodal_vals(f::Function, mesh::Mesh)
    # Returns the vector of nodal values of f
    nonodes = size(mesh.coord, 2)
    x = zeros(3)
    fvec = zeros(nonodes)
    for nd = 1:nonodes
        for j = 1:3
            x[j] = mesh.coord[j,nd]
        end
        fvec[nd] = f(x)
    end
    return fvec
end

function successive_refine(shape::ASCIIString, dimen::Integer, 
                           refinements::Integer, h0::Float64)
    # Read the file shape.geo and writes a sequence of files
    # shape0.msh, shape1.msh, ...,
    # by generating an initial mesh (shape0.msh) with maximum
    # element size h0, and performing the specified number of
    # refinements.  Thus, refinements+1 files are generated
    # altogether.
    infile = "$shape.geo"
    outfile = "$(shape)0.msh"
    run(`gmsh -optimize -$dimen -clmax $h0 -o $outfile $infile`)
    for k = 1:refinements
        infile = "$shape$(k-1).msh"
        outfile = "$shape$k.msh"
        run(`gmsh -v 0 -refine -o $outfile $infile`)
    end
    return
end
