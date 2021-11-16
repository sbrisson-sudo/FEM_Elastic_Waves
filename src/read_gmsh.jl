import PyPlot ; const plt = PyPlot
import PyCall: pyimport
mpl =pyimport("matplotlib")
mpl.use("Qt5Agg")

"""
Read a gmsh file (version 2) of 2D quadrangles.
Return an array of nodes positions and an array of the nodes composing each element.
"""
function readGmsh(meshFile)

    nodes = Array{Float64}[]
    elements = Array{Int}[]

    open(meshFile, "r") do f

        while ! eof(f)

            if readline(f) == "\$Nodes"
                s = readline(f) # numbre of nodes

                s = readline(f)
                while s != "\$EndNodes"
                    l = map(x->parse(Float64, x), split(s))
                    push!(nodes, l[2:3])
                    s = readline(f)
                end

                readline(f) ; readline(f); readline(f) # comments

                s = readline(f)
                while s != "\$EndElements" 
                    l = map(x->parse(Int64, x), split(s))
                    l[2] == 3 && push!(elements, l[6:9])
                    s = readline(f)
                end
            end
        end
    end

    return nodes, elements

end


function plotMesh(nodes, elements)

    fig, ax = plt.subplots()

    ax.scatter([n[1] for n in nodes], [n[2] for n in nodes], zorder=5)

    for el in elements

        n1,n2,n3,n4 = el
        n1,n2,n3,n4 = nodes[n1],nodes[n2],nodes[n3],nodes[n4]
        
        plt.plot([n1[1],n2[1]], [n1[2],n2[2]], "k")
        plt.plot([n2[1],n3[1]], [n2[2],n3[2]], "k")
        plt.plot([n3[1],n4[1]], [n3[2],n4[2]], "k")
        plt.plot([n4[1],n1[1]], [n4[2],n1[2]], "k")
    end

    ax.set_aspect("equal", "box")
    plt.show()
end

if abspath(PROGRAM_FILE) == @__FILE__ 

    nodes, elements = readGmsh("../gmsh_meshes/octogon.msh")
    println(length(nodes), " nodes")
    println(length(elements), " elements")
    plotMesh(nodes, elements)

end
nodes, elements = readGmsh("./gmsh_meshes/octogon.msh")
println(length(nodes), " nodes")
println(length(elements), " elements")
plotMesh(nodes, elements)