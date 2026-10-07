import py2dm
import os
cwd = os.getcwd()

def mesh_1D(n_elements, length, width, x_min, y_min, z_all):
# @brief Script to build a 1D mesh made of quadrangle elements
# @par n_elements: number of elements in the mesh
# @par length: x-length of the domain
# @par width: y-width of the domain (and of every quadrangle)
# @par x_min: x-coordinate of the bottom left node
# @par y_min: y-coordinate of the bottom left node
# @par z_all: z-coordinate of all nodes (and elements)
    dx = length / n_elements
    with py2dm.Writer(cwd + '/mesh_1D.2dm') as mesh:
        mesh.node(1, 0.0, y_min, z_all)  # First bottom node
        mesh.node(2, 0.0, y_min+width, z_all)  # First top node
        n_nodes = 2
        x = x_min
        for elm in range(n_elements):
            n_nodes += 1
            x += dx
            mesh.node(n_nodes, x, y_min, z_all)  # Right bottom node to close element #elm
            n_nodes += 1
            mesh.node(n_nodes, x, y_min+width, z_all)  # Right top node to close element #elm
            mesh.element('E4Q', elm+1, n_nodes-3, n_nodes-2, n_nodes, n_nodes-1)  # Element #elm
    print("Mesh file 'mesh_1D.2dm' written in current directory")

# Examples:
# mesh_1D(10, 10, 1, 0, 0, 1.23)