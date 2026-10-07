import struct
import sys

"""
This file reads binary file 'input/geometry.sw2d' and creates text file 'input/geometry_check.txt' to visualize the data
Command line :
python binaryReader.py <path_to_simulation_folder>
"""

# open binary file to read
ifile = open(sys.argv[1] + "/input/geometry.sw2d", mode="rb")
#ifile = open(sys.argv[1] + "/output/ResultHydro", mode="rb")

# open or create output txt file
ofile_name = sys.argv[1] + "/log/geometry.txt"
ofile = open(ofile_name, mode="w")

#-- write version number
version_lenght = struct.unpack("Q",ifile.read(8))[0]
version = ifile.read(version_lenght).decode("utf-8")
ofile.write(version + "\n")
#== write version number

#-- write number of elements
nb_cells = struct.unpack("Q",ifile.read(8))[0]
nb_edges = struct.unpack("Q",ifile.read(8))[0]
nb_nodes = struct.unpack("Q",ifile.read(8))[0]
nb_nodestrings = struct.unpack("Q",ifile.read(8))[0]

ofile.write(str(nb_cells) + "\t"
            + str(nb_edges) + "\t"
            + str(nb_nodes) + "\t"
            + str(nb_nodestrings) + "\n")
#== write number of elements

#-- write cells
nb_nodes_list = []
nb_edges_list = []
nb_neighs_list = []

for i in range(nb_cells) :
    id = struct.unpack("Q",ifile.read(8))[0]
    sms_id = struct.unpack("Q",ifile.read(8))[0]
    nodes_size = struct.unpack("Q",ifile.read(8))[0]
    edges_size = struct.unpack("Q",ifile.read(8))[0]
    neigh_size = struct.unpack("Q",ifile.read(8))[0]
    area = struct.unpack("d",ifile.read(8))[0]
    perimeter = struct.unpack("d",ifile.read(8))[0]
    x_g = struct.unpack("d",ifile.read(8))[0]
    y_g = struct.unpack("d",ifile.read(8))[0]
    zBottom = struct.unpack("d",ifile.read(8))[0]

    nb_nodes_list.append(nodes_size)
    nb_edges_list.append(edges_size)
    nb_neighs_list.append(neigh_size)

    ofile.write(str(id) + "\t"
                + str(sms_id) + "\t"
                + str(nodes_size) + "\t"
                + str(edges_size) + "\t"
                + str(neigh_size) + "\t"
                + "{:.15e}".format(area) + "\t"
                + "{:.15e}".format(perimeter) + "\t"
                + "{:.15e}".format(x_g) + "\t"
                + "{:.15e}".format(y_g) + "\t"
                + "{:.15e}".format(zBottom) + "\n")
#== write cells

#-- write nodes
nb_cells_list = []
nb_edges_per_node_list = []

for i in range(nb_nodes) :
    id = struct.unpack("Q",ifile.read(8))[0]
    sms_id = struct.unpack("Q",ifile.read(8))[0]
    x = struct.unpack("d",ifile.read(8))[0]
    y = struct.unpack("d",ifile.read(8))[0]
    z = struct.unpack("d",ifile.read(8))[0]
    cells_size = struct.unpack("Q",ifile.read(8))[0]
    edges_size = struct.unpack("Q",ifile.read(8))[0]

    nb_cells_list.append(cells_size)
    nb_edges_per_node_list.append(edges_size)
    
    ofile.write(str(id) + "\t"
                + str(sms_id) + "\t"
                + "{:.15e}".format(x) + "\t"
                + "{:.15e}".format(y) + "\t"
                + "{:.15e}".format(z) + "\t"
                + str(cells_size) + "\t"
                + str(edges_size) + "\n")
#== write nodes

#-- write edges
edges_boundary_size = 0

for i in range(nb_edges) :
    id = struct.unpack("Q",ifile.read(8))[0]
    ns = struct.unpack("Q",ifile.read(8))[0]
    type_lenght = struct.unpack("Q",ifile.read(8))[0]
    type = ifile.read(type_lenght).decode("utf-8")
    if (type=="boundary") :
            edges_boundary_size += 1
    cells_size = struct.unpack("Q",ifile.read(8))[0]
    node_1_id = struct.unpack("Q",ifile.read(8))[0]
    node_2_id = struct.unpack("Q",ifile.read(8))[0]
    x_m = struct.unpack("d",ifile.read(8))[0]
    y_m = struct.unpack("d",ifile.read(8))[0]
    z_m = struct.unpack("d",ifile.read(8))[0]
    lenght = struct.unpack("d",ifile.read(8))[0]
    d_norm = struct.unpack("d",ifile.read(8))[0]
    norm_x = struct.unpack("d",ifile.read(8))[0]
    norm_y = struct.unpack("d",ifile.read(8))[0]

    ofile.write(str(id) + "\t"
                + str(ns) + "\t"
                + type + "\t"
                + str(cells_size) + "\t"
                + str(node_1_id) + "\t"
                + str(node_2_id) + "\t"
                + "{:.15e}".format(x_m) + "\t"
                + "{:.15e}".format(y_m) + "\t"
                + "{:.15e}".format(z_m) + "\t"
                + "{:.15e}".format(lenght) + "\t"
                + "{:.15e}".format(d_norm) + "\t"
                + "{:.15e}".format(norm_x) + "\t"
                + "{:.15e}".format(norm_y) + "\n")
    
    first = True
    for j in range(cells_size) :
        if (first==False) :
            ofile.write("\t")
        first = False
        cell_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(cell_id))
    ofile.write("\n")

first = True
for i in range(edges_boundary_size) :
    if (first == False) :
        ofile.write("\t")
    first = False
    bc = struct.unpack("Q",ifile.read(8))[0]
    ofile.write(str(bc))
ofile.write("\n")
#== write edges

#-- convert nodestrings
for i in range(nb_nodestrings) :
    size = struct.unpack("Q",ifile.read(8))[0]
    ofile.write(str(size))
    for j in range(size) :
        edge_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(edge_id))
    ofile.write("\n")
#== convert nodestrings

#-- convert cell link with other elements
for i in range(nb_cells) :
    first = True
    for j in range(nb_nodes_list[i]) :
        if (first == False) :
            ofile.write("\t")
        first = False
        node_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(node_id))
    ofile.write("\n")

    first = True
    for j in range(nb_edges_list[i]) :
        if (first == False) :
            ofile.write("\t")
        first = False
        edge_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(edge_id))
    ofile.write("\n")

    first = True
    for j in range(nb_neighs_list[i]) :
        if (first == False) :
            ofile.write("\t")
        first = False
        neigh_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(neigh_id))
    ofile.write("\n")
#== convert cell link with other elements

#-- convert node link with other elements
for i in range(nb_nodes) :
    first = True
    for j in range(nb_cells_list[i]) :
        if (first == False) :
            ofile.write("\t")
        first = False
        cell_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(cell_id))
    ofile.write("\n")

    first = True
    for j in range(nb_edges_per_node_list[i]) :
        if (first == False) :
            ofile.write("\t")
        first = False
        edge_id = struct.unpack("Q",ifile.read(8))[0]
        ofile.write(str(edge_id))
    ofile.write("\n")
#== convert node link with other elements

#-- convert content of input.conv
nb_lines = struct.unpack("i",ifile.read(4))[0]

for i in range(nb_lines) :
    line_lenght = struct.unpack("Q",ifile.read(8))[0]
    line = ifile.read(line_lenght).decode("utf-8")
    ofile.write(line + "\n")
#== convert content of input.conv

# convert the position before input.conv
pos = struct.unpack("Q",ifile.read(8))[0]
ofile.write(str(pos))
ofile.write("\n")

#-- convert content from .2dm
nb_lines = struct.unpack("i",ifile.read(4))[0]

for i in range(nb_lines) :
    line_lenght = struct.unpack("Q",ifile.read(8))[0]
    line = ifile.read(line_lenght).decode("utf-8")
    ofile.write(line + "\n")
#== convert content from .2dm



print("File " + ofile_name + " created")
