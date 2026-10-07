import os
import sys

"""
This file reads nodestrings in a *.2dm file creates files qsection.txt and *.bc
A file primary_files/NS_conversion.txt is needed to specify bc code or qsection name for each nodestring
Command line :
python 2dmToBC.py <path_to_simulation_folder>
"""

# Check arguments
if (len(sys.argv) != 2) :
    print("sw2dNSPreProcess requires exactly 1 argument and received", len(sys.argv))
    print("The argument is the path to the simulation folder")
    exit(1)

# Check folders exist
sim_folder = sys.argv[1]
primary_folder = os.path.join(sim_folder, "primary_files")
input_folder = os.path.join(sim_folder, "input")

if not(os.path.exists(sim_folder)): # check simulation folder exists
    print("Simlation folder not found: " + sim_folder)
    print("Abort")
    exit(1)

if not(os.path.exists(primary_folder)): # check primary_files folder exists
    print("Primary files folder not found: " + primary_folder)
    print("Abort")
    exit(1)

if not(os.path.exists(input_folder)): # create input folder if it does not exist
    os.mkdir(input_folder)

# Open .2dm file 
primary_files = os.listdir(primary_folder)
mesh_filename = ""
for filename in primary_files :
    if filename[-4::] == ".2dm" :
        mesh_filename = os.path.join(primary_folder, filename)

if (mesh_filename == "") :
    print("2dm file not found in primary files folder : " + primary_folder)
    print("Abort")
    exit(1)

mesh_file = open(mesh_filename, mode="r")

# Read .2dm
nodestring = {}   # dict for nodestrings : keys = id, values = list of nodes id
edges = {}        # dict for edges : keys = tuple of nodes id, values = nb occurence of the edge

print("Reading mesh file :", mesh_filename)
for line in mesh_file :
    if ((line[0:3] == "E3T") or (line[0:3] == "E4Q")) : # Read edges
        line_split = line.split()
        nodes = [int(n) for n in line_split[2:-1]]
        for n in range(len(nodes)) :
            node_1 = nodes[n]
            node_2 = nodes[(n+1) % len(nodes)]
            # If edge already in edges dict, increment nb occurence
            if (node_1, node_2) in edges :
                edges[(node_1, node_2)] += 1
            elif (node_2, node_1) in edges :
                edges[(node_2, node_1)] += 1
            # else add it and nb occurence = 1
            else :
                edges[(node_1, node_2)] = 1
    if (line[0:2] == "NS") : # Read nodestrings
        line_split = line.split()
        id = int(line_split[-1])
        nodes = [abs(int(n)) for n in line_split[1:-1]]
        nodestring[id] = nodes
print("Done")

# Read conversion table and separate ns into qsections and bc
conversion_filename = os.path.join(primary_folder, "NS_conversion.txt")
conversion_file = open(conversion_filename, mode="r")

qsections = {}   # dict for qsections : keys = qsection name, values = list of nodes id
bc = {}          # dict for bc : keys = bc code, values = list of nodes id

print("Converting nodestrings to qsections and bc")
for line in conversion_file :
    if ((line[0] != "#") and (line[0] != "=") and (line.strip != "") and (line[0:2] != "NS")) :
        line_split = line.split()
        bc_code = line_split[1]
        qsection_name = line_split[2]
        if bc_code != "0" :
            bc[bc_code] = nodestring[int(line_split[0])]
        if qsection_name != "None" :
            qsections[qsection_name] = nodestring[int(line_split[0])]
print("  " + str(len(bc)) + " bc")
print("  " + str(len(qsections)) + " qsections")
print("Done")

# Write into .bc file
bc_filename = mesh_filename[:-4] + ".bc"
bc_file = open(bc_filename, mode="w")

print("Writing .bc file :", bc_filename)
bc_file.write(str(sum([len(nodes) - 1 for nodes in bc.values()])) + "\n") # write number of lines

for code, nodes in bc.items() :
    for i in range(len(nodes) - 1) :
        if (((nodes[i], nodes[i+1]) in edges) and (edges[(nodes[i], nodes[i+1])] > 1) or \
            ((nodes[i+1], nodes[i]) in edges) and (edges[(nodes[i+1], nodes[i])] > 1)) :
            print("[Warn] edge (" + str(nodes[i]) + ", " + str(nodes[i+1]) + ") from bc " + str(code) + " is internal")
            print("  It will be ignored")
        else :
            bc_file.write(str(nodes[i]) + "\t" + str(nodes[i+1]) + "\t" + code + "\n")
print("Done")

# Write into qsection.txt file
qsection_filename = os.path.join(input_folder, "qsection.txt")
qsection_file = open(qsection_filename, mode="w")

print("Writting qsections file :", qsection_filename)
qsection_file.write(str(sum([len(nodes)  - 1 for nodes in qsections.values()])) + "\n") # write number of lines

for name, nodes in qsections.items() :
    for i in range(len(nodes) - 1) :
        qsection_file.write(str(nodes[i]) + "\t" + str(nodes[i+1]) + "\t" + name + "\n")
print("Done")