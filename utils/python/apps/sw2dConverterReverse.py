import os
import struct
import sys

"""
This file reads binary file 'input/geometry.sw2d' and creates files that have been used to create this geometry : input.conv, *.2dm, *.bc (optionnal)
Command line :
python sw2dConverterReverse.py <path_to_project_folder>
"""

# Open binary file to read
ifile = open(sys.argv[1] + "/input/geometry.sw2d", mode="rb")

# Go to the end of file
ifile.seek(-8, 2)

# Read the size of the input files content
content_size = struct.unpack("Q",ifile.read(8))[0]

# Go to the begining of the input files
ifile.seek(-8 -content_size, 2)

#-- create log folder if it does not exist
if not os.path.exists(sys.argv[1] + "/log") :
    os.mkdir(sys.argv[1] + "/log")
#== create log folder if it does not exist

#-- write input.conv
conv_filename = sys.argv[1] + "/log/input.conv"
print("Writting " + conv_filename)
conv_file = open(conv_filename, mode="w")

nb_lines = struct.unpack("i",ifile.read(4))[0]

mesh_filename = ""
bc_filename = ""

for i in range(nb_lines) :
    line_lenght = struct.unpack("Q",ifile.read(8))[0]
    line = ifile.read(line_lenght).decode("utf-8")

    line_split = line.split()
    if ((len(line_split) == 2) and (line_split[0] == "mesh_file")) :
        mesh_filename = sys.argv[1] + "/log/" + line_split[1]
    elif ((len(line_split) == 2) and (line_split[0] == "bc_file")) :
        bc_filename = sys.argv[1] + "/log/" + line_split[1]

    conv_file.write(line)
    if (i < nb_lines - 1) :
         conv_file.write("\n")

print("...done")
#== write input.conv

#-- convert content from .2dm
print("Writting " + mesh_filename)
mesh_file = open(mesh_filename, mode="w")

nb_lines = struct.unpack("i",ifile.read(4))[0]

for i in range(nb_lines) :
    line_lenght = struct.unpack("Q",ifile.read(8))[0]
    line = ifile.read(line_lenght).decode("utf-8")
    mesh_file.write(line)
    if (i < nb_lines - 1) :
         mesh_file.write("\n")

print("...done")
#== convert content from .2dm

#-- convert content from .bc
if (bc_filename != "") :
    print("Writting " + bc_filename)
    bc_file = open(bc_filename, mode="w")

    nb_lines = struct.unpack("i",ifile.read(4))[0]

    for i in range(nb_lines) :
        line_lenght = struct.unpack("Q",ifile.read(8))[0]
        line = ifile.read(line_lenght).decode("utf-8")
        bc_file.write(line)
        if (i < nb_lines - 1) :
            bc_file.write("\n")

    print("...done")
#== convert content from .bc
