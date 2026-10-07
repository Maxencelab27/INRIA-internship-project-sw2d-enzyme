import os
import re
import shutil
import sys

from personal_path import sw2dRotation

"""
This file is a tool to test the Rotation tool.
It runs sw2dRotation on the current directory 6 times with an angle of 60° and check if the final test is the same than the original one
Command line :
python testRotationTool.py path_to_simulation
"""

original_folder = sys.argv[1]
current_folder = sys.argv[1]

print("************ PROCEEDING 6 ROTATIONS OF 60 DEGREES ************")
# Run the rotation 6 times with an angle of 60°
for i in range(6) :
    command = sw2dRotation + " "  + current_folder + " 60"
    rotationStatus = os.system(command)
    if (rotationStatus != 0):
        print("[ERR.] Error occured in sw2dRotation. Abort.")
        print(command)
        exit()
    else:
        if (i!=6) :
            current_folder = os.path.abspath(current_folder) + "_rotated_60"

print("\n")
print("************ COMPARING RESULT AND ORIGINAL .2dm FILE ************")
success = True
same_values = True

# Open original 2dm
pf_list = os.listdir(os.path.join(original_folder, "primary_files"))
for file in pf_list :
    if (file[-4:] == ".2dm") :
        origin_2dm_file = open(os.path.join(original_folder, "primary_files", file))

# Open transformed 2dm
pf_list = os.listdir(os.path.join(current_folder, "primary_files"))
for file in pf_list :
    if (file[-4:] == ".2dm") :
        final_2dm_file = open(os.path.join(current_folder, "primary_files", file))

# Read and compare 2dm files
or_line = re.split("\\s+", origin_2dm_file.readline())
fi_line = re.split("\\s+", final_2dm_file.readline())
i_line = 1

max_diff_2dm = 10**(-2)

while ((or_line[0] != "") and (fi_line[0] != "")) :
    # If node line, compare values with an error margin
    if (or_line[0] == "ND") :
        if ((or_line[0]!=fi_line[0]) or (or_line[1]!=fi_line[1]) or (or_line[4]!=fi_line[4])) :
            success = False
            break
        x_or = float(or_line[2])
        y_or = float(or_line[3])
        x_fi = float(fi_line[2])
        y_fi = float(fi_line[3])
        if ((abs(x_or-x_fi)>max_diff_2dm) or (abs(y_or-y_fi)>max_diff_2dm)) :
            success = False
            same_values = False
            break

    # Every other line should be exactly the same
    elif (fi_line != or_line) :
        success = False
        break
            
    or_line = re.split("\\s+", origin_2dm_file.readline())
    fi_line = re.split("\\s+", final_2dm_file.readline())
    i_line += 1

# Print output message
if (not success) :
    if (not same_values) :
        print ("TEST FAILED IN 2dm : ORIGINAL AND FINAL VALUES ARE DIFFERENT")
        print ("First different values at line " + str(i_line))
        print ("Original values : " + str(x_or) + " " + str(y_or))
        print ("Final values : " + str(x_fi) + " " + str(y_fi))
        print("\n")
    else :
        print ("TEST FAILED IN 2dm : FOUND A DIFFERENT LINE")
        print ("First difference found at line " + str(i_line))
        print ("Original line : " + str(or_line))
        #print (str(elt) for elt in or_line)
        print ("Final line : " + str(fi_line))
        #print (str(elt) for elt in fi_line)
        print("\n")
else :
    print("TEST WAS SUCCESSFUL IN 2dm WITH ERROR MARGIN OF " + str(max_diff_2dm))
    print("\n")

# Close files
origin_2dm_file.close()
final_2dm_file.close()


if (os.path.isfile(os.path.join(original_folder, "input", "probe_locations.txt"))) :

    print("************ COMPARING RESULT AND ORIGINAL probe_locations FILE ************")
    success = True
    same_values = True

    # Open original probe_locations
    origin_pl_file = open(os.path.join(original_folder, "input", "probe_locations.txt"))

    # Open transformed probe_locations
    final_pl_file = open(os.path.join(current_folder, "input", "probe_locations.txt"))

    # Read and compare pl files
    or_line = re.split("\\s+", origin_pl_file.readline())
    fi_line = re.split("\\s+", final_pl_file.readline())
    i_line = 1

    max_diff_pl = 10**(-3)

    while ((or_line[0] != "") and (fi_line[0] != "")) :
        # "x	y	name" lines should be exactly the same
        if ((or_line[0]=="x") or (fi_line[0]=="x")) :
            if (fi_line != or_line) :
                success = False
                break

        # Compare values in other lines with error margin
        else :
            x_or = float(or_line[0])
            y_or = float(or_line[1])
            x_fi = float(fi_line[0])
            y_fi = float(fi_line[1])
            # Compare names
            if (or_line[2] != fi_line[2]) :
                success = False
                break
            # Compare values
            if ((abs(x_or-x_fi)>max_diff_pl) or (abs(y_or-y_fi)>max_diff_pl)) :
                success = False
                same_values = False
                break
                
        or_line = re.split("\\s+", origin_pl_file.readline())
        fi_line = re.split("\\s+", final_pl_file.readline())
        i_line += 1

    # Print output message
    if (not success) :
        if (not same_values) :
            print ("TEST FAILED IN probe_locations : ORIGINAL AND FINAL VALUES ARE DIFFERENT")
            print ("First different values at line " + str(i_line))
            print ("Original values : " + str(x_or) + " " + str(y_or))
            print ("Final values : " + str(x_fi) + " " + str(y_fi))
            print("\n")
        else :
            print ("TEST FAILED IN probe_locations : FOUND A DIFFERENT LINE")
            print ("First difference found at line " + str(i_line))
            print("\n")
    else :
        print("TEST WAS SUCCESSFUL IN probe_locations WITH ERROR MARGIN OF " + str(max_diff_pl))
        print("\n")

    # Close files
    origin_pl_file.close()
    final_pl_file.close()

# Delete all directories created during the process
current_folder = os.path.abspath(sys.argv[1]) + "_rotated_60"
for i in range(6) :
    shutil.rmtree(current_folder)
    current_folder = os.path.abspath(current_folder) + "_rotated_60"