import os
import sys

sw2d_folder = os.getcwd()
sw2d_folder = sw2d_folder.split("test")[0]
sys.path.insert(0, os.path.join(sw2d_folder, "utils/python"))
from personal_path import cmd_clear, sw2dConverter, sw2dPorosityMapper, sw2dSolver
    
def runPorosityMapperTest(porosity_root_path, test):
    test_name = test.get("Folder")
    test_path = os.path.join(porosity_root_path, test_name)
    input_dir = os.path.join(test_path, "input")
    ref_dir = os.path.join(test_path, "reference")

    # Removing previous output file
    if (os.path.exists(os.path.join(input_dir, "geometry.sw2d"))):
        os.remove(os.path.join(input_dir, "geometry.sw2d"))
    # TO DO : remove old porosity maps 

    conv = "Fail"
    pmap = "Fail"
    solv = "Fail"
    output = "Fail"
    
    print("\n\n")
    # Running converter
    print("- Running converter for", test_name)
    print("current folder", os.getcwd())
    command = sw2dConverter + " " + test_path
    conv_status = os.system(command)
    if (conv_status == 0):
        conv = "Success"
        print("\n")

        # Running porosity mapper
        print("- Running porosity mapper for", test_name)    
        command = sw2dPorosityMapper + " " + test_path
        mapper_status = os.system(command)
        if (mapper_status == 0):
            pmap = "Success"
            print("\n")

            # Running solver
            print("- Running solver for", test_name)    
            command = sw2dSolver + " " + test_path
            solver_status = os.system(command)
            if (solver_status == 0):
                solv = "Success"
            print("\n")

            # Checking outputs
            for filename in os.listdir(ref_dir) :
                if (filename[-17:] == "_porosity_map.txt") :
                    ref_file = os.path.join(ref_dir, filename)
                    result_file = os.path.join(input_dir, filename)
                    if (os.path.exists(result_file) and (open(ref_file, "r").read() == open(result_file, "r").read())) :
                        output = "Success"
    
    return (conv, pmap, solv, output)    