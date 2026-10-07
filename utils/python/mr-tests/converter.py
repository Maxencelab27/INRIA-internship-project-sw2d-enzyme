import os
from personal_path import sw2dConverter, sw2dSolver

def runConverterTest(path, test):
    test_name = test.get("Folder")
    test_path = os.path.join(path, test_name)
    input_dir = os.path.join(test_path, "input")

    # Create folder for output file if needed
    if (os.path.exists(input_dir) == False):
        os.system("mkdir " + input_dir)

    # Removing previous output file
    if (os.path.exists(os.path.join(input_dir, "geometry.sw2d"))):
        os.remove(os.path.join(input_dir, "geometry.sw2d"))
    
    print("\n\n")
    # Running converter
    print("- Running converter for", test_name)
    command = sw2dConverter + " " + test_path
    conv_status = os.system(command)
    if not(conv_status == 0):
        return ("Err Converter")

    if (test_name != "no_input") :
        print("\n")
        # Running solver
        print("- Running solver for", test_name)    
        command = sw2dSolver + " " + test_path
        solver_status = os.system(command)
        if not(solver_status == 0):
            return ("Err Solver")
    
    return ("Success")