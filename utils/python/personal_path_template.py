# THIS FILE SHOULD BE RENAMED <personal_path.py>

import os

# define OS-dependent variables
if (os.name == 'nt'):
    # windows
    print("running on windows")
    tab = '\\'
    cmd_clear = "cls"
elif (os.name == 'posix'):
    # unix
    tab = "/"
    cmd_clear = "clear"
else:
    print("unknown OS")
    quit()

# Define path to main folder
folder_name = '<NAME_OF_YOUR_MAIN_SW2D_FOLDER>' # Example: folder_name = 'sw2d'
sw2d_folder = os.getcwd().split(folder_name)[0] + folder_name + tab

# define executable path
sw2dConverter = sw2d_folder + '<LOCAL_PATH_TO_CONVERTER>'
sw2dPostProcess = sw2d_folder + '<LOCAL_PATH_TO_POSTPROCESS>'
sw2dPorosityMapper = sw2d_folder + '<LOCAL_PATH_TO_POROSITY_MAPPER>'
sw2dSolver = sw2d_folder + '<LOCAL_PATH_TO_SOLVER>'
sw2dUnitTests = sw2d_folder + '<LOCAL_PATH_TO_UNIT_TESTS>'
sw2dZone = sw2d_folder + '<LOCAL_PATH_TO_ZONE>'