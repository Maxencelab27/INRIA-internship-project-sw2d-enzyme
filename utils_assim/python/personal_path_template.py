"""THIS FILE SHOULD BE RENAMED <personal_path.py>"""

import os

if (os.name == 'nt'):
    print("running on windows")
    tab = '\\'
    cmd_clear = "cls"
elif (os.name == 'posix'):
    tab = "/"
    cmd_clear = "clear"
else:
    print("unknown OS")
    quit()

folder_name = '<NAME_OF_YOUR_MAIN_SW2D_FOLDER>'
sw2d_folder = os.getcwd().split(folder_name)[0] + folder_name + tab

sw2dConverter = sw2d_folder + '<LOCAL_PATH_TO_CONVERTER>'
sw2dPostProcess = sw2d_folder + '<LOCAL_PATH_TO_POSTPROCESS>'
sw2dPorosityMapper = sw2d_folder + '<LOCAL_PATH_TO_POROSITY_MAPPER>'
sw2dSolver = sw2d_folder + '<LOCAL_PATH_TO_SOLVER>'
sw2dUnitTests = sw2d_folder + '<LOCAL_PATH_TO_UNIT_TESTS>'
sw2dZone = sw2d_folder + '<LOCAL_PATH_TO_ZONE>'
