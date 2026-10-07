import os
import numpy as np
import sys
import read
import sw2d

sw2dFolder = os.getcwd()
sw2dFolder = sw2dFolder.split("utils")[0]
sys.path.insert(0, sw2dFolder + "utils/python")

def rotate(field_x, field_y, alpha):
    # Rotate a given double list field=(field_x, field_y) with a given angle alpha (in degrees)
    ## Check sizes
    n = len(field_x)
    rotated_field_x = field_x
    rotated_field_y = field_y
    if (len(field_y) != n):
        print("Error")
        exit_success = 0
    else:
        for n in range(n):
            alpha_rad = alpha / 180 * np.pi
            x1 = field_x[n]
            y1 = field_y[n]
            rotated_field_x[n] = np.cos(alpha_rad) * x1 - np.sin(alpha_rad) * y1
            rotated_field_y[n] = np.sin(alpha_rad) * x1 + np.cos(alpha_rad) * y1
        exit_success = 1
    return exit_success, rotated_field_x, rotated_field_y

def read_fields(date, fields):
    # Go to output folder
    csvFileName = "output/ResultHydro_p" + date + ".txt"
    # Read file
    result = sw2dToolbox.classHydroMap.hydroMap()
    read.mapFileHydro(result, csvFileName, table = fields)
    return result

def rotate_result(result, angle):
    result_rotated = result
    x = result.getVar('x')
    y = result.getVar('y')
    qx = result.getVar('qx')
    qy = result.getVar('qy')
    exit_success, rotated_x, rotated_y = rotate(x, y, angle)
    if (exit_success != 1):
        print("Error")
        exit
    exit_success, rotated_qx, rotated_qy = rotate(qx, qy, angle)
    if (exit_success != 1):
        print("Error")
        exit
    result_rotated.setVar('x',rotated_x)
    result_rotated.setVar('y',rotated_y)
    result_rotated.setVar('qx',rotated_qx)
    result_rotated.setVar('qy',rotated_qy)
    return result_rotated
