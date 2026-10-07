import os
import shutil
import sys
import time

sw2d_folder = os.getcwd().split("test")[0]
sys.path.insert(0, os.path.join(sw2d_folder, "utils/python"))
import sw2dToolbox
from personal_path import sw2dConverter, sw2dSolver


def computeConverterSolverGraph(path, test_name):
    print("    Compute:", test_name)
    
    stop = False
    if not(os.path.exists(sw2dConverter)):
        print("sw2dConverter not found", sw2dConverter)
        stop = True
    if not(os.path.exists(sw2dSolver)):
        print("sw2dSolver not found", sw2dSolver)
        stop = True
    if (stop ==True):
        print("executables was not found. Abort.")
        quit() 

    converterStatus = -999
    sw2dStatus = -999
    PlotStatus = -999
    Info = ""
        
    # check if test exists
    if not(os.path.exists(os.path.join(path, test_name))):
        # test does not exist
        converterStatus = -1
        sw2dStatus = -1
        PlotStatus = -1        
        Info = "Test case not found"
    else:
        # test exists
        input_dir = os.path.join(path, test_name, "input")
        output_dir = os.path.join(path, test_name, "output")
        graphs_dir = os.path.join(path, test_name, "graph_scripts")

        # remove output dir if it exists
        if os.path.exists(output_dir):
            shutil.rmtree(output_dir)

        # clean previous png files if folder exists
        if (os.path.exists(graphs_dir)) :
            PngList = sw2dToolbox.getListPng(graphs_dir)
            for png in PngList:
                sw2dToolbox.deleteFile(os.path.join(graphs_dir, png))

        if not(os.path.exists(input_dir)):
            # no input folder
            converterStatus = -2
            sw2dStatus = -2
            Info = "Missing input folder"
        else:
            # input folder exists
            sw2dToolbox.deleteFile(os.path.join(input_dir, "geometry.sw2d")) # removing old geometry.sw2d file             
            
            # running converter
            Command = sw2dConverter + ' ' + os.path.join(path, test_name)
            converterStatus = os.system(Command)
            if not(converterStatus == 0):
                print("[ERR.] 2dm to geometry conversion failed. Abort.")
                print(Command)
                exit()
            else:
                print("[Info] Successfully converted 2dm file.")
            
                # running solver
                if not(os.path.exists(os.path.join(input_dir, "input.solv"))):
                    print(os.path.join(input_dir, "input.solv") + " missing")
                    sw2dStatus = -3
                    Info = "Input file not found"
                else:
                    Info = ""
                    Command = sw2dSolver + ' ' + os.path.join(path, test_name)
                    sw2dStatus = os.system(Command)	    #### equal to zero = SUCCESS
                    print("\n")
                    if not(sw2dStatus == 0):
                        Info = "Error in executing sw2dSolver"
    
                    else:
                        if not ((os.path.exists(graphs_dir) == True) and
                                (os.path.exists(os.path.join(graphs_dir, "graph_script.py")) == True)):
                            PlotStatus = -2
                            Info = "Graphic folder not found or missing python script: graph_script.py"
                        else:    
                            # run the graphic script
                            if (os.name == 'nt'):
                              # windows system
                                PlotStatus = os.system("python " + os.path.join(graphs_dir, "graph_script.py"))
                            else:
                                PlotStatus = os.system("python " + os.path.join(graphs_dir, "graph_script.py"))
                                    
                            if not (PlotStatus == 0):
                                Info = "Error during Python plot script"
                            else:
                                # Check at least one png has been generated
                                if len(sw2dToolbox.getListPng(graphs_dir)) > 0:
                                    PlotStatus = 0
                                else:
                                    PlotStatus = -3
                                    Info = "File(s) to plot not found"

    print("    ", sw2dToolbox.sep)

    return converterStatus, sw2dStatus, PlotStatus, Info
        
def runSolverTest(path, test_name):
    test_path = os.path.join(path, test_name)
    input_dir = os.path.join(test_path, "input")
    print("    ", sw2dToolbox.sep, "test:", test_name)
    
    CompTime = time.time()
    thisConverterStatus, thisSw2dStatus, thisPlotStatus, thisInfo = computeConverterSolverGraph(path, test_name)
    CompTime = time.time() - CompTime
    Sw2dName = os.path.join(path, test_name)[1 + len(os.path.join(sw2d_folder, "test/merge/solver")):]
    
    return Sw2dName, thisConverterStatus, thisSw2dStatus, thisPlotStatus, thisInfo, CompTime        