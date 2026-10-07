##from datetime import datetime
#import math
##import csv
import os
import shutil
import sys
import time

sw2d_folder = os.getcwd()
sw2d_folder = sw2d_folder.split("test")[0]
from personal_path import cmd_clear, sw2dConverter, sw2dZone
sys.path.insert(1, os.path.join(sw2d_folder,"utils/python"))
import sw2dToolbox

def processFolder(Dossier, recursive = True):
    # define the storing table
    ZoneNames = []
    ZoneStatus = []
    ConverterStatus = []
    OutputStatus = []
    PlotStatus = []
    Info = []
    CompTime = []
    total = 0

    FolderList = sw2dToolbox.getFolderList(Dossier)
    
    for Test in FolderList:
        if os.path.exists(os.path.join(Dossier, Test, "primary_files")):
            # primary_files subfolder exists => running sw2dZone
            print("    ", sw2dToolbox.sep, "test:", Test)
            ListZoneNames = []
            ListConverterStatus = []
            ListZoneStatus = []
            ListOutputStatus = []
            ListInfo = []
            ListCompTime = []
            t_0 = time.time()
            a, b, c, d = Compute(Dossier, Test)
            t = time.time() - t_0
            ListZoneNames.append(os.path.join(Dossier, Test)[1 + len(os.path.join(sw2d_folder, "test/merge/zone")):])
            ListConverterStatus.append(a)
            ListZoneStatus.append(b)
            ListOutputStatus.append(c)
            ListInfo.append(d)
            ListCompTime.append(t)

        else:
            # running Python
            ListZoneNames = []
            ListConverterStatus = []
            ListZoneStatus = []
            ListInfo = []
            ListCompTime = []
            src = os.path.join(sw2d_folder, "test/merge/zone", "testZone.py")
            dst = os.path.join(Dossier, Test, "testZone.py")
            shutil.copyfile(src, dst)
            print("======================")
            print("  subfolder:", os.path.join(Dossier, Test))
            ListZoneNames, ListConverterStatus, ListZoneStatus, ListOutputStatus, ListInfo, ListCompTime = processFolder(os.path.join(Dossier, Test))

        for name in ListZoneNames:
            ZoneNames.append(name)
            
        for name in ListZoneStatus:
            ZoneStatus.append(name)
            
        for name in ListConverterStatus:
            ConverterStatus.append(name)
            
        for name in ListOutputStatus:
            OutputStatus.append(name)

        for name in ListInfo:
            Info.append(name)
            
        for name in ListCompTime:
            CompTime.append(name)

    return ZoneNames, ConverterStatus, ZoneStatus, OutputStatus, Info, CompTime
    
def checkStatus(chemin):
    Status = 0
    RefFolder = os.path.join(chemin, "reference")
    InputFolder = os.path.join(chemin, "input")
    refFiles = sw2dToolbox.getRefFilesList(RefFolder)
    
    for fich in refFiles:
        refFile = os.path.join(RefFolder, fich)
        inputFile = os.path.join(InputFolder, fich)
        if not(os.path.isfile(os.path.join(InputFolder, fich))):
            Status -= 1
        else:
            refContent = open(refFile, "r").read()
            inputContent = open(inputFile, "r").read()
            if not(refContent == inputContent):
                Status -= 1
            
    return Status

def Compute(Root, Test):
    Folder = os.path.join(Root, Test)
    print("    Running sw2dZone for", Folder)

    stop = False
    if not(os.path.exists(sw2dConverter)):
        print("sw2dConverter not found", sw2dConverter)
        stop = True
    if not(os.path.exists(sw2dZone)):
        print("sw2dZone not found", sw2dZone)
        stop = True
    if (stop ==True):
        print("executables was not found. Abort.")
        quit()    
        
    converterStatus = -3
    ZoneStatus = -3
    OutputStatus = -3
    Info = ""

    # check if test exists
    if not(os.path.exists(Folder)):
        converterStatus = -2
        ZoneStatus = -2
        OutputStatus = -2
        Info = "Test case not found"
    else:
        primary_dir = os.path.join(Folder, "primary_files")
        input_dir = os.path.join(Folder, "input")

        if not(os.path.exists(primary_dir)):
            converterStatus = -1
            ZoneStatus = -1
            Info = "Missing input files"
        else:
            # Running converter
            ## Removing previous output files
            if (os.path.exists(os.path.join(input_dir, "geometry.sw2d")) == True):
                os.remove(os.path.join(input_dir, "geometry.sw2d"))

            Command = sw2dConverter + " " + Folder
            converterStatus = os.system(Command)
            if not(converterStatus == 0):
                Info = "sw2dConverter failed"
            else:
            # Running zone
                if not(os.path.isfile(sw2dZone)):
                    print("    sw2dZone was not found. Abort.")
                    quit()
                else:
                    print("    Start with running sw2dZone to generate parameters files in", Test)
                    Command = sw2dZone + " " + Folder
                    ZoneStatus = os.system(Command)	    #### equal to zero = SUCCESS    
                    
                    if (ZoneStatus == 0):
                        OutputStatus = checkStatus(Folder)

    print("    ", sw2dToolbox.sep)

    return converterStatus, ZoneStatus, OutputStatus, Info

def writeZoneConsole(ZoneNames, ConverterStatus, ZoneStatus, OutputStatus, Info, CompTime):
    ZoneNamesSorted = sorted(ZoneNames)
    ZoneStatusSorted = []
    ConverterStatusSorted = []
    OutputStatusSorted = []
    InfoSorted = []
    CompTimeSorted = []

    # print the summary
    MaxSizeZone = max([len(Test) for Test in ZoneNames])
    if not ((MaxSizeZone % 2) == 0): # non divisible par 2
        MaxSizeZone = MaxSizeZone + 1

    print("\n")
    print("***** TESTS SUMMARY *****")

    if (sum(ZoneStatus) == 0):
        print("- sw2d tests processed successfully")
    else:
        print("- sw2d tests encountered errors")

    if (sum(OutputStatus) == 0):
        print("- output files are correct")
    else:
        print("- error in the output files")

    print("- TOTAL TIME: " + str(sum(CompTime)))
    print(" ")

    print(sw2dToolbox.repeatString(" ", (MaxSizeZone - 4) / 2) + "TEST" +
          sw2dToolbox.repeatString(" ", (MaxSizeZone - 4) / 2) +
          "\tCONV.\tZONE\tOUPUT\tTIME\tINFO")

    nb_errors = 0

    for idx, Test in enumerate(ZoneNamesSorted):
        if (ConverterStatus[ZoneNamesSorted.index(Test)] == 0):
            ConverterStatusSorted.append("SUCCESS")
        else:
            ConverterStatusSorted.append("FAIL")
            nb_errors += 1

        if (ZoneStatus[ZoneNamesSorted.index(Test)] == 0):
            ZoneStatusSorted.append("SUCCESS")
        else:
            ZoneStatusSorted.append("FAIL")
            nb_errors += 1

        if (OutputStatus[ZoneNamesSorted.index(Test)] == 0):
            OutputStatusSorted.append("SUCCESS")
        else:
            OutputStatusSorted.append("FAIL")
            nb_errors += 1

        InfoSorted.append(Info[ZoneNamesSorted.index(Test)])
        CompTimeSorted.append(round(CompTime[ZoneNamesSorted.index(Test)],5))

        print(Test + sw2dToolbox.repeatString(" ", MaxSizeZone - len(Test)) +
              "\t" + ConverterStatusSorted[idx] +
                       "\t" + ZoneStatusSorted[idx] +
                       "\t" + OutputStatusSorted[idx] +
                       "\t" + str(CompTimeSorted[idx]) +
                       "\t" + InfoSorted[idx])

    print("*************************")
    print("")
    
    return nb_errors, ZoneNamesSorted, ConverterStatusSorted, ZoneStatusSorted, OutputStatusSorted, InfoSorted, CompTimeSorted