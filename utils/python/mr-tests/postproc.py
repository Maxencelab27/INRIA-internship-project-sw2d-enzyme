import os
import sys
##from datetime import datetime
#import math
##import csv
import re
import shutil
import time

sw2d_folder = os.getcwd().split("test")[0]
sys.path.insert(0, os.path.join(sw2d_folder,"utils/python"))
import sw2dToolbox

from personal_path import cmd_clear, sw2dConverter, sw2dSolver, sw2dPostProcess

def processFolder(Dossier, MaxDifference, recursive = True):
    # define the storing table
    Names = []
    ConverterStatus = []
    SolverStatus = []
    PostProcStatus = []
    OutputStatus = []
    PlotStatus = []
    Info = []
    CompTime = []
    total = 0

    FolderList = sw2dToolbox.getFolderList(Dossier)
    
    for Test in FolderList:
        parent = os.path.basename(os.path.dirname(os.path.join(Dossier, Test)))
        if (os.path.exists(os.path.join(Dossier, Test, "input"))):
            # primary_files subfolder exists => running sw2dZone
            print("    ", sw2dToolbox.sep, "test:", Test)
            ListNames = []
            ListConverterStatus = []
            ListSolverStatus = []
            ListPostProcStatus = []
            ListOutputStatus = []
            ListInfo = []
            ListCompTime = []
            t_0 = time.time()
            if (parent == "submersion_durations") :
                a, b, e, c, d = ComputeSubmersion(Dossier, Test, MaxDifference)
            else :
                a, b, e, c, d = Compute(Dossier, Test, MaxDifference)
            t = time.time() - t_0
            ListNames.append(os.path.join(Dossier, Test)[1 + len(os.path.join(sw2d_folder, "test/merge/postprocess")):])
            ListConverterStatus.append(a)
            ListSolverStatus.append(b)
            ListPostProcStatus.append(e)
            ListOutputStatus.append(c)
            ListInfo.append(d)
            ListCompTime.append(t)

        else:
            # running Python
            ListNames = []
            ListConverterStatus = []
            ListSolverStatus = []
            ListPostProcStatus = []
            ListInfo = []
            ListCompTime = []
            src = os.path.join(sw2d_folder, "test/merge/postprocess", "testPostProcess.py")
            dst = os.path.join(Dossier, Test, "testPostProcess.py")
            shutil.copyfile(src, dst)
            print("======================")
            print("  subfolder:", os.path.join(Dossier, Test))
            ListNames, ListConverterStatus, ListSolverStatus, ListPostProcStatus, ListOutputStatus, ListInfo, ListCompTime = processFolder(os.path.join(Dossier, Test), MaxDifference)

        for name in ListNames:
            Names.append(name)
            
        for name in ListConverterStatus:
            ConverterStatus.append(name)
            
        for name in ListSolverStatus:
            SolverStatus.append(name)
            
        for name in ListPostProcStatus:
            PostProcStatus.append(name)
            
        for name in ListOutputStatus:
            OutputStatus.append(name)

        for name in ListInfo:
            Info.append(name)
            
        for name in ListCompTime:
            CompTime.append(name)

    return Names, ConverterStatus, SolverStatus, PostProcStatus, OutputStatus, Info, CompTime


def runPostProcTest(Dossier, MaxDifference):
    os.system(cmd_clear)
    print("run Test on:", Dossier)
    Names, ConverterStatus, SolverStatus, PostProcStatus, OutputStatus, Info, CompTime = processFolder(Dossier, MaxDifference)

    NamesSorted = sorted(Names)
    ConverterStatusSorted = []
    SolverStatusSorted = []
    PostProcStatusSorted = []
    OutputStatusSorted = []
    InfoSorted = []
    CompTimeSorted = []

    # print the summary
    MaxSizeZone = max([len(Test) for Test in Names])
    if not ((MaxSizeZone % 2) == 0): # non divisible par 2
        MaxSizeZone = MaxSizeZone + 1

    print("\n")
    print("***** TESTS SUMMARY *****")

    if (sum(ConverterStatusSorted) + sum(SolverStatusSorted) == 0):
        print("- sw2d tests processed successfully")
    else:
        print("- sw2d tests encountered errors")

    if (sum(PostProcStatusSorted) == 0):
        print("- PostProcessor run without error")
    else:
        print("- PostProcessor run with error")
        
    if (sum(OutputStatus) == 0):
        print("- output files are correct considering a maximum gap of " + str(MaxDifference))
    else:
        print("- error in the output files considering a maximum gap of " + str(MaxDifference))
    
    print("- TOTAL TIME: " + str(sum(CompTime)))
    print(" ")

    print(sw2dToolbox.repeatString(" ", (MaxSizeZone - 4) / 2) + "TEST" +
          sw2dToolbox.repeatString(" ", (MaxSizeZone - 4) / 2) +
          "\tConv.\tSolv.\tPostPr.\tOutput\tTIME\tINFO")

    nb_errors = 0

    for idx, Test in enumerate(NamesSorted):
        nb_test_error = 0
        if (Test[0:20] == "submersion_durations") :
            ConverterStatusSorted.append("(OFF)")
            SolverStatusSorted.append("(OFF)")
                
            if (PostProcStatus[NamesSorted.index(Test)] == 0):
                PostProcStatusSorted.append("SUCCESS")
            else:
                PostProcStatusSorted.append("FAIL")
                nb_test_error += 1
                
            if (OutputStatus[NamesSorted.index(Test)] == 0):
                OutputStatusSorted.append("SUCCESS")
            else:
                OutputStatusSorted.append("FAIL")
                nb_test_error += 1

        else :
            if (ConverterStatus[NamesSorted.index(Test)] == 0):
                ConverterStatusSorted.append("SUCCESS")
            else:
                ConverterStatusSorted.append("FAIL")
                nb_test_error += 1
                
            if (SolverStatus[NamesSorted.index(Test)] == 0):
                SolverStatusSorted.append("SUCCESS")
            else:
                SolverStatusSorted.append("FAIL")
                nb_test_error += 1
                
            if (PostProcStatus[NamesSorted.index(Test)] == 0):
                PostProcStatusSorted.append("SUCCESS")
            else:
                PostProcStatusSorted.append("FAIL")
                nb_test_error += 1
                
            if (OutputStatus[NamesSorted.index(Test)] == 0):
                OutputStatusSorted.append("SUCCESS")
            else:
                OutputStatusSorted.append("FAIL")
                nb_test_error += 1
                
        if nb_test_error > 0:
            nb_errors += 1

        InfoSorted.append(Info[NamesSorted.index(Test)])
        CompTimeSorted.append(CompTime[NamesSorted.index(Test)])

        print(Test + sw2dToolbox.repeatString(" ", MaxSizeZone - len(Test)) +
                       "\t" + ConverterStatusSorted[idx] +
                       "\t" + SolverStatusSorted[idx] +
                       "\t" + PostProcStatusSorted[idx] +
                       "\t" + OutputStatusSorted[idx] +
                       "\t" + str(CompTimeSorted[idx])[0:6] +
                       "\t" + InfoSorted[idx])

    print(sw2dToolbox.sep)
    print("")

    return(nb_errors)


def getRefFilesList(Folder):
    filelist = []
    if os.path.exists(Folder):
        for file in os.listdir(Folder):
            if file.endswith(".txt"):
                filelist.append(file)        

    return filelist
    
def checkStatus(path, MaxDifference):
    Status = 0
    RefFolder = os.path.join(path, "reference")
    OutputFolder = os.path.join(path, "output")
    refFiles = getRefFilesList(RefFolder)
    for fich in refFiles: # could be 'ResultHydro_Max.txt' or 'ResultHydro_DateOfMax.txt' 
        # Load output and reference files
        refFile = os.path.join(RefFolder, fich)
        outputFile = os.path.join(OutputFolder, fich)
        if not(os.path.isfile(os.path.join(OutputFolder, fich))):
            Status -= 1
            print("Output file not found")
        else:
            refStream = open(refFile, "r")
            outputStream = open(outputFile, "r")
            ref_data = refStream.readlines()
            out_data = outputStream.readlines()

            if (len(ref_data) != len(out_data)) :
                print ("Different number of lines for files " + refFile + " and " + outputFile)
                Status -= 1

            else : 
                i = 0
                for refline in ref_data: # Browse every line of reference file
                    i += 1
                    j = 0
                    for outline in out_data: # Browse every line of output file
                        j += 1
                        if i==j & i>1: # When reading the same line, let us compare them
                            references = refline.split()
                            outputs = outline.split()
                            if (len(references) != len(outputs)): # Check that lines have the same number of entries
                                Status -= 1
                                print ("Different number of entry for line" + str(i))
                            else:
                                for c in range(3): # Browse columns 0, 1 and 2 (x, y and zb) first
                                    reference_value = float(references[c])
                                    output_value = float(outputs[c])
                                    if ( abs(reference_value-output_value) > MaxDifference ):
                                        Status -= 1
                                        print(fich,"KO for line", i,", column", c + 1, ":", references[c], outputs[c], outputs)
                                for c in range(3,len(references)): # Browse all other entries, skipping the x, y and zb columns
                                    if ((fich=="ResultHydro_DateOfMax.txt") and (references[c] != outputs[c])): 
                                        # If times are different, we need to check that max values are OK
                                        filenameref = FindFileAtTime(OutputFolder, float(references[c]))
                                        # Remark : could be filenameref = FindFileAtTime(OutputFolder, 'ResultHydro_Max.txt')
                                        if (filenameref != "") :
                                            filenameref = os.path.join(OutputFolder, filenameref)
                                            fileStreamref = open(filenameref, "r")
                                            file_data_ref = fileStreamref.readlines()
                                            k=0
                                            # Goto line i, column c of both reference & output files and check the corresponding entries are identical (hence a max)
                                            for linefile in file_data_ref:
                                                k += 1
                                                if k==i:
                                                    values = linefile.split()
                                                    reference_value = float(values[c])
                                            fileStreamref.close()
                                        else :
                                            reference_value = None
                                        filenameout = FindFileAtTime(OutputFolder, float(outputs[c]))
                                        if (filenameout != "") :
                                            filenameout = os.path.join(OutputFolder, filenameout)
                                            fileStreamout = open(filenameout, "r")
                                            file_data_out = fileStreamout.readlines()
                                            k=0
                                            for linefile in file_data_out:
                                                k += 1
                                                if k==i:
                                                    values = linefile.split()
                                                    output_value = float(values[c])
                                            fileStreamout.close()
                                        else :
                                            output_value = None
                                    else: # if files are 'ResultHydro_Max.txt' or if times are identical in 'ResultHydro_TimeOfMax.txt'
                                        reference_value = float(references[c])
                                        output_value = float(outputs[c])
                                    if ((output_value == None) or (reference_value == None) or (abs(reference_value-output_value) > MaxDifference)):
                                        Status -= 1
                                        print(fich,"KO for line", i,", column", c + 1, ":", reference_value, output_value)
                            break
            # closing files
            refStream.close()
            outputStream.close()
        
    return Status

def checkStatusSubmersion(path):
    Status = 0
    RefFolder = os.path.join(path, "reference")
    OutputFolder = os.path.join(path, "output")
    refFiles = getRefFilesList(RefFolder)
    
    refFile = os.path.join(RefFolder, "SubmersionDurations.txt")
    outputFile = os.path.join(OutputFolder, "SubmersionDurations.txt")
    if not(os.path.isfile(outputFile)):
        Status -= 1
    else:
        refContent = open(refFile, "r").read()
        inputContent = open(outputFile, "r").read()
        if not(refContent.split() == inputContent.split()):
            Status -= 1
            
    return Status

def FindFileAtTime(Folder, Time): # Returns filename corresponding to $Time$ in the $Folder$ directory (could be in ResultHydro or in ResultTime)
    result = ""
    iso_date = sw2dToolbox.timeToISO8601(Time)
    file_at_time = os.path.join(Folder, "ResultHydro_"+iso_date+".txt")
    if os.path.isfile(file_at_time):
        result = file_at_time
    else:
        searched_result = "ResultTime_"+iso_date
        for f in os.listdir(Folder):
            if re.match(searched_result, f):
                result = f
    return result


def Compute(Root, Test, MaxDifference):
    Folder = os.path.join(Root, Test)
    print("    Running sw2dPostProcessor for", Test)
    
    # define the sw2dSolver path
    stop = False
    if not(os.path.exists(sw2dConverter)):
        print("sw2dConverter not found", sw2dConverter)
        stop = True
    if not(os.path.exists(sw2dSolver)):
        print("sw2dSolver not found", sw2dSolver)
        stop = True
    if not(os.path.exists(sw2dPostProcess)):
        print("sw2dPostProcess not found", sw2dPostProcess)
        stop = True
    if (stop ==True):
        print("executables was not found. Abort.")
        quit() 
        
    converterStatus = -3
    solverStatus = -3
    postProcStatus = -3
    OutputStatus = -3
    Info = ""

    # check if test exists
    if not(os.path.exists(Folder)):
        converterStatus = -2
        solverStatus = -2
        postProcStatus = -2
        OutputStatus = -2
        Info = "Test case not found"
    else:
        primary_dir = os.path.join(Folder, "primary_files")
        input_dir = os.path.join(Folder, "input")

        if not(os.path.exists(primary_dir)):
            converterStatus = -1
            Info = "Missing input files"
        else:
            # Running converter
            ## Removing previous output files
            if (os.path.exists(os.path.join(input_dir, "geometry.sw2d")) == True):
                os.remove(os.path.join(input_dir, "geometry.sw2d"))
            ## Run
            Command = sw2dConverter + " " + Folder
            converterStatus = os.system(Command)
            
            if not(converterStatus == 0):
                Info = "sw2dConverter failed"
            else:
            # Running solver
                print("    Start with running sw2dSolver to generate results in", Test)
                Command = sw2dSolver + " " + Folder
                solverStatus = os.system(Command)	    #### equal to zero = SUCCESS    
                
                # Running post process
                if (solverStatus == 0):
                    Command = sw2dPostProcess + " " + Folder
                    postProcStatus = os.system(Command)	    #### equal to zero = SUCCESS
                    
                    if (postProcStatus == 0):
                        OutputStatus = checkStatus(Folder, MaxDifference)

    print("    ", sw2dToolbox.sep)

    return (converterStatus, solverStatus, postProcStatus, OutputStatus, Info)

def ComputeSubmersion(Root, Test, MaxDifference):
    Folder = os.path.join(Root, Test)
    print("    Running sw2dPostProcessor for", Test)
    
    # define the sw2dSolver path
    if not(os.path.exists(sw2dPostProcess)):
        print("sw2dPostProcess not found", sw2dPostProcess)
        print("executable was not found. Abort.")
        quit() 

    postProcStatus = -3
    OutputStatus = -3
    Info = ""

    # check if test exists
    if not(os.path.exists(Folder)):
        postProcStatus = -2
        OutputStatus = -2
        Info = "Test case not found"
    else:
        output_dir = os.path.join(Folder, "output")

        if not(os.path.exists(output_dir)):
            postProcStatus = -1
            Info = "Missing output_dir files"
        else:
            # Running post process
            ## Removing previous output files
            if (os.path.exists(os.path.join(output_dir, "SubmersionDurations.txt")) == True):
                os.remove(os.path.join(output_dir, "SubmersionDurations.txt"))
            ## Run
            Command = sw2dPostProcess + " " + Folder
            postProcStatus = os.system(Command)
            
            if not(postProcStatus == 0):
                Info = "sw2dPostProcess failed"
            else:
                OutputStatus = checkStatusSubmersion(Folder)

    print("    ", sw2dToolbox.sep)

    return (0, 0, postProcStatus, OutputStatus, Info)