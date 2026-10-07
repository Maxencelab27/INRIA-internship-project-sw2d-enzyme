import math
import os
import shutil
import sys
import time

sw2dFolder = os.getcwd().split("utils")[0]
sys.path.insert(0, sw2dFolder + "utils/python")

from personal_path import cmd_clear, sw2dConverter, sw2dSolver

sep = "*****************************"

def processFolder(Dossier, MaxDifference):
    # define the storing table
    Names = []
    ConverterStatus = []
    SolverStatus = []
    OutputConvStatus = []
    OutputStatus = []
    Info = []
    CompTime = []

    FolderList = GetFolderList(Dossier)
    
    for Test in FolderList:
        if os.path.exists(os.path.join(Dossier, Test, "primary_files")):
            # primary_files subfolder exists => running test
            print("    ", sep, "test:", Test)
            FolderNames = []
            FolderConverterStatus = []
            FolderSolverStatus = []
            FolderOutputConvStatus = []
            FolderOutputStatus = []
            FolderInfo = []
            FolderCompTime = []
            t_0 = time.time()
            a, b, c, d, e = Compute(Dossier, Test, MaxDifference)
            t = time.time() - t_0
            FolderNames.append(os.path.join(Dossier, Test)[1 + len(os.path.join(sw2dFolder, "test/merge/output_converter")):])
            FolderConverterStatus.append(a)
            FolderSolverStatus.append(b)
            FolderOutputConvStatus.append(c)
            FolderOutputStatus.append(d)
            FolderInfo.append(e)
            FolderCompTime.append(t)

        else:
            # running Python
            FolderNames = []
            FolderConverterStatus = []
            FolderSolverStatus = []
            FolderOutputConvStatus = []
            FolderOutputStatus = []
            FolderInfo = []
            FolderCompTime = []
            src = os.path.join(sw2dFolder, "test/merge/output_converter", "run_mr-tests.py")
            dst = os.path.join(Dossier, Test, "run_mr-tests.py")
            shutil.copyfile(src, dst)
            print("======================")
            print("  subfolder:", os.path.join(Dossier, Test))
            FolderNames, FolderConverterStatus, FolderSolverStatus, FolderOutputConvStatus, FolderOutputStatus, FolderInfo, FolderCompTime = processFolder(os.path.join(Dossier, Test), MaxDifference)

        for name in FolderNames:
            Names.append(name)
            
        for status in FolderConverterStatus:
            ConverterStatus.append(status)
            
        for status in FolderSolverStatus:
            SolverStatus.append(status)

        for status in FolderOutputConvStatus:
            OutputConvStatus.append(status)
            
        for status in FolderOutputStatus:
            OutputStatus.append(status)

        for info in FolderInfo:
            Info.append(info)
            
        for t in FolderCompTime:
            CompTime.append(t)

    return Names, ConverterStatus, SolverStatus, OutputConvStatus, OutputStatus, Info, CompTime


def RunTest(Dossier, MaxDifference):
    os.system(cmd_clear)
    print("run Test on:", Dossier)
    Names, ConverterStatus, SolverStatus, OutputConvStatus, OutputStatus, Info, CompTime = processFolder(Dossier, MaxDifference)

    NamesSorted = sorted(Names)
    ConverterStatusSorted = []
    SolverStatusSorted = []
    OutputConvStatusSorted = []
    OutputStatusSorted = []
    InfoSorted = []
    CompTimeSorted = []

    # print the summary
    MaxSizeZone = max([len(Test) for Test in Names])
    if not ((MaxSizeZone % 2) == 0): # non divisible par 2
        MaxSizeZone = MaxSizeZone + 1

    print("\n")
    print("***** TESTS SUMMARY *****")

    if (sum(OutputConvStatus) == 0):
        print("- output converter tests processed successfully")
    else:
        print("- output converter tests encountered errors")

    if (sum(OutputStatus) == 0):
        print("- output files are correct considering a maximum gap of " + str(MaxDifference))
    else:
        print("- error in the output files considering a maximum gap of " + str(MaxDifference))
    
    print("- TOTAL TIME: " + str(sum(CompTime)))
    print(" ")

    print(RepeatString(" ", (MaxSizeZone - 4) / 2) + "TEST" +
          RepeatString(" ", (MaxSizeZone - 4) / 2) + "\t" +
          '{0: <10}'.format("| CONV.") +
          '{0: <10}'.format("| SOLVER") +
          '{0: <15}'.format("| OUTPUT CONV.") +
          '{0: <10}'.format("| OUTPUT") +
          '{0: <13}'.format("| COMP. TIME") +
          "| INFO")

    
    nb_errors = 0

    for idx, Test in enumerate(NamesSorted):
        if (ConverterStatus[NamesSorted.index(Test)] == 0):
            ConverterStatusSorted.append("SUCCESS")
        else:
            ConverterStatusSorted.append("FAIL")
            nb_errors += 1
            
        if (SolverStatus[NamesSorted.index(Test)] == 0):
            SolverStatusSorted.append("SUCCESS")
        else:
            SolverStatusSorted.append("FAIL")
            nb_errors += 1

        if (OutputConvStatus[NamesSorted.index(Test)] == 0):
            OutputConvStatusSorted.append("SUCCESS")
        else:
            OutputConvStatusSorted.append("FAIL")
            nb_errors += 1
            
        if (OutputStatus[NamesSorted.index(Test)] == 0):
            OutputStatusSorted.append("SUCCESS")
        else:
            OutputStatusSorted.append("FAIL")
            nb_errors += 1

        InfoSorted.append(Info[NamesSorted.index(Test)])
        CompTimeSorted.append(CompTime[NamesSorted.index(Test)])
        
        print(Test + RepeatString(" ", MaxSizeZone - len(Test)) + "\t"
              '{0: <10}'.format("| " + ConverterStatusSorted[idx]) +
              '{0: <10}'.format("| " + SolverStatusSorted[idx]) +
              '{0: <15}'.format("| " + OutputConvStatusSorted[idx]) +
              '{0: <10}'.format("| " + OutputStatusSorted[idx]) +
              '{0: <13}'.format("| " + str(CompTimeSorted[idx])[0:6]) +
              "| " + InfoSorted[idx])

    print("*************************")
    print("")

    return(nb_errors)
    

def getRefFilesList(Folder):
    if os.path.exists(Folder):
        filelist = os.listdir(Folder)
    else:
        filelist = []

    return filelist

def GetFolderList(Folder):
    # get list of tests
    SubFolderList = [ f.path[len(Folder):] for f in os.scandir(Folder) if f.is_dir() ]
    SubFolderList.sort()

    ## remove possible slashs
    for idx, test in enumerate(SubFolderList):
        # ending slash
        if ((test[-1] == "/") or
            (test[-1] == "\\")):
            test = test[:-1]
            SubFolderList[idx] = test
        # starting slash
        if ((test[0] == "/") or
            (test[0] == "\\")):
            test = test[1:]
            SubFolderList[idx] = test

    return SubFolderList

def RepeatString(String, Nb):
    return (String * math.floor(Nb))

def checkStatus(path, MaxDifference):
    Status = 0
    RefFolder = os.path.join(path, "reference")
    OutputFolder = os.path.join(path, "output")
    refFiles = getRefFilesList(RefFolder)

    for fich in refFiles:
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
                        if (i==j and i>1) : # When reading the same line, let us compare them
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
                                        print(fich,"difference for line", i,", column", c, ":", references[c], outputs[c], abs(reference_value-output_value))
                                        break
            # closing files
            refStream.close()
            outputStream.close()
        
    return Status

def Compute(Root, Test, MaxDifference):
    Folder = os.path.join(Root, Test)
    print("    Running sw2dZone for", Folder)

    stop = False
    sw2dOutputConverter = sw2dFolder + "utils/python/apps/sw2dOutputConverter/main.py"
    if not(os.path.exists(sw2dConverter)):
        print("sw2dConverter not found", sw2dConverter)
        stop = True
    if not(os.path.exists(sw2dSolver)):
        print("sw2dSolver not found", sw2dSolver)
        stop = True
    if not(os.path.exists(sw2dOutputConverter)):
        print("sw2dOutputConverter not found", sw2dOutputConverter)
        stop = True
    if (stop ==True):
        print("executables were not found. Abort.")
        quit()    
        
    ConverterStatus = -3
    SolverStatus = -3
    OuputConvStatus = -3
    OutputStatus = -3
    Info = ""

    # check if test exists
    if not(os.path.exists(Folder)):
        ConverterStatus = -2
        SolverStatus = -2
        OuputConvStatus = -2
        OutputStatus = -2
        Info = "Test case not found"
    else:
        primary_dir = os.path.join(Folder, "primary_files")
        input_dir = os.path.join(Folder, "input")

        if not(os.path.exists(primary_dir)):
            ConverterStatus = -1
            SolverStatus = -1
            OuputConvStatus = -1
            Info = "Missing primary files"
        else:
            # Running converter
            ## Removing previous output files
            if (os.path.exists(os.path.join(input_dir, "geometry.sw2d")) == True):
                os.remove(os.path.join(input_dir, "geometry.sw2d"))

            Command = sw2dConverter + " " + Folder
            ConverterStatus = os.system(Command)
            if not(ConverterStatus == 0):
                Info = "sw2dConverter failed"
            else:
            # Running solver
                Command = sw2dSolver + " " + Folder
                SolverStatus = os.system(Command)

                if not(SolverStatus == 0):
                    Info = "sw2dSolver failed"
                else:
                # Running sw2dOutputConverter
                    Command = "python " + sw2dOutputConverter + " " + Folder +  " text"
                    OuputConvStatus = os.system(Command)
                
                    if not(OuputConvStatus == 0):
                        Info = "sw2dOuputConv failed"
                    else :
                        OutputStatus = checkStatus(Folder, MaxDifference)
                        if not(OutputStatus == 0) :
                            Info = "Output files are different from reference"

    print("    ", sep)

    return ConverterStatus, SolverStatus, OuputConvStatus, OutputStatus, Info