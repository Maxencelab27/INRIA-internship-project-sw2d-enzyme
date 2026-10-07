import os
import shutil
import sys

sw2d_folder = os.getcwd().split("test")[0]

doc_folder = os.path.join(sw2d_folder, "doc/pdf/src")
ValidationDoc_folder = os.path.join(doc_folder, "validation_doc")

sys.path.insert(0, os.path.join(sw2d_folder, "utils/python"))
from personal_path import tab
import sw2dToolbox

slashTeX = "/"

def buildValidationDoc(Sw2dNames, ConverterStatus, Sw2dStatus, CompTime, PlotStatus, Info):
    print("***** BUILD VALIDATION DOC *****")
    
    # remove subfiles dir if it exists        
    if os.path.exists(os.path.join(ValidationDoc_folder, "subfiles")):
        shutil.rmtree(os.path.join(ValidationDoc_folder, "subfiles"))

    os.makedirs(os.path.join(ValidationDoc_folder, "subfiles"))

    fileSummaryTable = open(os.path.join(ValidationDoc_folder, "subfiles", "ValidationSummaryTable.tex"), 'w')
    fileSummaryTable.write("\\begin{center}\n")
    fileSummaryTable.write("\\begin{tabular}{l|c|c|c|c|c}\n")
    fileSummaryTable.write("Test & Conv. & SW2D & TIME & PLOT & INFO\\tabularnewline\n")
    fileSummaryTable.write("\\hline\n")
    for idx, Test in enumerate(Sw2dNames):
        if (idx % 59) == 58: # Newpage for the summary table
            fileSummaryTable.write("\\hline\n")
            fileSummaryTable.write("\\end{tabular} \n")
            fileSummaryTable.write("\\end{center}\n")
            fileSummaryTable.write("\\pagebreak\n")
            fileSummaryTable.write("\\begin{center}\n")
            fileSummaryTable.write("\\begin{tabular}{l|c|c|c|c|c}\n")
            fileSummaryTable.write("Test & Conv. & SW2D & TIME & PLOT & INFO\\tabularnewline\n")
            fileSummaryTable.write("\\hline\n")
        ref = reference(Test)
        Test = Test.replace('\\', slashTeX)
        Test = Test.replace('_', '\\_')
        fileSummaryTable.write("\\hyperref[" + ref + "]{" + Test + "} & " +
                                ConverterStatus[idx] + " & " +
                                Sw2dStatus[idx] + " & " +
                                str(CompTime[idx])[0:7] + " & " +
                                PlotStatus[idx] + " & " +
                                Info[idx] + "\\tabularnewline\n")

    fileSummaryTable.write("\\hline\n")
    fileSummaryTable.write("\\end{tabular} \n")
    fileSummaryTable.write("\\end{center}\n")
    fileSummaryTable.write("TOTAL TIME: " + str(sum(CompTime))[0:7])
    fileSummaryTable.close()


    Chapter = []
    LabelLvl = []
    Section = []
    SubSection = []
    SubSubSection = []

    for Test in Sw2dNames:
        SubFolders = Test.split(tab)
        #if (len(SubFolders) >= 1):
        #    if (SubFolders[0] == "fully_documented"):
        #        SubFolders.pop(0)
             
        Chap = ""
        Sec = ""
        SubSec = ""
        SubSubSec = ""
        
        if (len(SubFolders) >= 1):
            Chap = "Model: " + SubFolders[0].replace("\\_", "_")
            LabelLevel = 1
            
        if (len(SubFolders) >= 2):
            Sec = "Feature: " + SubFolders[1].replace("\\_", "_")
            LabelLevel = 2
            
        if (len(SubFolders) >= 3):
            SubSec = "Test case: " + SubFolders[2].replace("\\_", "_")
            LabelLevel = 3
            
        if (len(SubFolders) >= 4):
            SubSubSec = SubFolders[3].replace("\\_", "_")
            LabelLevel = 4
            
        if (len(SubFolders) >= 5):
            print("not implemented")
            exit()            
        
        Chapter.append(Chap)
        Section.append(Sec)
        SubSection.append(SubSec)
        SubSubSection.append(SubSubSec)
        LabelLvl.append(LabelLevel)
        
    writeDocLatex("FullTests.tex", "validation", Sw2dNames, LabelLvl, Chapter, Section, SubSection, SubSubSection)

    # Generate the pdf
    CurrentFolder = os.getcwd()
    dst = os.path.join(sw2d_folder, "doc/pdf", "sw2d_validationDoc.pdf")
    ## delete pdf
    deleteFile(dst)

    ## produce pdf
    os.chdir(doc_folder)
    os.system("pdflatex --interaction=batchmode --halt-on-error " + os.path.join(doc_folder, "validationDoc-main.tex"))
    os.system("pdflatex --interaction=batchmode --halt-on-error " + os.path.join(doc_folder, "validationDoc-main.tex")) # run twice for hyperlink generation

    ## delete compilation files
    deleteLatexCompilationFile("validationDoc-main")

    src = os.path.join(doc_folder, "validationDoc-main.pdf")
    shutil.move(src, dst)    
    os.chdir(CurrentFolder)
    print("*************************")
    
def cleanValidationDocSubFiles(version):
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
        
    if (version == "validation"):
        folder_name = "subfiles_fulldoc"
        subfiles_list = ["converterFullDoc.tex", \
                         "porosityMapperFullDoc.tex", \
                         "postProcessFullDoc.tex", \
                         "solverFullDoc.tex", \
                         "solverSummaryTable.tex", \
                         "solverTests.tex", \
                         "zoneFullDoc.tex"]
    elif (version == "merge"):
        folder_name = "subfiles_shortdoc"
        subfiles_list = ["converterShortDoc.tex", \
                         "porosityMapperShortDoc.tex", \
                         "postProcessShortDoc.tex", \
                         "solverShortDoc.tex", \
                         "solverSummaryTable.tex", \
                         "solverTests.tex", \
                         "zoneShortDoc.tex"]
        
    output_dir = os.path.join(ValidationDoc_folder, folder_name)

    # remove output dir if it exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # Create folder for output file if needed
    if (os.path.exists(output_dir) == False):
        os.makedirs(output_dir)
        
    ## check that compulsory files exist (otherwise create empty files)
    for file in subfiles_list:
        dst = os.path.join(ValidationDoc_folder, folder_name, file)
        
        sw2dToolbox.deleteFile(dst)
        file = open(dst, 'w')
        file.close()     

def compileValidationDoc(version):
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
        
    if (version == "validation"):
        print("***** BUILD VALIDATION DOC *****")
        pdf_file_name = "sw2d_validationDoc.pdf"
        latex_file_name = "validationDoc-main"
    elif (version == "merge"):
        print("***** BUILD SHORT VALIDATION DOC *****")
        pdf_file_name = "sw2d_mergeRequestDoc.pdf"
        latex_file_name = "mergeRequestDoc-main"
        
    
    # Generate the pdf
    CurrentFolder = os.getcwd()
    dst = os.path.join(sw2d_folder, "doc", "pdf", pdf_file_name)
    ## delete pdf
    sw2dToolbox.deleteFile(dst)

    ## produce pdf
    os.chdir(doc_folder)
    os.system("pdflatex --interaction=batchmode " + os.path.join(doc_folder, latex_file_name + ".tex"))
    os.system("pdflatex --interaction=batchmode " + os.path.join(doc_folder, latex_file_name + ".tex")) # run twice for hyperlink generation

    ## delete compilation files
    deleteLatexCompilationFile(latex_file_name)
    
    src = os.path.join(doc_folder, latex_file_name + ".pdf")
    shutil.move(src, dst)    
    os.chdir(CurrentFolder)
    print("*************************")
    
def deleteLatexCompilationFile(name):
    sw2dToolbox.deleteFile(name + ".aux")
    sw2dToolbox.deleteFile(name + ".bbl")
    sw2dToolbox.deleteFile(name + ".blg")
    sw2dToolbox.deleteFile(name + ".log")
    sw2dToolbox.deleteFile(name + ".out")
    sw2dToolbox.deleteFile(name + ".synctex")
    sw2dToolbox.deleteFile(name + ".toc")
        
def writeDocLatex(filename, version, Sw2dNames, LabelLvl, Chapter, Section, SubSection, SubSubSection):
    if not((version == "validation") or (version == "merge")):
        print("defined value for version" + version)
        print("authorized values for 'version' are:")
        print("  - full : for the complete validation doc")
        print("  - short : for the short validation doc (with only the figure)")
        exit()

    if (version == "validation"):
        fileLatex = open(os.path.join(ValidationDoc_folder, "subfiles_fulldoc", filename), "w")
    elif (version == "merge"):
        fileLatex = open(os.path.join(ValidationDoc_folder, "subfiles_shortdoc", filename), "w")            
    

    i = 0
    prec_chapter = slashTeX # a folder cannot be named with the tabulation symbol
    prec_section = slashTeX # a folder cannot be named with the tabulation symbol
    prec_Ssection = slashTeX # a folder cannot be named with the tabulation symbol
    prec_SSsection = slashTeX # a folder cannot be named with the tabulation symbol
    while i < len(Sw2dNames):
        ref = sw2dToolbox.reference(Sw2dNames[i])    
        label_chap = ""
        label_sec = ""
        label_Ssec = ""
        label_SSsec = ""
        if LabelLvl[i] == 1:
            label_chap = "\\label{" + ref + "}"
        if LabelLvl[i] == 2:
            label_sec = "\\label{" + ref + "}"
        if LabelLvl[i] == 3:
            label_Ssec = "\\label{" + ref + "}"
        if LabelLvl[i] == 4:
            label_SSsec = "\\label{" + ref + "}"
           
        if (Chapter[i] != prec_chapter):
            # new chapter
            prec_chapter = Chapter[i]
            prec_section = slashTeX # a folder cannot be named with the tabulation symbol
            prec_Ssection = slashTeX # a folder cannot be named with the tabulation symbol
            prec_SSsection = slashTeX # a folder cannot be named with the tabulation symbol
            fileLatex.write("\\chapter{" + label_chap + Chapter[i].replace("_", "\\_") + "}\n")
        
        if ((Section[i] != prec_section) and (len(Section[i]) > 0)):
            # new section
            prec_section = Section[i]
            prec_Ssection = slashTeX # a folder cannot be named with the tabulation symbol
            prec_SSsection = slashTeX # a folder cannot be named with the tabulation symbol
            fileLatex.write("\\section{" + label_sec + Section[i].replace("_", "\\_") + "}\n")
    
        if ((SubSection[i] != prec_Ssection) and (len(SubSection[i]) > 0)):
            # new subsection
            prec_Ssection = SubSection[i]
            prec_SSsection = slashTeX # a folder cannot be named with the tabulation symbol
            fileLatex.write("\\subsection{" + label_Ssec + SubSection[i].replace("_", "\\_") + "}\n")
    
        if ((SubSubSection[i] != prec_SSsection) and (len(SubSubSection[i]) > 0)):
            # new subsubsection
            prec_SSsection = SubSubSection[i]
            fileLatex.write("\\subsubsection{" + label_SSsec + SubSubSection[i].replace("_", "\\_") + "}\n")
        
        Sw2dNames[i] = Sw2dNames[i].replace(tab,"/")
        PngList = sw2dToolbox.getListPng(os.path.join(sw2d_folder, "test/merge/solver", Sw2dNames[i], "graph_scripts"))

        printed = 0
        for png in PngList:
            #print(".." + "/" + Sw2dNames[i] + "/" + "graph_scripts" + "/" + png)
            if (version == "merge"):
                path = os.path.join(sw2d_folder, "test/merge/solver", Sw2dNames[i], "doc/short-mr.tex").replace("\\","/")
            if (version == "validation"):
                path = os.path.join(sw2d_folder, "test/merge/solver", Sw2dNames[i], "doc/full-doc.tex").replace("\\","/")

            if os.path.exists(path) and printed == 0:
                texstring = "\\input{" + path + "}\n"
                fileLatex.write(texstring)
                printed = 1
                    
            fileLatex.write("\\begin{figure}[!h]\n")
            fileLatex.write("\\centering\n")
            fileLatex.write("\\includegraphics[width=0.8\\columnwidth]{" +
                            os.path.join(sw2d_folder, "test/merge/solver", Sw2dNames[i], "graph_scripts"+ slashTeX + png).replace("\\","/") + "}\n")
            fileLatex.write("\\caption{Test \\texttt{" + Sw2dNames[i].replace("_", "\\_") + "}}\n")
            fileLatex.write("\\end{figure}\n")
            fileLatex.write("\\clearpage\n")

        fileLatex.write("\n\n\n")
        i = i + 1
    
    fileLatex.close()

def writeLatexConverter(version, list_tests, execution_result):
    stop = False
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
    if (stop == True):
        quit()         

    if (version == "validation"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_fulldoc", "converterFullDoc.tex"), "w")
    elif (version == "merge"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_shortdoc", "converterShortDoc.tex"), "w")          
    nb_errors = 0

    file.write("\\part{sw2dConverter tests}\n")
    
    file.write("\\begin{center}\n")
    file.write("\\begin{tabular}{l|c|l}\n")
    file.write("Test & EXECUTION RESULT & INFO\\tabularnewline\n")
    file.write("\\hline\n")    
    
    for i, test in enumerate(list_tests):
        test_name = test.get("Folder").replace('\\',slashTeX)
        test_name = test_name.replace('_', '\\_')        
        file.write(test_name + " & " + execution_result[i] + " & " + test.get("Info") + "\\tabularnewline\n")
        res = 0
        if (execution_result[i] != "Success"):
            res = test.get("Result")
        else :
            res = 1 - test.get("Result")
        nb_errors += res
    
    file.write("\\hline\n")
    file.write("\\end{tabular} \n")
    file.write("\\end{center}\n")
    file.write("NUMBER OF ERRORS: " + str(nb_errors))
    
    file.close()
    
def writeLatexPorosityMapper(version, list_tests, converter_status, mapper_status, solver_status, output_status):
    stop = False
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
    if (stop == True):
        quit() 
        
    if (version == "validation"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_fulldoc", "porosityMapperFullDoc.tex"), "w")
    elif (version == "merge"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_shortdoc", "porosityMapperShortDoc.tex"), "w")  
        
    file.write("\\part{sw2dPorosityMapper tests}\n")
    
    file.write("\\begin{center}\n")
    file.write("\\begin{tabular}{l|c|c|c|c|l}\n")
    file.write("TEST & CONVERTER & MAPPER & SOLVER & OUTPUT & INFO\\tabularnewline\n")
    file.write("\\hline\n")
    
    nb_errors = 0
    for i, test in enumerate(list_tests):
        test_name = test.get("Folder").replace('\\',slashTeX)
        test_name = test_name.replace('_', '\\_')  
        file.write(test_name + " & " + converter_status[i] + " & " + mapper_status[i] + " & " + solver_status[i] + " & " + output_status[i] + " & " + test.get("Info") + "\\tabularnewline\n")
        res = 0
        if (not (converter_status[i] == "Success" and solver_status[i] == "Success" and mapper_status[i] == "Success" and output_status[i] == "Success")):
            res = test.get("Result")
        else :
            res = 1 - test.get("Result")
        nb_errors += res
    
    file.write("\\hline\n")
    file.write("\\end{tabular} \n")
    file.write("\\end{center}\n")
    file.write("NUMBER OF ERRORS: " + str(nb_errors))    
    
    file.close()
    
def writeLatexPostProcess(version):
    stop = False
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
    if (stop == True):
        quit() 
        
    if (version == "validation"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_fulldoc", "postProcessFullDoc.tex"), "w")
    elif (version == "merge"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_shortdoc", "postProcessShortDoc.tex"), "w")  
        
    file.write("\\part{sw2dPostProcess tests}\n")
    file.write("This part has not been implemented yet.\n")
    
    file.close()
        
        
def writeLatexSolver(Sw2dNames, ConverterStatus, Sw2dStatus, CompTime, PlotStatus, Info, version):
    if (version == "validation"):
        subfiles_folder = "subfiles_fulldoc"
        root_latex_file = "solverFullDoc.tex"
    elif (version == "merge"):
        subfiles_folder = "subfiles_shortdoc"
        root_latex_file = "solverShortDoc.tex"
        
    file = open(os.path.join(ValidationDoc_folder, subfiles_folder, root_latex_file), 'w')
    file.write("\\part{sw2dSolver tests}\n")
    file.close()    
    
    fileSummaryTable = open(os.path.join(ValidationDoc_folder, subfiles_folder, "solverSummaryTable.tex"), 'w')
    fileSummaryTable.write("\\begin{center}\n")
    fileSummaryTable.write("\\begin{tabular}{l|c|c|c|c|c}\n")
    fileSummaryTable.write("Test & Conv. & SW2D & TIME & PLOT & INFO\\tabularnewline\n")
    fileSummaryTable.write("\\hline\n")
    for idx, Test in enumerate(Sw2dNames):
        if (idx % 59) == 58: # Newpage for the summary table
            fileSummaryTable.write("\\hline\n")
            fileSummaryTable.write("\\end{tabular} \n")
            fileSummaryTable.write("\\end{center}\n")
            fileSummaryTable.write("\\pagebreak\n")
            fileSummaryTable.write("\\begin{center}\n")
            fileSummaryTable.write("\\begin{tabular}{l|c|c|c|c|c}\n")
            fileSummaryTable.write("Test & Conv. & SW2D & TIME & PLOT & INFO\\tabularnewline\n")
            fileSummaryTable.write("\\hline\n")
        ref = sw2dToolbox.reference(Test)
        Test = Test.replace('\\',slashTeX)
        Test = Test.replace('_', '\\_')
        fileSummaryTable.write("\\hyperref[" + ref + "]{" + Test + "} & " +
                                ConverterStatus[idx] + " & " +
                                Sw2dStatus[idx] + " & " +
                                str(CompTime[idx])[0:7] + " & " +
                                PlotStatus[idx] + " & " +
                                Info[idx] + "\\tabularnewline\n")

    fileSummaryTable.write("\\hline\n")
    fileSummaryTable.write("\\end{tabular} \n")
    fileSummaryTable.write("\\end{center}\n")
    fileSummaryTable.write("TOTAL TIME: " + str(sum(CompTime))[0:7])
    fileSummaryTable.close()   
    
    Chapter = []
    LabelLvl = []
    Section = []
    SubSection = []
    SubSubSection = []
    for Test in Sw2dNames:
        SubFolders = Test.split(tab)
        if (len(SubFolders) >= 1):
            if (SubFolders[0] == "in_progress"):
                SubFolders.pop(0)
                
        Chap = ""
        Sec = ""
        SubSec = ""
        SubSubSec = ""
        
        if (len(SubFolders) >= 1):
            Chap = "Model: " + SubFolders[0].replace("\\_", "_")
            LabelLevel = 1
            
        if (len(SubFolders) >= 2):
            Sec = "Feature: " + SubFolders[1].replace("\\_", "_")
            LabelLevel = 2
            
        if (len(SubFolders) >= 3):
            SubSec = "Test case: " + SubFolders[2].replace("\\_", "_")
            LabelLevel = 3
            
        if (len(SubFolders) >= 4):
            SubSubSec = SubFolders[3].replace("\\_", "_")
            LabelLevel = 4
            
        if (len(SubFolders) >= 5):
            print("not implemented")
            exit()            
            
        Chapter.append(Chap)
        Section.append(Sec)
        SubSection.append(SubSec)
        SubSubSection.append(SubSubSec)
        LabelLvl.append(LabelLevel)

    writeDocLatex("solverTests.tex", version, Sw2dNames, LabelLvl, Chapter, Section, SubSection, SubSubSection)  
    
def writeLatexZone(version, nb_errors, ZoneNamesSorted, ConverterStatusSorted, ZoneStatusSorted, OutputStatusSorted, InfoSorted, CompTimeSorted):
    stop = False
    if ((version != "validation") and (version != "merge")):
        print("parameter 'version' have to be: 'full' ou 'short'")
        stop = True
    if (stop == True):
        quit() 
        
    if (version == "validation"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_fulldoc", "zoneFullDoc.tex"), "w")
    elif (version == "merge"):
        file = open(os.path.join(ValidationDoc_folder, "subfiles_shortdoc", "zoneShortDoc.tex"), "w")  
        
    file.write("\\part{sw2dZone tests}\n")
    
    file.write("\\begin{center}\n")
    file.write("\\begin{tabular}{l|c|c|c|c|l}\n")
    file.write("TEST & CONVERTER & ZONE & OUTPUT & TIME & INFO\\tabularnewline\n")
    file.write("\\hline\n")
    
    for i, test_name in enumerate(ZoneNamesSorted):
        test_name = test_name.replace('\\',slashTeX)
        test_name = test_name.replace('_', '\\_')  
        file.write(test_name + " & " + ConverterStatusSorted[i] + " & " + ZoneStatusSorted[i] + " & " + OutputStatusSorted[i] + " & " + str(CompTimeSorted[i]) + " & " + InfoSorted[i] + "\\tabularnewline\n")
    
    file.write("\\hline\n")
    file.write("\\end{tabular} \n")
    file.write("\\end{center}\n")
    file.write("NUMBER OF ERRORS: " + str(nb_errors))    
    
    file.close()