import csv
import matplotlib.pyplot as plt
import plotCommon
import numpy as np
from matplotlib import cm


def XYmap_FilledContour(**kwargs):
   # declare default value
    DefinedSim = False
    DefinedRef = False
    DefineAspect = False
    DefineContour = False
    DefineXLimits = False
    DefineLimits = False
    file = ""
    Language = "EN"
    Variable = []
    contour_number = 50

    # read argument
    for key, value in kwargs.items():
        if (key == "aspect"):
            # Define the output filename
            aspect = value
            DefineAspect = True
        elif (key == "filename"):
            # Define the output filename
            file = value
        elif (key == "ref"):
            # reference value
            ResultRef = value
            DefinedRef = True
        elif (key == "sim"):
            # SW2D computed value
            ResultSim = value
            DefinedSim = True
        elif (key == "variable"):
            # Table of the plotted variables
            Variable = value
        elif (key == "number_of_contours"):
            # Number of contour levels in the 2D plots
            contour_number = value
            DefineContour = True
        elif (key == "xlimits"):
            # Min and max values for the 2D plots
            xmin = value[0]
            xmax = value[1]
            DefineXLimits = True
        elif (key == "zlimits"):
            # Min and max values for the 2D plots
            zmin = value[0]
            zmax = value[1]
            DefineContour = True
            DefineLimits = True
        else:
            print("unknown argument:", key)
            exit()

    # check compulsory variables
    if not ((DefinedSim == True) and (DefinedRef == True)):
        print("Data have to be provided using 'ref' and 'sim'")
        exit()
    if (len(Variable) == 0):
        print("'variable' is a compulsory parameter")
        exit()

    # Clear figure
    plt.clf()

    # plot
    NbPlot = len(Variable)
    
    for i in range(NbPlot):
        fig, (plt1, plt2) = plt.subplots(2)
        plt.subplots_adjust(hspace = 0.5)
        # Plots
        if (DefineAspect == True):
            plt1.set_aspect(aspect)
            plt2.set_aspect(aspect)
        if (DefineXLimits == True):
            # Pick the min and max from the reference data
            plt1.set_xlim([xmin,xmax])
            plt2.set_xlim([xmin,xmax])
        if (DefineLimits == False):
            # Pick the min and max from the reference data
            zmin = np.min(ResultRef.getVar(Variable[i]))
            zmax = np.max(ResultRef.getVar(Variable[i]))
        if (DefineContour == False):
            # use automatic limits from plotting
            temp1 = plt1.tricontourf(ResultRef.getVar("x"), ResultRef.getVar("y"), ResultRef.getVar(Variable[i]), cmap=cm.rainbow)
            temp2 = plt2.tricontourf(ResultSim.getVar("x"), ResultSim.getVar("y"), ResultSim.getVar(Variable[i]), cmap=cm.rainbow)
        else:
            # use imposed number of contour lines, with z limits forced by reference data
            temp1 = plt1.tricontourf(ResultRef.getVar("x"), ResultRef.getVar("y"), ResultRef.getVar(Variable[i]), np.arange(zmin, zmax, (zmax - zmin) / contour_number), cmap=cm.rainbow)
            temp2 = plt2.tricontourf(ResultSim.getVar("x"), ResultSim.getVar("y"), ResultSim.getVar(Variable[i]), np.arange(zmin, zmax, (zmax - zmin) / contour_number), cmap=cm.rainbow)
        # Set axes labels
        plt1.set(xlabel=plotCommon.sw2d_Label("x", Language))
        plt1.set(ylabel=plotCommon.sw2d_Label("y", Language))
 
        plt2.set(xlabel=plotCommon.sw2d_Label("x", Language))
        plt2.set(ylabel=plotCommon.sw2d_Label("y", Language))

        # Set colorbars
        cbar1 = fig.colorbar(temp1, ax=plt1, format='%.3f')
        cbar1.ax.set_ylabel(plotCommon.sw2d_Legend('ref', Language) + "\n" + plotCommon.sw2d_Label(Variable[i], Language))
        cbar2 = fig.colorbar(temp2, ax=plt2, format='%.3f')
        cbar2.ax.set_ylabel(plotCommon.sw2d_Legend('sim', Language) + "\n" + plotCommon.sw2d_Label(Variable[i], Language))

        # Save in png file
        if (len(file) > 0):
            plt.savefig(file)

def plot_2d_compare_temp(simu, reference):
    """
    Plotting a subplot with temp(x,y) from the <simu> and <reference> files.
    """
    # Clear figure
    plt.clf()

    x1 = []
    y1 = []
    t1 = []
    x2 = []
    y2 = []
    t2 = []
    # Read the files
    with open(simu, 'r') as csvfile:
        # The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        # Remove dummy first row in <filename>
        next(rows)
        # Fill the x and y arrays
        for row in rows:
            # Remove empty strings
            row = [string for string in row if string != ""]
            # Set values from file rows
            x1.append(float(row[0]))
            y1.append(float(row[1]))
            t1.append(float(row[2]))
    with open(reference, 'r') as csvfile:
        # The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        # Remove dummy first row in <filename>
        next(rows)
        # Fill the x and y arrays
        for row in rows:
            # Remove empty strings
            row = [string for string in row if string != ""]
            # Set values from file rows
            x2.append(float(row[0]))
            y2.append(float(row[1]))
            t2.append(float(row[2]))

    # Manage suplotting
    fig, (plt1, plt2) = plt.subplots(2)
    plt.subplots_adjust(hspace = 0.5)
    # Plots
    temp1 = plt1.tricontourf(x1, y1, t1)
    temp2 = plt2.tricontourf(x2, y2, t2)
    # Set axes labels
    plt1.set(ylabel='y(m)')
    plt2.set(xlabel='x(m)')
    plt2.set(ylabel='y(m)')
    # Set colorbars
    fig.colorbar(temp1, ax=plt1, format='%.3f')
    fig.colorbar(temp2, ax=plt2, format='%.3f')
    # Save figure
    plt.savefig('figure_python_results.png')
