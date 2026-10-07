import matplotlib.pyplot as plt
import plotCommon


def profile(**kwargs):
    # declare default value
    file = ""
    Language = "EN"
    DefinedSim = False
    DefinedRef = False
    Variable = ["h"]
    PlotZb = False
    Direction = "x"
    SetXmin = False
    Xmin = 0
    SetXmax = False
    Xmax = 0
    SetYmin = False
    Ymin = 0
    SetYmax = False
    Ymax = 0

    # read argument
    for key, value in kwargs.items():
        if (key == "direction"):
            # Define the plotting direction (x or y)
            if (value == "x"):
                Direction = "x"
            elif (value == "y"):
                Direction = "y"
            else:
                print("unallowed parameter value for key:", key, "\n",
                      "allowed value: 'x' or 'y'")
        elif (key == "filename"):
            # Define the output filename
            file = value
        elif (key == "language"):
            # Output language
            Language = value
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
        elif (key == "plot_zb"):
            # Boolean to know if we plot zb on z
            PlotZb = value
        elif (key == "xmax"):
            # Maximal value for X axis
            Xmax = value
            SetXmax = True
        elif (key == "xmin"):
            # Minimal value for X axis
            Xmin = value
            SetXmin = True
        elif (key == "ymax"):
            # Minimal value for X axis
            Ymax = value
            SetYmax = True
        elif (key == "ymin"):
            # Minimal value for X axis
            Ymin = value
            SetYmin = True
        else:
            print("unknown argument:", key)
            exit()

    # check compulsory variables
    if ((DefinedSim == False) and (DefinedRef == False)):
        print("Data have to be provided using 'ref' or 'sim'")
        exit()

    # Clear figure
    plt.clf()

    # plot
    NbPlot = len(Variable)

    for i in range(NbPlot):
        # create the subplot
        if (i == 0):
            plt.subplots(NbPlot, 1, figsize=(6, 2*NbPlot), sharex=True)
            plt.subplots_adjust(hspace = 0.5)
            CurrentPlot = plt.subplot(NbPlot, 1, i + 1)
        else:
            CurrentPlot = plt.subplot(NbPlot, 1, i + 1)

        if (DefinedRef == True):
            if ((Direction == "x") | (Direction == "y")):
                x = ResultRef.getVar(Direction)
                TitleAxis = Direction
            else:
                print("not handled")
                exit()
            var = ResultRef.getVar(Variable[i])
            CurrentPlot.plot(x, var, label=plotCommon.sw2d_Legend('ref', Language), linestyle='solid')
            CurrentPlot.set(ylabel=plotCommon.sw2d_Label(Variable[i], Language))
            CurrentPlot.set(xlabel=plotCommon.sw2d_Label(TitleAxis, Language))

        if (DefinedSim == True):
            if ((Direction == "x") | (Direction == "y")):
                x = ResultSim.getVar(Direction)
                TitleAxis = Direction
            else:
                print("not handled")
                exit()
            var = ResultSim.getVar(Variable[i])
            CurrentPlot.plot(x, var, label=plotCommon.sw2d_Legend('sim', Language), linestyle='none', marker='+')
            yLabel = plotCommon.sw2d_Label(Variable[i], Language)
            if ((Variable[i] == 'z') and (PlotZb == True)):
                var_2 = ResultSim.getVar('zb')
                CurrentPlot.plot(x, var_2, label=plotCommon.sw2d_Legend('zb', Language), linestyle='solid', color='black', zorder=0, lw=0.8)
                yLabel = yLabel + "\n" + plotCommon.sw2d_Label('zb', Language)
                
            CurrentPlot.set(ylabel = yLabel)
            CurrentPlot.set(xlabel = plotCommon.sw2d_Label(TitleAxis, Language))

        CurrentPlot.legend()

    # Defines the axes limits
    if (SetXmax == True): plt.xlim(xmax=Xmax)
    if (SetXmin == True): plt.xlim(xmin=Xmin)
    if (SetYmax == True): plt.ylim(ymax=Ymax)
    if (SetYmin == True): plt.ylim(ymin=Ymin)

    # Save in png file
    if (len(file) > 0):
        plt.savefig(file)

def timeseries(**kwargs):
    # declare default value
    addname = False
    file = ""
    Language = "EN"
    DefinedSim = False
    DefinedRef = False
    Variable = ["h"]
    SetTmin = False
    Tmin = 0
    SetTmax = False
    Tmax = 0
    SetYmin = False
    Ymin = 0
    SetYmax = False
    Ymax = 0
    probe = ""
    
    # read argument
    for key, value in kwargs.items():
        if (key == "addname"):
            # Add the probe name
            addname = value
        elif (key == "filename"):
            # Define the output filename
            file = value
        elif (key == "language"):
            # Output language
            Language = value
        elif (key == "probe"):
            probe = value
        elif (key == "ref"):
            # reference value
            ResultRef = value
            DefinedRef = True
        elif (key == "sim"):
            # SW2D computed value
            ResultSim = value
            DefinedSim = True
        elif (key == "t_max"):
            # Maximal value for X axis
            Tmax = value
            SetTmax = True
        elif (key == "tmin"):
            # Minimal value for X axis
            Tmin = value
            SetTmin = True
        elif (key == "variable"):
            # Table of the plotted variables
            Variable = value
        elif (key == "ymax"):
            # Minimal value for X axis
            Ymax = value
            SetYmax = True
        elif (key == "ymin"):
            # Minimal value for X axis
            Ymin = value
            SetYmin = True
        else:
            print("unknown argument:", key)
            exit()
            
    # check compulsory variables
    if ((DefinedSim == False) and (DefinedRef == False)):
        print("Data have to be provided using 'ref' or 'sim'")
        exit()

    # Clear figure
    plt.clf()
   
    # plot
    NbPlot = len(Variable)
    NbProbes = len(probe)
    if (NbProbes == 0):
        if (DefinedSim == True):
            probe = ResultSim.getNames()
            NbProbes = len(probe)
        else:
            probe = ResultRef.getNames()
        NbProbes = len(probe)

    for i in range(NbPlot):
        # create the subplot
        for j in range(NbProbes):
            if ((i == 0) & (j == 0)):
                plt.subplots(NbPlot, NbProbes, figsize=(6, 2*NbPlot))
                plt.subplots_adjust(hspace = 0.5)
                CurrentPlot = plt.subplot(NbPlot, NbProbes, 1 + i * NbProbes + j)
            else:
                CurrentPlot = plt.subplot(NbPlot, NbProbes, 1 + i * NbProbes + j)
        
            if (DefinedRef == True):
                probe_nb = ResultRef.checkName(probe[j])
                t = ResultRef.getTime()
                TitleAxis = "t"
                var = ResultRef.getVar(probe_nb, Variable[i])
                CurrentPlot.plot(t, var, label = plotCommon.sw2d_Legend('ref', Language), linestyle='solid')
                CurrentPlot.set(ylabel=plotCommon.sw2d_Label(Variable[i], Language))
                CurrentPlot.set(xlabel=plotCommon.sw2d_Label(TitleAxis, Language))
        
            if (DefinedSim == True):
                probe_nb = ResultSim.checkName(probe[j])
                t = ResultSim.getTime()
                TitleAxis = "t"
                var = ResultSim.getVar(probe_nb, Variable[i])
                CurrentPlot.plot(t, var, label = plotCommon.sw2d_Legend('sim', Language), linestyle='none', marker='+')
                CurrentPlot.set(ylabel=plotCommon.sw2d_Label(Variable[i], Language))
                CurrentPlot.set(xlabel=plotCommon.sw2d_Label(TitleAxis, Language))
                if ((i == 0) & (addname == True)):
                    plt.title(probe[j])
        
        CurrentPlot.legend()
    
    # Defines the axes limits
    if (SetTmax == True): plt.xlim(xmax = Tmax)
    if (SetTmin == True): plt.xlim(xmin = Tmin)
    if (SetYmax == True): plt.ylim(ymax = Ymax)
    if (SetYmin == True): plt.ylim(ymin = Ymin)
    

    # Save in png file
    if (len(file) > 0):
        plt.savefig(file)
        
def timeseriesQ(**kwargs):
    # declare default value
    file = ""
    Language = "EN"
    DefinedSim = False
    DefinedRef = False
    SetTmin = False
    Tmin = 0
    SetTmax = False
    Tmax = 0
    SetYmin = False
    Ymin = 0
    SetYmax = False
    Ymax = 0
    bc = ""
    
    # read argument
    for key, value in kwargs.items():
        if (key == "filename"):
            # Define the output filename
            file = value
        elif (key == "language"):
            # Output language
            Language = value
        elif (key == "bc"):
            bc = value
        elif (key == "ref"):
            # reference value
            ResultRef = value
            DefinedRef = True
        elif (key == "sim"):
            # SW2D computed value
            ResultSim = value
            DefinedSim = True
        elif (key == "t_max"):
            # Maximal value for X axis
            Tmax = value
            SetTmax = True
        elif (key == "tmin"):
            # Minimal value for X axis
            Tmin = value
            SetTmin = True
        elif (key == "ymax"):
            # Minimal value for Y axis
            Ymax = value
            SetYmax = True
        elif (key == "ymin"):
            # Minimal value for Y axis
            Ymin = value
            SetYmin = True
        else:
            print("unknown argument:", key)
            exit()
            
    # check compulsory variables
    if ((DefinedSim == False) and (DefinedRef == False)):
        print("Data have to be provided using 'ref' or 'sim'")
        exit()

    # Clear figure
    plt.clf()
   
    # plot
    NbBC = len(bc)

    for i in range(NbBC):
        # create the subplot
        #for j in range(NbProbes):
        if (i == 0):
            CurrentPlot = plt.subplot(1, 1, 1)
            #else:
                #CurrentPlot = plt.subplot(NbPlot, NbProbes, 1 + i * NbProbes + j, sharex = CurrentPlot)
        
        if (DefinedRef == True):
            probe_nb = ResultRef.checkName(bc[i])
            t = ResultRef.getTime()
            TitleAxis = "t"
            var = ResultRef.getVar(probe_nb, 'Q')
            CurrentPlot.plot(t, var, label = plotCommon.sw2d_Legend('ref', Language) + " " + bc[i], linestyle='solid')
            CurrentPlot.set(ylabel=plotCommon.sw2d_Label('Q', Language))
            CurrentPlot.set(xlabel=plotCommon.sw2d_Label(TitleAxis, Language))
    
        if (DefinedSim == True):
            probe_nb = ResultSim.checkName(bc[i])
            t = ResultSim.getTime()
            TitleAxis = "t"
            var = ResultSim.getVar(probe_nb, 'Q')
            CurrentPlot.plot(t, var, label = plotCommon.sw2d_Legend('sim', Language) + " " + bc[i], linestyle='none', marker='+')
            CurrentPlot.set(ylabel=plotCommon.sw2d_Label('Q', Language))
            CurrentPlot.set(xlabel=plotCommon.sw2d_Label(TitleAxis, Language))
        
        CurrentPlot.legend()
    
    # Defines the axes limits
    if (SetTmax == True): plt.xlim(xmax = Tmax)
    if (SetTmin == True): plt.xlim(xmin = Tmin)
    if (SetYmax == True): plt.ylim(ymax = Ymax)
    if (SetYmin == True): plt.ylim(ymin = Ymin)
    

    # Save in png file
    if (len(file) > 0):
        plt.savefig(file)