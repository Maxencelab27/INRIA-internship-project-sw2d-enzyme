import numpy as np
import sys

import readBinResult

def create_text(result, folder) :
    """Create a text files in sw2d format containing the sw2d simulation result
    
    Keyword arguments :
    result -- an object of readBinResult.BinResult class containing the simulation result
    folder -- string containing the name of the simulation folder where the file will be created
    """
    
    for t in range(len(result.getTimesMaps())) :
        time = result.getTimesMaps()[t]
        
        filename = folder + "/output/ResultHydro_" + date_to_ISO8601(time) + ".txt"
        file = open(filename, mode="w")

        # write header
        file.write("x\ty\tzb")
        vars = result.getVarNames()
        for i in range(1, len(vars)) : # don't write bottom elevation
            file.write("\t" + vars[i])
        file.write("\n")

        # write values
        result_hydro = result.getResultHydro()[t]
        nb_var = len(result.getVarNames())
        nb_cells = result.getResultHydro().shape[1]

        for c in range(nb_cells) :
            file.write(f"{result.getCellsX()[c]:.{10}e}"+ "\t" + f"{result.getCellsY()[c]:.{10}e}")
            for v in range(nb_var) :
                file.write("\t" + f"{result_hydro[c,v]:.{10}e}")
            file.write("\n")

        file.close()


def date_to_ISO8601(date) :
    """Transform a date into a string under ISO8601 norm : p****d**h**m**s*** from days to milliseconds
    
    Keyword argument :
    date -- double number specifying the date in seconds
    """
    
    time = int(np.round(date * 1000.)) # convert time into milliseconds

    milliseconds = time % 1000
    time = time // 1000
    seconds = time % 60
    time = time // 60
    minutes = time % 60
    time = time // 60
    hours = time % 24
    time = time // 24
    days = time

    iso_date = "p" + str(days).rjust(4,'0') + "d" \
                   + str(hours).rjust(2,'0') + "h" \
                   + str(minutes).rjust(2,'0') + "m" \
                   + str(seconds).rjust(2,'0') + "s" \
                   + str(milliseconds).ljust(3,'0')
    
    return(iso_date)
