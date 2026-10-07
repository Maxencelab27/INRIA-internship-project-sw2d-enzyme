import csv
import os

from classHydroProbes import hydroProbe
from classOpenings import opening
from classHydroBC import hydroBC

header_table = {
    'fr': 'froude',
    'h': 'water_depth', 'hs': 'specific_head', 'ht': 'total_head',
    'q': 'unit_discharge_norm', 'qx': 'unit_discharge_x', 'qy': 'unit_discharge_y', 
    'theta': 'theta', 'theta_x': 'theta_velocity_x', 'theta_y': 'theta_velocity_y',
    'u': 'velocity_norm', 'ux': 'velocity_x', 'uy': 'velocity_y',
    'x': 'x',
    'y': 'y', 
    'z': 'surface_elevation', 'zb': 'zb', 
}

def allowedHeader(Var):
    if not(Var in header_table):
        print("unknown header for output maps:", Var)
        keys = ''
        for i, key in enumerate(header_table.keys()):
            if (i == len(header_table.keys()) - 1):
                keys += key
            else:
                keys += key + ", "
        
        print("\tavailable Var:", keys)
        return False
    else:
        return True      
    
def getHeader(Var):
    if allowedHeader(Var):
        return header_table[Var]
    else:
        exit()

def mapFileHydro(result, csvFileName, **kwargs):
    """
    Reading a formatted hydro CSV file <csvFilename> into a <result> mapHydroclass object. Optional arguments may indicate specific columns to be read. Example table=['x','h'].
    """
    # Set dictionary
    dictionary = ['x', 'y', 'zb']
    # Set the columns to be read
    if len(kwargs) == 0:  # read all 6 default fields
        print('Use default dictionary [' + ','.join(dictionary) + '] to fill mapHydro class from ' + csvFileName)
        table = dictionary
    else:  # read only required fields
        for key, value in kwargs.items():
            if (key == "table"):  # read only required fields
                table = dictionary
                for var in value:
                    if not(allowedHeader(var)):
                        exit()
                    table.append(var)
                print('Use provided dictionary [' + ','.join(table) + '] to fill mapHydro class from ' + csvFileName)
            else:  # Error
                print('Error reading columns to be read. Aborting... ')
                quit()

    # Count number of records in the file
    with open(csvFileName, 'r') as csvfile:
        # The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        # get column header
        csv_headings = next(rows)
        for i in range(len(csv_headings)) :
            csv_headings[i] = csv_headings[i].replace(" ", "")
        # Remove hashtag and spaces from the begining of the line (if any)
        csv_headings[0] = csv_headings[0].replace("#","")
        csv_headings[0] = csv_headings[0].replace(" ","")
        # Number of records in the sim file
        nbRows = sum(1 for row in rows)
        result.set_nbCells(nbRows)
        
        # check result map contains the desired variables
        for var in table:
            if not(getHeader(var) in csv_headings):
                print("looking for", getHeader(var), "in",csvFileName)
                print("\tthis variable is not stored")
                print("\tstored variables:", csv_headings)
                exit()

    # Open file and fill values
    with open(csvFileName, 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        # Remove dummy first row in <filename>
        next(rows)
        rowNumber = 0
        # Parse file lines
        for row in rows:
            # Parse file columns
            for var in table:
                # print(field)
                col = csv_headings.index(getHeader(var))
                result.setVar(var, float(row[col]), rowNumber)
            rowNumber += 1


def mapFileBio(result, csvFileName, *args):
    """
    Reading a formatted CSV biology file <csvFilename> into a <Result> mapBio class object. Optional arguments may indicate specific columns to be read. Example table=['x','y','heat'].
    """
    # Set dictionary
    dictionary = ['x', 'y', 'heat']

    # Set the columns to be read
    if len(args) == 0:  # read all 6 default fields
        print('Use default dictionary [' + ','.join(dictionary) + '] to fill mapBio class from ' + csvFileName)
        table = dictionary
    else:  # read only required fields
        table = args
        print(args)
        print('Use provided dictionary [' + ','.join(table) + '] to fill mapBio class from ' + csvFileName)

    # Count number of records in the file
    with open(csvFileName, 'r') as csvfile:
        # The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        # Number of records in the sim file
        nbRows = sum(1 for row in rows) - 1
        result.set_nbCells(nbRows)

    # Fields initialization
    for field in table:
        # print('Initialize '+field)
        result.setVar(field, [0] * nbRows)

    # Open file and fill values
    with open(csvFileName, 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        # Remove dummy first row in <filename>
        next(rows)
        rowNumber = 0
        # Parse file lines
        for row in rows:
            # Parse file columns
            for field in table:
                # print(field)
                col = dictionary.index(field)
                result.setVar(field, float(row[col]), rowNumber)
            rowNumber += 1

def BcFileHydro(result, BcFileName, **kwargs):
    ## Count number of records in the file
    if not(os.path.exists(BcFileName)):
        print('inexisting file:', BcFileName)
        exit()
        
    with open(BcFileName, 'r') as csvfile:
        ## The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        
        temps = []

        csv_headings = next(rows)
        for i in range(len(csv_headings)) :
            csv_headings[i] = csv_headings[i].replace(" ", "")
        NbBC = len(csv_headings) - 1
        for i in range(NbBC):
            bc = hydroBC(csv_headings[i + 1])
            result.addBC(bc)
        
        for row in rows:
            i += 1
            temps.append(float(row[0]))
            for i in range(NbBC):
                Q = float(row[1 + i])
                result.addBCValues(i, Q)
                            
    result.setTime(temps)    
    

def orificeFile(result, probesFileName, **kwargs):
    ## Count number of records in the file
    if not(os.path.exists(probesFileName)):
        print('inexisting file:', probesFileName)
        exit()
        
    with open(probesFileName, 'r') as csvfile:
        rows = csv.reader(csvfile)

        first = True
        NbOpenings = 0
        temps = []

        for row in rows:
            if (len(row) != 0):
                # skip empty line
                line = row[0].split("\t")
                
                if (first == True):
                    for i in range(int((len(line) - 1) / 3)):
                        name = line[1 + i * 3]
                        name = name[2:]
                        # add probes
                        p = opening(name)
                        result.addOpening(p)
                        NbOpenings = NbOpenings + 1
                        
                    first = False
                
                else:            
                    temps.append(float(line[0]))
                    for pp in range(NbOpenings):
                        Q = float(line[1 + pp * 3])
                        hStreet = float(line[2 + pp * 3])
                        hBuilding = float(line[3 + pp * 3])
                        result.addOpeningValues(pp, Q, hStreet, hBuilding)
                                    
    result.setTime(temps)                            

def probesFileHydro(result, probesFileName, **kwargs):
    ## Count number of records in the file
    if not(os.path.exists(probesFileName)):
        print('inexisting file:', probesFileName)
        exit()
        
    with open(probesFileName, 'r') as csvfile:
        rows = csv.reader(csvfile)

        i = 0
        NbProbes = 0
        temps = []

        for row in rows:
            i += 1
            if (len(row) != 0):
                # skip empty line
                line = row[0].split()
            
                if (len(line) == 3):
                    # add probes
                    line = row[0].split()
                    p = hydroProbe(line[0], line[1], line[2])
                    result.addProbe(p)
                    NbProbes = NbProbes + 1
                elif (len(line) == 1 + NbProbes * 4):
                    if (line[0] != 'time'):
                        temps.append(float(line[0]))
                        for pp in range(NbProbes):
                            zb = float(line[1 + pp * 4])
                            h = float(line[2 + pp * 4])
                            qx = float(line[3 + pp * 4])
                            qy = float(line[4 + pp * 4])
                            result.addProbeValues(pp, zb, h, qx, qy)  

    result.setTime(temps)
