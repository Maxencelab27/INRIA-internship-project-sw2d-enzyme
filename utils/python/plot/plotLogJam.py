import csv
import glob
import numpy as np
from classTransportType import *
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.cm as cm

header_table = {
    'fr': 'froude',
    'h': 'water_depth', 'hs': 'specific_head', 'ht': 'total_head',
    'q': 'unit_discharge_norm', 'qx': 'unit_discharge_x', 'qy': 'unit_discharge_y', 
    'theta': 'theta', 'theta_x': 'theta_velocity_y', 'theta_y': 'theta_velocity_y',
    'u': 'velocity_norm', 'ux': 'velocity_x', 'uy': 'velocity_y',
    'x': 'x',
    'y': 'y', 
    'z': 'surface_elevation', 'zb': 'zb', 
}

def myStrToFloat(my_string):
    """
    Convert a string into a float
    If my_string is M_PI_**, the output is pi/**. Can use M_PI for pi
    """
    if (my_string[0:4]=="M_PI"):
        if(len(my_string)==4):
            return np.pi
        elif(len(my_string)>= 6):
            return np.pi / float(my_string[5:])
        else:
            raise "Error in converting PI"
    else:
        return float(my_string)


def readResultHydroFile(csvFileName, list_header):
    """
    Read the columns which header is in list_header from a formatted ResultHydro CSV file <csvFilename>.
    Store all the rows in a np.array((nbHeader, nbRows)) and return it
    """
    nbHeaders = len(list_header)

    # Count number of records in the file
    with open(csvFileName, 'r') as csvfile:
        # The separator name is expected to be a tab (\t) caracter.
        rows = csv.reader(csvfile, delimiter='\t')
        # get column header
        csv_headings = next(rows)
        for i in range(len(csv_headings)) :
            csv_headings[i] = csv_headings[i].replace(" ", "")
        # Number of records in the sim file
        nbRows = sum(1 for row in rows)

    rslt_table = np.zeros((nbHeaders, nbRows), dtype='float')

    # Open file and fill values
    with open(csvFileName, 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        # Remove dummy first row in <filename>
        csv_headings = next(rows)
        for i in range(len(csv_headings)) :
            csv_headings[i] = csv_headings[i].replace(" ", "")
        rowNumber = 0
        # Parse file lines
        for row in rows:
            # Parse file columns
            for ind, header in enumerate(list_header):
                col = csv_headings.index(header_table[header])
                rslt_table[ind, rowNumber] = float(row[col])
            rowNumber += 1

    return rslt_table

def readResultLogJamFile(csvFileName, dict = {}):
    """
    Read a ResultLogJam_p*.txt file, and fill the dictionary dict with all the Cylinder present in the file.
    Keys are segments id, value are Cylinders
    If an Id is already present inside the dictionary, the function add new positions, velocities to
    the corresponding Cylinder
    """
    with open(csvFileName, 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        first_row = next(rows)
        if (first_row == ["### LogJam ###"]):
            logjam_headings = next(rows)
        for current_row in rows:
            if (current_row == ["### Wall ###"]):
                break
            new_id = int(current_row[0])
            
            if not (new_id in dict):
                param = [myStrToFloat(current_row[k]) for k in range(9)]
                dict[new_id] = Cylinder(param)
            
            else:
                dict[new_id].position.append(np.array([myStrToFloat(current_row[1]), myStrToFloat(current_row[2])]))
                dict[new_id].velocity.append(np.array([myStrToFloat(current_row[3]), myStrToFloat(current_row[4])]))
                dict[new_id].theta.append(myStrToFloat(current_row[5]))
                dict[new_id].omega.append(myStrToFloat(current_row[6]))
            
    return dict

def readWall(csvFileName, dict = {}):
    """
    Read a ResultLogJam_p*.txt file, and fill the dictionary dict with Wall present in the file.
    Keys are segments id, value are Wall
    If an Id is already present inside the dictionary, the function add new positions, velocities to
    the corresponding Wall
    """
    with open(csvFileName, 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        row = next(rows)
        while (row != ["### Wall ###"]):
            row = next(rows)
        
        wall_headings = next(rows)
        for current_row in rows:
            new_id = int(current_row[0])
            
            if not (new_id in dict):
                param = [myStrToFloat(current_row[k]) for k in range(5)]
                dict[new_id] = Wall(param)
            
            else:
                dict[new_id].position.append(np.array([myStrToFloat(current_row[1]), myStrToFloat(current_row[2])]))
                dict[new_id].theta.append(myStrToFloat(current_row[3]))
            
    return dict


def readAllResultLogJamFile(path):
    """
    Read all the files ResultLogJam*.txt and which are present in path.
    Return a dict containing all the Cylinders described inside those files. 
    """

    dict = {}
    name_files = glob.glob(path+"ResultLogJam*.txt")
    name_files.sort()
    for name_file in name_files:
        readResultLogJamFile(name_file, dict)
    
    return dict

def readAllWall(path):
    """
    Read all the files ResultLogJam*.txt and which are present in path.
    Return a dict containing all the Walls described inside those files. 
    """

    dict = {}
    name_files = glob.glob(path+"ResultLogJam*.txt")
    name_files.sort()
    for name_file in name_files:
        readWall(name_file, dict)
    
    return dict

def readAllResultHydroFile(path):
    """
    Read all the files ResultHydro*.txt and which are present in path.
    Return a list containing the hydrodynamic variables 'x', 'y', 'ux', 'uy', 'u' for each cell
    Format : [[cell_i]_time]
    """
    return_list = []
    name_files = glob.glob(path+"ResultHydro*.txt")
    name_files.sort()
    for name_file in name_files:
        return_list.append(readResultHydroFile(name_file, ['x', 'y', 'ux', 'uy', 'u', 'h']))
    
    return return_list

def readAllTime(path):
    """
    Read all the name_files ResultHydro_p_time.txt which are present in path.
    Return a list containing the time, under format ****d**m**h**s***
    Format : [time]
    """
    return_list = []
    length_path = len(path)
    name_files = glob.glob(path+"ResultHydro*.txt")
    ind_start_time = length_path + len("ResultHydro_p")
    name_files.sort()
    for name_file in name_files:
        return_list.append(name_file[ind_start_time : -4])
    
    return return_list

def addTwoPi(x_tab, y_tab):
    phi_tab = []
    compteur = 0
    add = False
    for (x, y) in zip(x_tab, y_tab):
        phi = np.arctan2(y, x)
        if (x < 0) and (y < 0) and (add):
            compteur += 1
            add = False
        if (x > 0) and (y> 0):
            add = True
        phi_tab.append(phi + 2*np.pi*compteur)
    return phi_tab


def plot1DCylindersRefVersusSim(RefCylinder, SimCylinder, variables, filename, **kwargs):
    """
    Take two lists of Cylinder and plot their variables (x, y, v_x, v_y, theta) given by variables versus time.
    """
    number_variables = len(variables)
    fig, list_ax = plt.subplots(number_variables, 1)
    plt.subplots_adjust(hspace = 0.5)
    
    if (number_variables == 1):
        list_ax = [list_ax]

    if "time" in kwargs:
        time = kwargs["time"]
        time_ref = time[0]
        time_sim = time[1]
    else:
        time_ref = range(0,len(RefCylinder.position))
        time_sim = range(0, len(SimCylinder.position))

    for ind in range(number_variables):
        if (variables[ind] == "x"):
            list_ax[ind].plot(time_ref, [pos[0] for pos in RefCylinder.position])
            list_ax[ind].plot(time_sim, [pos[0] for pos in SimCylinder.position], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$x (m)$")
        elif (variables[ind] == "y"):
            list_ax[ind].plot(time_ref, [pos[1] for pos in RefCylinder.position])
            list_ax[ind].plot(time_sim, [pos[1] for pos in SimCylinder.position], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$y (m)$")
        elif (variables[ind] == "v_x"):
            list_ax[ind].plot(time_ref, [vel[0] for vel in RefCylinder.velocity])
            list_ax[ind].plot(time_sim, [vel[0] for vel in SimCylinder.velocity], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$v_x (ms^{-1})$")
        elif (variables[ind] == "v_y"):
            list_ax[ind].plot(time_ref, [vel[1] for vel in RefCylinder.velocity])
            list_ax[ind].plot(time_sim, [vel[1] for vel in SimCylinder.velocity], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"v_y (ms^{-1}")
        elif (variables[ind] == "theta"):
            list_ax[ind].plot(time_ref, [180/np.pi * theta for theta in RefCylinder.theta])
            list_ax[ind].plot(time_sim, [180/np.pi *theta for theta in SimCylinder.theta], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$\theta (deg)$")
        elif (variables[ind] == "r"):
            r_0 = RefCylinder.position[0][0]
            list_ax[ind].plot(time_ref, [pos[0]/r_0 for pos in RefCylinder.position])
            list_ax[ind].plot(time_sim, [round(np.sqrt(pos[0]**2 + pos[1]**2)/r_0,4) for pos in SimCylinder.position], marker ='+', linestyle='none', color='red', markersize=7)
            list_ax[ind].set_ylabel(r"$r_b(t) \ (meter)$")
        elif (variables[ind] == "phi"):
            phi_sim = addTwoPi([pos[0] for pos in SimCylinder.position], [pos[1] for pos in SimCylinder.position])
            list_ax[ind].plot(time_ref, [pos[1] for pos in RefCylinder.position])
            list_ax[ind].plot(time_sim, phi_sim, marker ='+', linestyle='none', color='red', markersize=7)
            list_ax[ind].set_ylabel(r"$\phi_b(t) \ (rad)$")
        list_ax[ind].legend(["Reference", "Simulation"])

        if "xlim" in kwargs:
            list_ax[ind].set_xlim(kwargs["xlim"])
        list_ax[ind].set_xlabel(r'$Time \ (s)$')
    plt.savefig(filename)

def plot1DCylindersMaxVersusSim(DicMax, SimCylinder, variables, filename, **kwargs):
    """
    Take one list of Cylinder and plot their variables (x, y, v_x, v_y, theta) given by variables versus time.
    DicMax indicates both maximum and minimum for each variables. Constant lines equal at those extremum will be plot. 
    """
    number_variables = len(variables)
    fig, list_ax = plt.subplots(number_variables, 1)
    plt.subplots_adjust(hspace = 0.5)
    if (number_variables == 1):
        list_ax = [list_ax]

    number_time_step_sim = len(SimCylinder.position)

    for ind in range(number_variables):
        if (variables[ind] == "x"):
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["x_min"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["x_max"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [pos[0] for pos in SimCylinder.position], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$x (m)$")
        elif (variables[ind] == "y"):
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["y_min"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["y_max"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [pos[1] for pos in SimCylinder.position], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$y (m)$")
        elif (variables[ind] == "v_x"):
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["v_x_min"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["v_x_max"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [vel[0] for vel in SimCylinder.velocity], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$v_x (ms^{-1})$")
        elif (variables[ind] == "v_y"):
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["v_y_min"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["v_y_max"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [vel[1] for vel in SimCylinder.velocity], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"v_y (ms^{-1}")
        elif (variables[ind] == "theta"):
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["theta_min"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [DicMax["theta_max"] for _ in SimCylinder.position])
            list_ax[ind].plot(range(0, number_time_step_sim), [180/np.pi *theta for theta in SimCylinder.theta], marker ='+', linestyle='none')
            list_ax[ind].set_ylabel(r"$\theta (deg)$")
        list_ax[ind].legend(["Minimum", "Maximum", "Simulation"])

        if "xlim" in kwargs:
            list_ax[ind].set_xlim(kwargs["xlim"])
        list_ax[ind].set_xlabel(r'$time \ (s)$')

    plt.savefig(filename)

def plot2DCylinders(list_Cylinder, **kwargs):
    fig, ax = plt.subplots()
    plt.subplots_adjust(hspace = 0.5)
    ax.set_aspect('equal')
    for ind, cyl in enumerate(list_Cylinder):
        x_center = []
        y_center = []
        # for ind_pos in range(0, len(cyl.position)-60, 4):
        #     x = cyl.position[ind_pos] + cyl.length/2*np.array([np.cos(cyl.theta[ind_pos]), np.sin(cyl.theta[ind_pos])])
        #     y = cyl.position[ind_pos] - cyl.length/2*np.array([np.cos(cyl.theta[ind_pos]), np.sin(cyl.theta[ind_pos])])
            # ax.plot([x[0], y[0]], [x[1], y[1]], color = 'brown', linewidth = 5)
            
        for ind_pos in range(0, len(cyl.position)):
            x_center.append(cyl.position[ind_pos][0])
            y_center.append(cyl.position[ind_pos][1])
        
        ax.plot(x_center, y_center, '-', color='red')
    
    if "hydro" in kwargs:
        hydro = kwargs["hydro"]
        ax.quiver(hydro[0], hydro[1], hydro[2], hydro[3], hydro[4], cmap='winter')
    ax.set_xlabel(r'$x (m)$')
    ax.set_ylabel(r'$y (m)$')


def animationCylinder(list_Cylinder, list_Wall, list_time, speed = 1e3, filename=None, xlim = None, ylim = None):
    """
    Make a plt animation using list of Cylinder, Wall and Hydrodinamic velocities fields.
    Each Cylinder contains all the positions it takes. All those position match with a time and with a Hydro fields.
    Ex : list_Cylinder = [Cylinder([x,y]=[(x0,y0), (x0.48, y0.48), ...], [theta=theta0, theta0.48,...]), Cylinder(), ...], list_Hydro = [Hydro_fileds0, HydroFields0.48, ...]
    Speed is the time step between each update, in ms.
    """
    fig, ax = plt.subplots()
    fig.set_size_inches(18.5, 10.5)
    if (xlim != None):
        ax.set_xlim(xlim)
    if (ylim != None):
        ax.set_ylim(ylim)
    ax.set_aspect('equal')
    cylinders = []
    centers = []
    walls = []
    for _ in list_Cylinder:
        cylinders.append(ax.plot([], [], '-', color = 'brown', linewidth = 5)[0])
        centers.append(ax.plot([], [], 'x', color='red')[0])
    for _ in list_Wall:
        walls.append(ax.plot([], [], '-', color='black')[0])
    
    range_time = len(list_Cylinder[0].position)

    def animate(ind_t):
        c_ind_t = ind_t%range_time
        # print(c_ind_t)
        for ind, cyl in enumerate(list_Cylinder):
            x = cyl.position[c_ind_t] + cyl.length/2*np.array([np.cos(cyl.theta[c_ind_t]), np.sin(cyl.theta[c_ind_t])])
            y = cyl.position[c_ind_t] - cyl.length/2*np.array([np.cos(cyl.theta[c_ind_t]), np.sin(cyl.theta[c_ind_t])])
            cylinders[ind].set_data([x[0], y[0]], [x[1], y[1]])
            centers[ind].set_data(cyl.position[c_ind_t][0], cyl.position[c_ind_t][1])
        
        for ind, wall in enumerate(list_Wall):
            x = wall.position[-1] + wall.length/2*np.array([np.cos(wall.theta[-1]), np.sin(wall.theta[-1])])
            y = wall.position[-1] - wall.length/2*np.array([np.cos(wall.theta[-1]), np.sin(wall.theta[-1])])
            walls[ind].set_data([x[0], y[0]], [x[1], y[1]])
        
        ax.lines[-1].set_label("current_time = {}".format(list_time[c_ind_t]))
        ax.legend()
        ax.set_xlabel(r"$x (m)$")
        ax.set_ylabel(r"$y (m)$")
        
        return centers+cylinders+walls

    anim = animation.FuncAnimation(fig, animate, save_count=range_time, interval=speed, blit=False, repeat=False)
    
    fig.tight_layout()
    if (filename != None):
        writergif = animation.PillowWriter(fps=1000)
        anim.save(filename, writergif)

    plt.show()