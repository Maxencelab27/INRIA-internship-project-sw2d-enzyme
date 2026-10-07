def sw2d_Legend(Var, Lang):
    Legend_EN = ["Reference", "Simulation", "Bottom"]

    index = -1
    if (Var == "ref"):
        index = 0
    elif (Var == "sim"):
        index = 1
    elif (Var == "zb"):
        index = 2
    else:
        print("sw2d_Legend - unknown variable:", Var)

    if (Lang == "EN"):
        return Legend_EN[index]
    # elif (Lang == "FR"):
    #    return Legend_FR[index]
    else:
        print("unknown language:", Lang)


def sw2d_Label(Var, Lang):
    Label_EN = ["Water depth (m)", "Unit-discharge (m^2/s)", "x (m)", "y (m)", "heat (°C)",
                "Time (s)", "Surface elevation (m)", "Transverse\n unit-discharge (m^2/s)",
                "x velocity (m/s)", "y velocity (m/s)", "velocity norm (m/s)",
                "Discharge (m3/s)", "Street water depth (m)", "Building water depth (m)", 
                "theta (-)", "x theta velocity (m/s)", "y theta velocity (m/s)",
                "Froude (-)", "Bottom elevation (m)"]

    index = -1
    if (Var == "h"):
        index = 0
    elif (Var == "qx"):
        index = 1
    elif (Var == "x"):
        index = 2
    elif (Var == "y"):
        index = 3
    elif (Var == "heat"):
        index = 4
    elif (Var == "t"):
        index = 5
    elif (Var == "z"):
        index = 6
    elif (Var == "qy"):
        index = 7
    elif (Var == "ux"):
        index = 8
    elif (Var == "uy"):
        index = 9
    elif (Var == "u"):
        index = 10
    elif (Var == "Q"):
        index = 11
    elif (Var == "hStreet"):
        index = 12
    elif (Var == "hBuilding"):
        index = 13
    elif (Var == "theta"):
        index = 14
    elif (Var == "theta_x"):
        index = 15
    elif (Var == "theta_y"):
        index = 16
    elif (Var == "fr"):
        index = 17
    elif (Var == "zb"):
        index = 18
    else:
        print("sw2d_Label - unknown variable:", Var)

    if (Lang == "EN"):
        return Label_EN[index]
    # elif (Lang == "FR"):
    #    return Label_FR[index]
    else:
        print("unknown language:", Lang)
