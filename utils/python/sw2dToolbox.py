import math
import os
import sys

sw2d_folder = os.getcwd().split("test")[0]

sys.path.insert(1, os.path.join(sw2d_folder,"utils/python/classes"))
import classBioMap
import classHydroBC
import classHydroMap
import classHydroProbes
import classOpenings
sys.path.insert(1, os.path.join(sw2d_folder,"utils/python/plot"))
import plot1D
import plot2D
import plotLogJam
sys.path.insert(1, os.path.join(sw2d_folder,"utils/python/mr-tests"))
import read

sep = "*****************************"

def deleteFile(name):
    if os.path.exists(name):
        os.remove(name)

def getFolderList(folder):
    # get list folders in a specified path
    subfolder_list = [f.path[len(folder):] for f in os.scandir(folder) if f.is_dir()]
    subfolder_list.sort()

    ## remove possible slashs
    for idx, name in enumerate(subfolder_list):
        # ending slash
        if ((name[-1] == "/") or
            (name[-1] == "\\")):
            name = name[:-1]
            subfolder_list[idx] = name
        # starting slash
        if ((name[0] == "/") or
            (name[0] == "\\")):
            name = name[1:]
            subfolder_list[idx] = name

    return subfolder_list

def getListPng(Folder):
    if os.path.exists(Folder):
        filelist = os.listdir(Folder)
    else:
        filelist = []

    for fichier in filelist[:]: # filelist[:] makes a copy of filelist.
        if not(fichier.endswith(".png")):
            filelist.remove(fichier)

    return filelist

def getRefFilesList(Folder):
    if os.path.exists(Folder):
        filelist = os.listdir(Folder)
    else:
        filelist = []

    return filelist

def getTestsFolderList(folder, skipped_list):
    folder_list = []
    print(folder)
    #exit()
    if (os.name == 'nt'):
        # windows system
        folder_split = folder.split("\\")
    else:
        folder_split = folder.split("/")
    
    if not (skipped_list.__contains__(folder_split[-1])):
        this_folder_list = getFolderList(folder)
    
        if (this_folder_list.__contains__("input") == True):
            folder_list.append(folder)
        else:
            for current_folder in this_folder_list:
                sub_folder_list = getTestsFolderList(os.path.join(folder, current_folder), skipped_list)
                folder_list += sub_folder_list
                
    return folder_list        

def reference(name):
    return name.replace("_", "-").replace("/", "-").replace("\\", "-")

def repeatString(String, Nb):
    return (String * math.floor(Nb))


def timeToISO8601(Time): # Transforms a floatting number of seconds into the 'p????d??h??m??s???' format 
    iso_date = "pDDDDdHHhMMmSSsMIL"
    days = math.floor(Time / 86400)
    if days < 10:
        days_string = "000"+str(days)
    elif days<100:
        days_string = "00"+str(days)
    elif days<1000:
        days_string = "0"+str(days)
    elif days<10000:
        days_string = str(days)
    else:
        print("Error, date too large")
    hours = math.floor( (Time - days * 86400) / 3600)
    if hours < 10:
        hours_string = "0"+str(hours)
    else:
        hours_string = str(hours)
    minutes = math.floor( (Time - days * 86400 - hours * 3600) / 60)
    if minutes < 10:
        minutes_string = "0"+str(minutes)
    else:
        minutes_string = str(minutes)
    seconds = math.floor( (Time - days * 86400 - hours * 3600 - minutes * 60) )
    if seconds < 10:
        seconds_string = "0"+str(seconds)
    else:
        seconds_string = str(seconds)
    milliseconds = round( (Time - days * 86400 - hours * 3600 - minutes * 60 - seconds) * 1000)
    if milliseconds < 10:
        milliseconds_string = "00"+str(milliseconds)
    elif milliseconds < 100:
        milliseconds_string = "0" + str(milliseconds)
    else:
        milliseconds_string = str(milliseconds)

    iso_date = iso_date.replace("DDDD", days_string)
    iso_date = iso_date.replace("HH", hours_string)
    iso_date = iso_date.replace("MM", minutes_string)
    iso_date = iso_date.replace("SS", seconds_string)
    iso_date = iso_date.replace("MIL", milliseconds_string)
    return iso_date