import numpy as np

class opening():
    def __init__(self, nom):
        self.Name = nom
        self.Q = []
        self.hStreet = []
        self.hBuilding = []
        
    def addValue(self, Q, hStreet, hBuilding):
        self.Q.append(Q)
        self.hStreet.append(hStreet)
        self.hBuilding.append(hBuilding)
    
class openings():
    def __init__(self,**kwargs):
        self.nb_openings = 0
        self.time = []
        self.openings = []
        
    def checkName(self, nom):
        exist = False        
        for idx, p in enumerate(self.openings):
            if p.Name == nom:
                exist = True
                return idx
        if (exist == False):
            print('opening "' + nom + '" does not exists')
            liste = ""
            for n in self.getNames():
                liste = liste + n + "; "
            print('existing names: ' + liste[:len(liste)-2])
            print("STOPPING\n")
            raise ValueError("Stopping traceback")
        
    def addOpening(self, p):
        self.openings.append(p)
        self.nb_openings += 1
        
    def setTime(self, time_value):
        self.time = time_value
        
    def addOpeningValues(self, i, Q, hStreet, hBuilding):
        self.openings[i].addValue(Q, hStreet, hBuilding)
        
    def getNames(self):
        names = []
        for p in self.openings:
            names.append(p.Name)
            
        return names
        
    def getTime(self):
        return self.time
        
    def getVar(self, probe, var):
        method_name = 'get_' + var
        method = getattr(self, method_name)
        return method(probe)
        
    def get_Q(self, opening):
        return self.openings[opening].Q
            
    def get_hStreet(self, opening):
        return self.openings[opening].hStreet

    def get_hBuilding(self, opening):
        return self.openings[opening].hBuilding
