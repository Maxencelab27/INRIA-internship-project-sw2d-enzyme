import numpy as np

class hydroBC():
    def __init__(self, nom):
        self.Name = nom
        self.Q = []
        
    def addValue(self, Q):
        self.Q.append(Q)
    
class hydroBCs():
    def __init__(self,**kwargs):
        self.nb_bc = 0
        self.time = []
        self.BCs = []
        
    def checkName(self, nom):
        exist = False        
        for idx, p in enumerate(self.BCs):
            if p.Name == nom:
                exist = True
                return idx
        if (exist == False):
            print('bc "' + nom + '" does not exists')
            print('stop')
            exit()
            
    def addBC(self, p):
        self.BCs.append(p)
        self.nb_bc += 1
        
    def setTime(self, time_value):
        self.time = time_value
        
    def addBCValues(self, i, Q):
        self.BCs[i].addValue(Q)
        
    def getTime(self):
        return self.time
        
    def getVar(self, i, var):
        method_name = 'get_' + var
        method = getattr(self, method_name)
        return method(i)
        
    def get_Q(self, i):
        return self.BCs[i].Q
