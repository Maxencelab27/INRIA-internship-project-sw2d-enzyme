import numpy as np

class hydroProbe():
    def __init__(self, nom, x, y):
        self.Name = nom
        self.x = x
        self.y = y
        self.zb = []
        self.h = []
        self.qx = []
        self.qy = []
        
    def addValue(self, zb, h, qx, qy):
        self.zb.append(zb)
        self.h.append(h)
        self.qx.append(qx)
        self.qy.append(qy)
    
class hydroProbes():
    def __init__(self,**kwargs):
        self.nb_probes = 0
        self.time = []
        self.probes = []
        
    def checkName(self, nom):
        exist = False        
        for idx, p in enumerate(self.probes):
            if p.Name == nom:
                exist = True
                return idx
        if (exist == False):
            raise ValueError('probe "' + nom + '" does not exists' + '\n' + 'stop')
            
        
    def addProbe(self, p):
        self.probes.append(p)
        self.nb_probes += 1
        
    def setTime(self, time_value):
        self.time = time_value
        
    def addProbeValues(self, i, zb, h, qx, qy):
        self.probes[i].addValue(zb, h, qx, qy)
        
    def getTime(self):
        return self.time
        
    def getVar(self, probe, var):
        method_name = 'get_' + var
        method = getattr(self, method_name)
        return method(probe)
        
    def get_h(self, probe):
        return self.probes[probe].h
            
    def get_qx(self, probe):
        return self.probes[probe].qx

    def get_qy(self, probe):
        return self.probes[probe].qy
    
    def get_zb(self, probe):
        return self.probes[probe].zb
