import numpy as np

class variable():
    def __init__(self, nbCell, name):
        self.name = name
        self.val = [np.nan] * nbCell
        
    def setVar(self, value, *position):
        if len(position) == 0:
            self.val = value
        else:
            self.val[position[0]] = value
        return
    
    def getVar(self):
        return self.val

class hydroMap():
    def __init__(self, **kwargs):
        # If no size provided, init with 1 cell
        self.nbCells = 1
        self.names = []
        self.var = []
        # If size provided thanks to nbCells keyword, use it
        for key, value in kwargs.items():
            if (key == "nbCells"):
                self.nbCells = value
            # Otherwise, abort
            else:
                print('Bad keyword argument when creating ' + self + '. Aborting...')
                quit()
                
    def setVar(self, name, value, *position):
        if not(name in self.names):
            # create new data
            self.names.append(name)
            self.var.append(variable(self.nbCells, name))
            
        if len(position) == 0:
            self.var[self.names.index(name)].setVar(value)
        else:
            self.var[self.names.index(name)].setVar(value, position[0])

        return
    
    def getVar(self, name):
        if not(name in self.names):
            print("variable", name, "is not available")
            exit()
            
        return self.var[self.names.index(name)].getVar()

    #def getVar(self, var):
        #method_name = 'get_' + var
        #method = getattr(self, method_name)
        #return method()

    #def setVar(self, var, value, *position):
        #if var in ['x', 'y', 'z', 'zb', 'h', 'qx', 'qy']:
            #method_name = 'set_' + var
            #method = getattr(self, method_name)
            #return method(value, *position)
        #else:
            #print('Trying to set unknown field <' + var + '>. Abort...')
            #quit()

    #def get_x(self):
        #return self.x

    #def set_x(self, value, *position):
        #if len(position) == 0:
            #self.x = value
        #else:
            #self.x[position[0]] = value
        #return

    #def get_y(self):
        #return self.y

    #def set_y(self, value, *position):
        #if len(position) == 0:
            #self.y = value
        #else:
            #self.y[position[0]] = value
        #return
    
    #def get_z(self):
        #self.update_z()
        #return self.z
    
    #def set_z(self, value, *position):
        #if len(position) == 0:
            #self.z = value
        #else:
            #self.z[position[0]] = value

        #return    
    
    #def update_z(self):
        #self.z = [np.nan] * self.nbCells
        #for i, h in enumerate(self.h):
            #self.z[i] = self.zb[i] + h

        #return

    #def get_zb(self):
        #return self.zb

    #def set_zb(self, value, *position):
        #if len(position) == 0:
            #self.zb = value
        #else:
            #self.zb[position[0]] = value
            
        #return

    #def get_h(self):
        #return self.h

    #def set_h(self, value, *position):
        #if len(position) == 0:
            #self.h = value
        #else:
            #self.h[position[0]] = value
        
        #return

    def get_nbCells(self):
        return self.nbCells

    def set_nbCells(self, nb):
        self.nbCells = nb
        return

    #def get_qx(self):
        #return self.qx

    #def set_qx(self, value, *position):
        #if len(position) == 0:
            #self.qx = value
        #else:
            #self.qx[position[0]] = value
        #return

    #def get_qy(self):
        #return self.qy

    #def set_qy(self, value, *position):
        #if len(position) == 0:
            #self.qy = value
        #else:
            #self.qy[position[0]] = value
        #return

    #def get_u(self):
        #return self.u

    #def set_u(self, value, *position):
        #if len(position) == 0:
            #self.u = value
        #else:
            #self.u[position[0]] = value
        #return

    #def get_v(self):
        #return self.v

    #def set_v(self, value, *position):
        #if len(position) == 0:
            #self.v = value
        #else:
            #self.v[position[0]] = value
        #return
