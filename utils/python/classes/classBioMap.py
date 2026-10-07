import numpy as np

class bioMap():
    def __init__(self,**kwargs):
        # If no size provided, init with 1 cell
        self.nbCells = 1
        # If size provided thanks to nbCells keyword, use it
        for key, value in kwargs.items():
            if (key == "nbCells"):
                self.nbCells = value
            # Otherwise, abort
            else:
                print('Bad keyword argument when creating ' + self + '. Aborting...')
                quit()
        # Fill all fields with nan values
        for var in ['x', 'y', 'heat']:
            self.setVar(var, [np.nan] * self.nbCells)

    def getVar(self, var):
        method_name = 'get_' + var
        method = getattr(self, method_name)
        return method()

    def setVar(self, var, value, *position):
        if var in ['x', 'y', 'heat']:
            method_name = 'set_' + var
            method = getattr(self, method_name)
            return method(value, *position)
        else:
            print('Trying to set unknown field <' + var + '>. Abort...')
            quit()

    def get_nbCells(self):
        return self.nbCells

    def set_nbCells(self, nb):
        self.nbCells = nb
        return

    def get_x(self):
        return self.x

    def set_x(self, value, *position):
        if len(position) == 0:
            self.x = value
        else:
            self.x[position[0]] = value
        return

    def get_y(self):
        return self.y

    def set_y(self, value, *position):
        if len(position) == 0:
            self.y = value
        else:
            self.y[position[0]] = value
        return

    def get_heat(self):
        return self.heat

    def set_heat(self, value, *position):
        if len(position) == 0:
            self.heat = value
        else:
            self.heat[position[0]] = value
        return
