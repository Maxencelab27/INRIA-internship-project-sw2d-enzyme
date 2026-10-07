import numpy as np

class Cylinder:   
    def __init__(self, list_arg):
        self.id = list_arg[0]
        self.position = [np.array([list_arg[1], list_arg[2]])]
        self.velocity = [np.array([list_arg[3], list_arg[4]])]
        self.theta = [list_arg[5]]
        self.omega = [list_arg[6]]
        self.rho = list_arg[7]
        self.length = list_arg[8]

    def __str__(self):
        my_str = "## Id : " + str(self.id) + "\n\t Coord : "+ str(self.position) + "\n\t Theta : " + str([180/np.pi*theta for theta in self.theta]) +\
            "\n\tVelocity : " + str(self.velocity[-1]) + "\n\tOmega : " + str(self.omega[-1]) +"\n\t L : " + str(self.length) +"\n"
        return my_str

    def compute_top_bottom(self):
        T = self.position[-1] + self.length/2*np.array([np.cos(self.theta[-1]), np.sin(self.theta[-1])])
        B = self.position[-1] - self.length/2*np.array([np.cos(self.theta[-1]), np.sin(self.theta[-1])])
        return T, B

    def plot_cylinder(self, ax, color='brown'):
        x, y = self.position[-1][0], self.position[-1][1]
        T, B = self.compute_top_bottom()
        ax.plot((T[0], B[0]), (T[1], B[1]), ls='-', color=color)
        ax.plot(x, y, 'x', color='blue')


class Wall:
    def __init__(self, list_arg) -> None:
        self.id = list_arg[0]
        self.position = [np.array([list_arg[1], list_arg[2]])]
        self.theta = [list_arg[3]]
        self.length = list_arg[4]

    def __str__(self):
        my_str = "## Id : " + str(self.id) + "\n\t Coord : "+ str(self.position) + "\n\t Theta : " + str([180/np.pi*theta for theta in self.theta]) +\
            "\n\t L : " + str(self.length) +"\n"
        return my_str

    def compute_top_bottom(self):
        T = self.position[-1] + self.length/2*np.array([np.cos(self.theta[-1]), np.sin(self.theta[-1])])
        B = self.position[-1] - self.length/2*np.array([np.cos(self.theta[-1]), np.sin(self.theta[-1])])
        return T, B

    def plot_wall(self, ax, color='black'):
        x, y = self.position[-1][0], self.position[-1][1]
        T, B = self.compute_top_bottom()
        ax.plot((T[0], B[0]), (T[1], B[1]), ls='-', color=color)
