import struct
import numpy as np

class BinResult() :
    """This class is used to manage all the information contained in the results of an sw2d simulation"""

    def __init__(self) :
        #-- Geometry
        # number of elements
        self.__nb_cells = 0
        self.__nb_edges_bnd = 0
        self.__nb_edges_int = 0
        self.__nb_nodes = 0
        self.__nb_boundaries = 0

        # cells geometry
        self.__l_cell_id = []
        self.__l_cell_sms_id = []
        self.__l_cell_nb_node = []
        self.__l_cell_nb_edge_bnd = []
        self.__l_cell_nb_edge_int = []
        self.__l_cell_nb_neigh = []
        self.__l_cell_area = []
        self.__l_cell_perimeter = []
        self.__l_cell_x_g = []
        self.__l_cell_y_g = []
        self.__l_cell_zBottom = []
        self.__l_cell_l_node_id = []
        self.__l_cell_l_edge_bnd_id = []
        self.__l_cell_l_edge_int_id = []
        self.__l_cell_l_neigh_id = []

        # nodes geometry
        self.__l_node_id = []
        self.__l_node_sms_id = []
        self.__l_node_x = []
        self.__l_node_y = []
        self.__l_node_z = []
        self.__l_node_nb_cell = []
        self.__l_node_nb_edge = []
        self.__l_node_l_cell_id = []

        # boundary edges geometry
        self.__l_edge_bnd_id = []
        self.__l_edge_bnd_node_1_id = []
        self.__l_edge_bnd_node_2_id = []
        self.__l_edge_bnd_x_m = []
        self.__l_edge_bnd_y_m = []
        self.__l_edge_bnd_z_m = []
        self.__l_edge_bnd_lenght = []
        self.__l_edge_bnd_d_norm = []
        self.__l_edge_bnd_norm_x = []
        self.__l_edge_bnd_norm_y = []
        self.__l_edge_bnd_l_cell_id = []

        # internal edges geometry
        self.__l_edge_int_id = []
        self.__l_edge_int_node_1_id = []
        self.__l_edge_int_node_2_id = []
        self.__l_edge_int_x_m = []
        self.__l_edge_int_y_m = []
        self.__l_edge_int_z_m = []
        self.__l_edge_int_lenght = []
        self.__l_edge_int_d_norm = []
        self.__l_edge_int_norm_x = []
        self.__l_edge_int_norm_y = []
        self.__l_edge_int_l_cell_id = []

        # boundaries
        self.__l_boundaries = []
        self.__l_bc_code = []
        #== Geometry

        # Variables
        self.__nb_var = 0
        self.__l_var_names = []

        # Time
        self.__t_0 = 0
        self.__t_max = 0
        self.__dt_map = 0
        self.__times_maps = []

        # Results
        self.__result_hydro = np.array([])

    def getNbCells(self) :
        return self.__nb_cells
    
    def getNbEdges(self) :
        return self.__nb_edges_bnd + self.__nb_edges_int
    
    def getNbNodes(self) :
        return self.__nb_nodes
    
    def getCellsX(self) :
        return self.__l_cell_x_g
    
    def getCellsY(self) :
        return self.__l_cell_y_g
    
    def getCellsNodes(self) :
        return self.__l_cell_l_node_id
    
    def getEdgesX(self) :
        return self.__l_edge_bnd_x_m + self.__l_edge_int_x_m
    
    def getEdgesY(self) :
        return self.__l_edge_bnd_y_m + self.__l_edge_int_y_m
    
    def getEdgesCells(self) :
        return self.__l_edge_bnd_l_cell_id + self.__l_edge_int_l_cell_id
    
    def getEdgesNode1(self) :
        return self.__l_edge_bnd_node_1_id + self.__l_edge_int_node_1_id
    
    def getEdgesNode2(self) :
        return self.__l_edge_bnd_node_2_id + self.__l_edge_int_node_2_id
    
    def getNodesX(self) :
        return self.__l_node_x
    
    def getNodesY(self) :
        return self.__l_node_y
    
    def getNodesZ(self) :
        return self.__l_node_z
    
    def getNbVar(self) :
        return self.__nb_var
    
    def getVarNames(self) :
        return self.__l_var_names
    
    def getResultHydro(self) :
        return self.__result_hydro
    
    def getDtMap(self) :
        return self.__dt_map
    
    def getTimesMaps(self) :
        return self.__times_maps

    def readFile(self, filename) :
        """Assigns values read from an sw2d binary output file to attributes"""

        # open binary file to read
        ifile = open(filename, mode="rb")

        #--- Read geometry
        #-- read version number
        version_lenght = struct.unpack("Q",ifile.read(8))[0]
        version = ifile.read(version_lenght).decode("utf-8")
        #== read version number

        #-- read number of elements
        self.__nb_cells = struct.unpack("Q",ifile.read(8))[0]
        self.__nb_edges_bnd = struct.unpack("Q",ifile.read(8))[0]
        self.__nb_edges_int = struct.unpack("Q",ifile.read(8))[0]
        self.__nb_nodes = struct.unpack("Q",ifile.read(8))[0]
        self.__nb_boundaries = struct.unpack("Q",ifile.read(8))[0]
        #== read number of elements

        #-- read cells
        self.__l_cell_id = []
        self.__l_cell_sms_id = []
        self.__l_cell_nb_node = []
        self.__l_cell_nb_edge_bnd = []
        self.__l_cell_nb_edge_int = []
        self.__l_cell_nb_neigh = []
        self.__l_cell_area = []
        self.__l_cell_perimeter = []
        self.__l_cell_x_g = []
        self.__l_cell_y_g = []
        self.__l_cell_zBottom = []

        for i in range(self.__nb_cells) :
            id = struct.unpack("Q",ifile.read(8))[0]
            sms_id = struct.unpack("Q",ifile.read(8))[0]
            nodes_size = struct.unpack("Q",ifile.read(8))[0]
            edges_bnd_size = struct.unpack("Q",ifile.read(8))[0]
            edges_int_size = struct.unpack("Q",ifile.read(8))[0]
            neigh_size = struct.unpack("Q",ifile.read(8))[0]
            area = struct.unpack("d",ifile.read(8))[0]
            perimeter = struct.unpack("d",ifile.read(8))[0]
            x_g = struct.unpack("d",ifile.read(8))[0]
            y_g = struct.unpack("d",ifile.read(8))[0]
            zBottom = struct.unpack("d",ifile.read(8))[0]

            self.__l_cell_id.append(id)
            self.__l_cell_sms_id.append(sms_id)
            self.__l_cell_nb_node.append(nodes_size)
            self.__l_cell_nb_edge_bnd.append(edges_bnd_size)
            self.__l_cell_nb_edge_int.append(edges_int_size)
            self.__l_cell_nb_neigh.append(neigh_size)
            self.__l_cell_area.append(area)
            self.__l_cell_perimeter.append(perimeter)
            self.__l_cell_x_g.append(x_g)
            self.__l_cell_y_g.append(y_g)
            self.__l_cell_zBottom.append(zBottom)
        #== read cells

        #-- read nodes
        self.__l_node_id = []
        self.__l_node_sms_id = []
        self.__l_node_x = []
        self.__l_node_y = []
        self.__l_node_z = []
        self.__l_node_nb_cell = []
        self.__l_node_nb_edge = []

        for i in range(self.__nb_nodes) :
            id = struct.unpack("Q",ifile.read(8))[0]
            sms_id = struct.unpack("Q",ifile.read(8))[0]
            x = struct.unpack("d",ifile.read(8))[0]
            y = struct.unpack("d",ifile.read(8))[0]
            z = struct.unpack("d",ifile.read(8))[0]
            cells_size = struct.unpack("Q",ifile.read(8))[0]
            edges_size = struct.unpack("Q",ifile.read(8))[0]

            self.__l_node_id.append(id)
            self.__l_node_sms_id.append(sms_id)
            self.__l_node_x.append(x)
            self.__l_node_y.append(y)
            self.__l_node_z.append(z)
            self.__l_node_nb_cell.append(cells_size)
            self.__l_node_nb_edge.append(edges_size)
        #== read nodes

        #-- read boundary edges
        self.__l_edge_bnd_id = []
        self.__l_edge_bnd_node_1_id = []
        self.__l_edge_bnd_node_2_id = []
        self.__l_edge_bnd_x_m = []
        self.__l_edge_bnd_y_m = []
        self.__l_edge_bnd_z_m = []
        self.__l_edge_bnd_lenght = []
        self.__l_edge_bnd_d_norm = []
        self.__l_edge_bnd_norm_x = []
        self.__l_edge_bnd_norm_y = []
        self.__l_edge_bnd_l_cell_id = []

        for i in range(self.__nb_edges_bnd) :
            id = struct.unpack("Q",ifile.read(8))[0]
            node_1_id = struct.unpack("Q",ifile.read(8))[0]
            node_2_id = struct.unpack("Q",ifile.read(8))[0]
            x_m = struct.unpack("d",ifile.read(8))[0]
            y_m = struct.unpack("d",ifile.read(8))[0]
            z_m = struct.unpack("d",ifile.read(8))[0]
            d_norm = struct.unpack("d",ifile.read(8))[0]
            lenght = struct.unpack("d",ifile.read(8))[0]
            norm_x = struct.unpack("d",ifile.read(8))[0]
            norm_y = struct.unpack("d",ifile.read(8))[0]
            bc_code = struct.unpack("i",ifile.read(4))[0]
            cell_id = struct.unpack("Q",ifile.read(8))[0]
            cells_id = [cell_id]

            self.__l_edge_bnd_id.append(id)
            self.__l_edge_bnd_node_1_id.append(node_1_id)
            self.__l_edge_bnd_node_2_id.append(node_2_id)
            self.__l_edge_bnd_x_m.append(x_m)
            self.__l_edge_bnd_y_m.append(y_m)
            self.__l_edge_bnd_z_m.append(z_m)
            self.__l_edge_bnd_lenght.append(lenght)
            self.__l_edge_bnd_d_norm.append(d_norm)
            self.__l_edge_bnd_norm_x.append(norm_x)
            self.__l_edge_bnd_norm_y.append(norm_y)
            self.__l_edge_bnd_l_cell_id.append(cells_id)
        #== read boundary edges

        #-- read internal edges
        self.__l_edge_int_id = []
        self.__l_edge_int_node_1_id = []
        self.__l_edge_int_node_2_id = []
        self.__l_edge_int_x_m = []
        self.__l_edge_int_y_m = []
        self.__l_edge_int_z_m = []
        self.__l_edge_int_lenght = []
        self.__l_edge_int_d_norm = []
        self.__l_edge_int_norm_x = []
        self.__l_edge_int_norm_y = []
        self.__l_edge_int_l_cell_id = []

        for i in range(self.__nb_edges_int) :
            id = struct.unpack("Q",ifile.read(8))[0]
            node_1_id = struct.unpack("Q",ifile.read(8))[0]
            node_2_id = struct.unpack("Q",ifile.read(8))[0]
            x_m = struct.unpack("d",ifile.read(8))[0]
            y_m = struct.unpack("d",ifile.read(8))[0]
            z_m = struct.unpack("d",ifile.read(8))[0]
            d_norm = struct.unpack("d",ifile.read(8))[0]
            lenght = struct.unpack("d",ifile.read(8))[0]
            norm_x = struct.unpack("d",ifile.read(8))[0]
            norm_y = struct.unpack("d",ifile.read(8))[0]
            
            cells_id = []
            for j in range(2) :
                cell_id = struct.unpack("Q",ifile.read(8))[0]
                cells_id.append(cell_id)

            self.__l_edge_int_id.append(id)
            self.__l_edge_int_node_1_id.append(node_1_id)
            self.__l_edge_int_node_2_id.append(node_2_id)
            self.__l_edge_int_x_m.append(x_m)
            self.__l_edge_int_y_m.append(y_m)
            self.__l_edge_int_z_m.append(z_m)
            self.__l_edge_int_lenght.append(lenght)
            self.__l_edge_int_d_norm.append(d_norm)
            self.__l_edge_int_norm_x.append(norm_x)
            self.__l_edge_int_norm_y.append(norm_y)
            self.__l_edge_int_l_cell_id.append(cells_id)
        #== read internal edges

        #-- read boundaries
        self.__l_boundaries = []
        self.__l_bc_code = []

        for i in range(self.__nb_boundaries) :
            bnd_size = struct.unpack("Q",ifile.read(8))[0]
            boundary = []
            for j in range(bnd_size) :
                edge_idx = struct.unpack("Q",ifile.read(8))[0]
                boundary.append(edge_idx)
            self.__l_boundaries.append(boundary)

        for i in range(self.__nb_boundaries) :
            code = struct.unpack("i",ifile.read(4))[0]
            self.__l_bc_code.append(code)
        #== read boundaries

        #-- read cell link with other elements
        self.__l_cell_l_node_id = []
        self.__l_cell_l_edge_id = []
        self.__l_cell_l_neigh_id = []

        for i in range(self.__nb_cells) :
            nodes = []
            for j in range(self.__l_cell_nb_node[i]) :
                node_id = struct.unpack("Q",ifile.read(8))[0]
                nodes.append(node_id)
            self.__l_cell_l_node_id.append(nodes)

            edges_bnd = []
            for j in range(self.__l_cell_nb_edge_bnd[i]) :
                edge_id = struct.unpack("Q",ifile.read(8))[0]
                edges_bnd.append(edge_id)

            edges_int = []
            for j in range(self.__l_cell_nb_edge_int[i]) :
                edge_id = struct.unpack("Q",ifile.read(8))[0]
                edges_int.append(edge_id)
        
            self.__l_cell_l_edge_bnd_id.append(edges_bnd)
            self.__l_cell_l_edge_int_id.append(edges_int)

            neigh = []
            for j in range(self.__l_cell_nb_neigh[i]) :
                neigh_id = struct.unpack("Q",ifile.read(8))[0]
                neigh.append(neigh_id)
            self.__l_cell_l_neigh_id.append(neigh)
        #== read cell link with other elements

        #-- read node link with other elements
        self.__l_node_l_cell_id = []

        for i in range(self.__nb_nodes) :
            cells = []
            for j in range(self.__l_node_nb_cell[i]) :
                cell_id = struct.unpack("Q",ifile.read(8))[0]
                cells.append(cell_id)
            self.__l_node_l_cell_id.append(cells)
        #== read node link with other elements
        #=== Read geometry

        #--- read variable names
        self.__nb_var = struct.unpack("Q",ifile.read(8))[0]
        self.__l_var_names = []

        for i in range(self.__nb_var) :
            name_lenght = struct.unpack("Q",ifile.read(8))[0]
            name = ifile.read(name_lenght).decode("utf-8")
            self.__l_var_names.append(name)
        #=== read variable names

        #--- read time information
        self.__t_0 = struct.unpack("d",ifile.read(8))[0]
        self.__t_max = struct.unpack("d",ifile.read(8))[0]
        self.__dt_map = struct.unpack("d",ifile.read(8))[0]
        #=== read time information

        #--- read simulation results
        self.__times_maps = []

        l_hydro = []

        time_bin = ifile.read(8)

        while (time_bin) : 
            time = struct.unpack("d", time_bin)[0]
            self.__times_maps.append(time)
            #-- read result map
            l_result = []

            for i in range(self.__nb_cells) :
                l_cell_result = []
                for j in range(self.__nb_var) :
                    value_bytes = ifile.read(8)
                    value = struct.unpack("d",value_bytes)[0]
                    l_cell_result.append(value)
                l_result.append(l_cell_result)

            l_hydro.append(l_result)
            #== read result map

            time_bin = ifile.read(8)
        #=== read simulation results

        self.__result_hydro = np.array(l_hydro)