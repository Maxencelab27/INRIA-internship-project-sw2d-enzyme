import numpy as np
import csv 
import sys

### Table of elements : can be completed ###
Table_Element_Type = {1 : "2-nodes line", 2 : "3-nodes triangle", 3 : "4-nodes quad"} #Table implement as memo
Table_MshElement_2dmElement = {1 : "E2L", 2 : "E3T", 3 : "E4Q"}    #Match with msh code and 2dm code, see
# http://www.manpagez.com/info/gmsh/gmsh-2.2.6/gmsh_63.php and https://www.xmswiki.com/wiki/SMS:2D_Mesh_Files_*.2dm
Table_Element_NbrNodes = {1: 2, 2 : 3, 3 : 4}      #Number of node for each type of elements

class Mesh:
    """
    Store the nodes and elements

    Attributes
    ----------
    Nodes : array
        array of 3d coordinates (npts x 3)
    Id2Ind_Nodes : dict
        Store the correspondance between msh nodes id and their positions in the array Nodes
        key : id, value : ind
    Elments : dict
        dictionary of array. Each column of array corresponds to id of nodes which form an element. One column per element.
        key : type of element, value : array of nodes id
    NbrElements : dict
        Store the amounts of elements per type.
        key : type, value : number of element

    /!\ Physical Entity are ignored. I don't know what they are used for.
    """
    def __init__(self) -> None:
        self.msh_filename = ""

        self.Nodes = np.array([])
        self.Id2Ind_Nodes = {}

        self.Elements = {}
        self.NbrElements = {}

    def read_msh(self, mshfile):
        """
        Input : a mshfile name, in format ASCII 2.2
        Fill self with the data find inside mshfile.
        /!\ Only the elements present inside Table_Element_Type are stored. Other are ignored.
        """
        self.msh_filename = mshfile
        print(self.msh_filename)

        with open(self.msh_filename, 'r') as csvfile:
            rows = csv.reader(csvfile, delimiter=' ')

            ###--- Check type and format of file ---###
            correct_file = False

            header = next(rows)
            print(header)
            if (header[0] == "$MeshFormat"):
                format = next(rows)
                print(format)
                if (format[0] == "2.2"):
                    correct_file = True
            
            if (not correct_file):
                raise ValueError("The mesh file is not a msh file in format 2.2 ASCI")

            end_header = next(rows)
            not_eof = True
            try:
                row = next(rows)
            except StopIteration:
                print("The file contains no mesh data")
                not_eof = False
            ###=== Check type and format of file ===###

            while not_eof:
                
                ###--- Store Nodes ---###
                if(row[0] == "$Nodes"):
                    nbr_nodes = int(next(rows)[0])
                    print("Number of Nodes : ", nbr_nodes)
                    self.Nodes = np.zeros((nbr_nodes, 3), dtype=float)

                    for ind in range(0, nbr_nodes):
                        node = next(rows)
                        id_node = int(node[0])
                        self.Id2Ind_Nodes[id_node] = ind
                        self.Nodes[ind, :] = [float(d) for d in node[1:]]

                    row = next(rows)
                    if(row[0] != "$EndNodes"):
                        raise ValueError("Expecting '$EndNodes'")
                ###=== Store Nodes ===###
                
                ###--- Store Elements ---###
                if (row[0] == "$Elements"):
                    nbr_elements = int(next(rows)[0])
                    print("Number of Elements : ", nbr_nodes)

                    for ind in range(0, nbr_elements):
                        element = next(rows)
                        id_element = int(element[0])
                        self.Id2Ind_Nodes[id_element] = ind

                        element_type = int(element[1])
                        if (element_type not in Table_Element_Type):
                            print("Type of element : ", element_type, " unknown ; This element is ignored")
                        else:
                            nbr_tags = int(element[2])
                            element_nodes = [self.Id2Ind_Nodes[int(id_node)] for id_node in element[3+nbr_tags:]]

                            if (element_type not in self.Elements):
                                self.Elements[element_type] = element_nodes
                                self.NbrElements[element_type] = 1
                            else:
                                self.Elements[element_type] = np.vstack((self.Elements[element_type], element_nodes))
                                self.NbrElements[element_type] += 1
                            
                    row = next(rows)
                    if(row[0] != "$EndElements"):
                        raise ValueError("Expecting '$EndElements'")
                ###=== Store Elements ===###
                    
                try:
                    row = next(rows)
                except StopIteration:
                    print("Reach Last Line")
                    not_eof = False

    
    def write_2dm(self, element_to_write = ["E3T", "E4Q"], filename=None):
        """
        Input : Eventually a list of 2dm element we want. by default, only E3T elements will be written. eventually the name of output file. The default name is mshfilename.2dm
        Use all the data stored to write them into a .2dm file.
        """
        ###--- Find filename ---###
        if filename is None:
            filename = self.msh_filename.split('.')[0] + '.2dm'
            print("### Write in : ", filename)
        ###=== Find filename ===###

        with open(filename, 'w', newline='') as csvfile:
            spamwriter = csv.writer(csvfile, delimiter='\t')
            ###--- Write Header ---###
            spamwriter.writerow(["MESH2D"])
            ###=== Write Header ===###

            ###--- Write Elements ---###
            print(self.Elements)
            print(element_to_write)
            for element_type in self.Elements:

                ###--- If this type of element has to be written
                if Table_MshElement_2dmElement[element_type] in element_to_write:
                    ### iterates through all elements of this type
                    for element_ind in range(self.NbrElements[element_type]):
                        element_data = [element_ind + 1] #id of this element
                        element_data.extend((self.Elements[element_type][element_ind] + 1).tolist()) #id of nodes which composed this elements
                        element_data.append(1)
                        row = [Table_MshElement_2dmElement[element_type]]   #2dm code
                        row.extend([str(d) for d in element_data])          #id of nodes
                        spamwriter.writerow(row)
            ###=== Write Elements ===###

            ###--- Write Nodes ---###
            for ind, node in enumerate(self.Nodes):
                row = ["ND", str(ind + 1)] #type and id of node; 2dm indices must start at 1
                row.extend([str(coord) for coord in node]) #coordinates of each node
                spamwriter.writerow(row)
            ###=== Write Nodes ===###


def main(filename):
    my_mesh = Mesh()
    print("### Begin read ###")
    my_mesh.read_msh(filename)
    print("### End read ###")
    print("### Begin write ###")
    my_mesh.write_2dm()
    print("### End write ###")

if __name__ == "__main__":
    args = sys.argv

    if (len(args) == 2):
        main(args[1])
    elif (len(args) == 3):
        main(args[1], args[2])
    else:
        print("Use : converterMsh_2dm.py input.msh [optionnal : output.2dm]")
            



