from netCDF4 import Dataset
import numpy as np
import sys

import readBinResult

def create_netCDF(result, folder) :
    """Create a NetCDF file containing the sw2d simulation result
    
    Keyword arguments :
    result -- an object of readBinResult.BinResult class containing the simulation result
    folder -- string containing the name of the simulation folder where the file will be created
    """
    # Create file
    filename = folder + "/output/ResultHydro.nc"
    file = Dataset(filename, 'w', format='NETCDF4')

    #--- Add dimensions
    nb_timestep = len(result.getTimesMaps())
    file.createDimension('time', nb_timestep)
    file.createDimension('mesh2d_nFaces', result.getNbCells())
    file.createDimension('mesh2d_nEdges', result.getNbEdges())
    file.createDimension('Two', 2)
    file.createDimension('mesh2d_nMax_face_nodes', 4)
    file.createDimension('mesh2d_nNodes', result.getNbNodes())
    #=== Add dimensions

    #--- Add geometry variables
    #-- mesh
    mesh2d = file.createVariable('mesh2d', 'i4', ())
    mesh2d.cf_role = 'mesh_topology'
    mesh2d.long_name = 'Topology data of 2D mesh'
    mesh2d.topology_dimension = 2
    mesh2d.node_coordinates = 'mesh2d_node_x mesh2d_node_y'
    mesh2d.node_dimension = 'mesh2d_nNodes'
    mesh2d.max_face_nodes_dimension = 'mesh2d_nMax_face_nodes'
    mesh2d.edge_node_connectivity = 'mesh2d_edge_nodes'
    mesh2d.edge_dimension = 'mesh2d_nEdges'
    mesh2d.edge_coordinates = 'mesh2d_edge_x mesh2d_edge_y'
    mesh2d.face_node_connectivity = 'mesh2d_face_nodes'
    mesh2d.face_dimension = 'mesh2d_nFaces'
    mesh2d.edge_face_connectivity = 'mesh2d_edge_faces'
    mesh2d.face_coordinates = 'mesh2d_face_x mesh2d_face_y'
    #== mesh

    #-- nodes
    mesh2d_node_x = file.createVariable('mesh2d_node_x', 'f8', ('mesh2d_nNodes',))
    mesh2d_node_x.units = 'm'
    mesh2d_node_x.standard_name = 'projection_x_coordinate'
    mesh2d_node_x.long_name = 'X-coordinate of mesh nodes'
    mesh2d_node_x[:] = result.getNodesX()

    mesh2d_node_y = file.createVariable('mesh2d_node_y', 'f8', ('mesh2d_nNodes',))
    mesh2d_node_y.units = 'm'
    mesh2d_node_y.standard_name = 'projection_y_coordinate'
    mesh2d_node_y.long_name = 'Y-coordinate of mesh nodes'
    mesh2d_node_y[:] = result.getNodesY()

    mesh2d_node_z = file.createVariable('mesh2d_node_z', 'f8', ('mesh2d_nNodes',))
    mesh2d_node_z.mesh = 'mesh2d'
    mesh2d_node_z.location = 'node'
    mesh2d_node_z.coordinates = 'mesh2d_node_x mesh2d_node_y'
    mesh2d_node_z.standard_name = 'altitude'
    mesh2d_node_z.long_name = 'Z-coordinate of mesh nodes'
    mesh2d_node_z.units = 'm'
    mesh2d_node_z[:] = result.getNodesZ()
    #== nodes

    #-- edges
    mesh2d_edge_faces = file.createVariable('mesh2d_edge_faces', 'i4', ('mesh2d_nEdges', 'Two'))
    mesh2d_edge_faces.cf_role = 'edge_face_connectivity'
    mesh2d_edge_faces.long_name = 'Neighboring faces of mesh edges'
    mesh2d_edge_faces.start_index = 1
    edge_faces = []
    for cells_id in result.getEdgesCells() :
        cells_id_nc = [id+1 for id in cells_id]
        if (len(cells_id_nc) == 1) :
            edge_faces.append([cells_id_nc[0], -1])
        else :
            edge_faces.append(cells_id_nc)
    mesh2d_edge_faces[:] = np.array(edge_faces)

    l_edge_node_1_id_nc = [id+1 for id in result.getEdgesNode1()]
    l_edge_node_2_id_nc = [id+1 for id in result.getEdgesNode2()]
    mesh2d_edge_nodes = file.createVariable('mesh2d_edge_nodes', 'i4', ('mesh2d_nEdges', 'Two'))
    mesh2d_edge_nodes.cf_role = 'edge_node_connectivity'
    mesh2d_edge_nodes.long_name = 'Start and end nodes of mesh edges'
    mesh2d_edge_nodes.start_index = 1
    mesh2d_edge_nodes[:] = np.array([l_edge_node_1_id_nc, l_edge_node_2_id_nc]).transpose()

    mesh2d_edge_x = file.createVariable('mesh2d_edge_x', 'f8', ('mesh2d_nEdges',))
    mesh2d_edge_x.units = 'm'
    mesh2d_edge_x.standard_name = 'projection_x_coordinate'
    mesh2d_edge_x.long_name = 'X-coordinate of mesh edges'
    mesh2d_edge_x[:] = result.getEdgesX()

    mesh2d_edge_y = file.createVariable('mesh2d_edge_y', 'f8', ('mesh2d_nEdges',))
    mesh2d_edge_y.units = 'm'
    mesh2d_edge_y.standard_name = 'projection_y_coordinate'
    mesh2d_edge_y.long_name = 'Y-coordinate of mesh edges'
    mesh2d_edge_y[:] = result.getEdgesY()
    #== edges

    #-- cells
    mesh2d_face_nodes = file.createVariable('mesh2d_face_nodes', 'i4', ('mesh2d_nFaces', 'mesh2d_nMax_face_nodes'))
    mesh2d_face_nodes.cf_role = 'face_node_connectivity'
    mesh2d_face_nodes.long_name = 'Vertex nodes of mesh faces (counterclockwise)'
    mesh2d_face_nodes.start_index = 1
    face_nodes = []
    for nodes_id in result.getCellsNodes() :
        nodes_id_nc = [id+1 for id in nodes_id]
        if (len(nodes_id_nc) == 3) :
            face_nodes.append([nodes_id_nc[0], nodes_id_nc[1], nodes_id_nc[2], -1])
        elif (len(nodes_id_nc) == 4):
            face_nodes.append(nodes_id_nc)
    mesh2d_face_nodes[:] = np.array(face_nodes)

    mesh2d_face_x = file.createVariable('mesh2d_face_x', 'f8', ('mesh2d_nFaces',))
    mesh2d_face_x.units = 'm'
    mesh2d_face_x.standard_name = 'projection_x_coordinate'
    mesh2d_face_x.long_name = 'X-coordinate of mesh faces'
    mesh2d_face_x[:] = result.getCellsX()

    mesh2d_face_y = file.createVariable('mesh2d_face_y', 'f8', ('mesh2d_nFaces',))
    mesh2d_face_y.units = 'm'
    mesh2d_face_y.standard_name = 'projection_y_coordinate'
    mesh2d_face_y.long_name = 'Y-coordinate of mesh faces'
    mesh2d_face_y[:] = result.getCellsY()
    #== cells
    #=== Add geometry variables

    #--- Add result variables
    #-- result maps
    for i in range(result.getNbVar()) :
        var = result.getVarNames()[i]
        if (var[-2:] == "_x") :
            mesh2d_var = file.createVariable('mesh2d_' + var, 'f8', ('mesh2d_nFaces','time'))
            mesh2d_var.mesh = 'mesh2d'
            mesh2d_var.location = 'face'
            mesh2d_var.coordinates = 'mesh2d_face_x mesh2d_face_y'
            mesh2d_var.cell_methods = 'mesh2d_' + var + ': mean'
            mesh2d_var.long_name = var[:-2] + ", x-component"
            mesh2d_var[:] = result.getResultHydro()[:,:,i].transpose()
        elif (var[-2:] == "_y") :
            mesh2d_var = file.createVariable('mesh2d_' + var, 'f8', ('mesh2d_nFaces','time'))
            mesh2d_var.mesh = 'mesh2d'
            mesh2d_var.location = 'face'
            mesh2d_var.coordinates = 'mesh2d_face_x mesh2d_face_y'
            mesh2d_var.cell_methods = 'mesh2d_' + var + ': mean'
            mesh2d_var.long_name = var[:-2] + ", y-component"
            mesh2d_var[:] = result.getResultHydro()[:,:,i].transpose()
        elif (var[-5:] != "_norm") :
            mesh2d_var = file.createVariable('mesh2d_' + var, 'f8', ('mesh2d_nFaces','time'))
            mesh2d_var.mesh = 'mesh2d'
            mesh2d_var.location = 'face'
            mesh2d_var.coordinates = 'mesh2d_face_x mesh2d_face_y'
            mesh2d_var.cell_methods = 'mesh2d_' + var + ': mean'
            mesh2d_var.long_name = var
            mesh2d_var[:] = result.getResultHydro()[:,:,i].transpose()
    #== result maps
    #=== Add result variables

    #--- Add time variable
    time = file.createVariable('time', 'f8', ('time',))
    time.standard_name = 'time'
    time.units = 'seconds since 2025-01-01 00:00:00'
    time[:] = [i*result.getDtMap() for i in range(nb_timestep)]
    #=== Add time variables

    # Close file
    file.close()