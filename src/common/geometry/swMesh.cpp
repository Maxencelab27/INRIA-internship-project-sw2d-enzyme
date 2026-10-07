#include "swMesh.h"

#include <algorithm>
#include <cmath>
#include <filesystem>

#include "../../common/swFile.h"
#include "swCell.h"
#include "swEdgeBoundary.h"
#include "swEdgeIntern.h"
#include "../../common/swIO.h"
#include "../../common/swLog.h"
#include "swNode.h"
#include "../../common/param.h"
#include "../../common/swString.h"

swMesh::swMesh(swIO *io)
{
    m_io = io;
}

bool swMesh::assignBC(std::string filename_bc, size_t default_bc)
/**
 * @brief Read a .bc file, set list_bc_codes attribute and boundary edges bc number
 *
 * @param[in] filename_bc name of the .bc file
 * @param[in] default_bc default value of boundary edges bc number
 * @return bool EXIT_SUCCES/EXIT_FAILURE
 */
{
    //-- read th .bc file
    //---- Open .bc file
    std::ifstream file(m_io->absolutePath(filename_bc));
    if (file.is_open() == false)
    {
        logError() << "  Unable to open file" << m_io->absolutePath(filename_bc);
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "  Reading bc file";
    //==== Open .bc file

    int line_nb = 1;
    bool conversion_success = true;

    //---- Read number of lines
    std::string line_txt = swIO::readLine(file);
    const size_t nb_lines = swString::stringToUInt(line_txt, conversion_success);
    if (!conversion_success)
    {
        logError() << "    line is misformatted: expecting a positive integer (nb of lines to read)";
        logFatal(__FILE__, __func__, __LINE__);
    }
    //==== Read number of lines

    //---- Read the user-defined boundary data
    std::vector<int> read_codes; // vector of the bc_code for user-defined boundaries
    std::vector<std::vector<size_t>> read_node_pairs; // vector of the 2 points defining a boundary
    for (size_t k = 0; k < nb_lines; ++k)
    {
        ++line_nb;
        //-- Read line
        line_txt = swIO::readLine(file);
        const std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
        if (line_split.size() != 3)
        {
            logError() << "    line" << line_nb << "is misformatted : expecting 2 positive integers (nodes id) and and a 1 negative (bc code)";
            logFatal(__FILE__, __func__, __LINE__);
        }

        size_t node_1 = swString::stringToUInt(line_split[0], conversion_success);
        size_t node_2 = swString::stringToUInt(line_split[1], conversion_success);
        int code_neg = swString::stringToInt(line_split[2], conversion_success);
        if ((!conversion_success) || (code_neg >= 0))
        {
            logError() << "    line" << line_nb << "is misformatted : expecting 2 positive integers (nodes id) and a 1 negative (bc code)";
            logFatal(__FILE__, __func__, __LINE__);
        }
        //== Read line

        // Fill lists with read elements
        read_node_pairs.push_back({node_1, node_2});
        read_codes.push_back(code_neg);
    }
    //==== Read the user-defined boundary data
    file.close();
    //== read th .bc file


    //-- assign bc code to the boundary edges
    std::vector<bool> bc_assigned = std::vector<bool>(nb_lines, false); // vector of the bc_code for user-defined boundaries
    size_t edge_assigned = 0;
    for (swEdgeBoundary *edge : list_edges_boundary)
    {
        for (size_t l = 0; l < nb_lines; l++)
        {
            std::vector<swNode *> nodes = edge->getNodes();
            if (((nodes[0]->getSMSIdx() == read_node_pairs[l][0]) &&
                 (nodes[1]->getSMSIdx() == read_node_pairs[l][1])) ||
                ((nodes[0]->getSMSIdx() == read_node_pairs[l][1]) &&
                 (nodes[1]->getSMSIdx() == read_node_pairs[l][0])))
            {
                // boundary is user-defined
                ++edge_assigned;
                bc_assigned[l] = true;
                int bc_code = read_codes[l];
                edge->setBCCode(bc_code);
                break;
            }
        }
    }

    if (edge_assigned < nb_lines)
    {
        logError() << "User defined boundaries (in" << filename_bc << ") does not correspond to existing boundary edges";
        for (size_t l = 0; l < nb_lines; l++)
        {
            if (bc_assigned[l] == false)
            {
                logError() << "- No boundary edge found with nodes"
                           << read_node_pairs[l][0] << "and" << read_node_pairs[l][1]
                           << "\t(at line" << l + 2 << ")";
            }
        }
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== assign bc code to the boundary edges

    //-- Fill boundary list and set BC code
    for (swEdgeBoundary *edge : list_edges_boundary)
    {
        int bc_code = edge->getBCCode();
        if (std::find(list_bc_codes.begin(), list_bc_codes.end(), bc_code) == list_bc_codes.end())
        {
            // bc_code not in the list 'list_bc_codes'
            list_bc_codes.push_back(bc_code);
            std::vector<swEdgeBoundary *> boundary;
            boundary.push_back(edge);
            list_boundaries.push_back(boundary);
        }
        else
        {
            size_t rank = std::distance(std::begin(list_bc_codes), std::find(list_bc_codes.begin(), list_bc_codes.end(), bc_code)); // get the rank of bc_code in list_bc_codes
            list_boundaries[rank].push_back(edge);
        }
    }
    //== Fill boundary list and set BC code

    size_t nb_bc = list_bc_codes.size();
    logInfo() << "   " << nb_bc << "BC identified";
    logInfo() << "  ...done";

    return EXIT_SUCCESS;
}

bool swMesh::build()
/**
 * @brief read the geometry file and set the mesh attributes
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- import existing geometry file
    logInfo() << "Reading mesh from:" << swFile::input_geo_file;
    if (readGeometryFile() != EXIT_SUCCESS)
    {
        logError() << "Unable to read mesh from:" << swFile::input_geo_file;
        return EXIT_FAILURE;
    }

    logInfo() << "   " << list_cells.size() << "cells";
    logInfo() << "   " << list_nodes.size() << "nodes";
    logInfo() << "   " << list_edges_intern.size() << "internal edges";
    logInfo() << "   " << list_edges_boundary.size() << "boundary edges";
    logInfo() << "...done";
    //== import existing geometry file
    return EXIT_SUCCESS;
}

bool swMesh::readGeometryFile()
/**
 * @brief read the mesh from file geometry.sw2d
 * @return boolean indicating successful read or not
 */
{
    // Open the file
    std::ifstream file;
    file.open(m_io->absolutePath(swFile::input_geo_file), std::ios::binary);
    if (file.is_open() == false)
    {
        logError() << "The file does not exist" << m_io->absolutePath(swFile::input_geo_file);
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    // Read the first line (version number)
    size_t version_lenght;
    file.read(reinterpret_cast<char*>(&version_lenght), sizeof(version_lenght));
    std::string read_version(version_lenght, '\0');
    file.read(&read_version[0], version_lenght);
    m_version = read_version;
    // Read the content with the reader corresponding to the version
    if (m_version == "2")
    {
        return readGeometryFileV2(file); /** @todo Check wether this version numbering is necessary  */
    }
    else
    {
        logError() << "Unknown version number" << m_version;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swMesh::readGeometryFileV2(std::ifstream &file)
/**
 * @brief read the mesh from file geometry.sw2d (version 1 of geometry.sw2d)
 * @param[in] file: geometry file stream
 * @return boolean indicating successful read or not
 */
{
    //-- read number of elements
    size_t nb_cell;
    file.read(reinterpret_cast<char*>(&nb_cell), sizeof(nb_cell));
    size_t nb_edges_boundary;
    file.read(reinterpret_cast<char*>(&nb_edges_boundary), sizeof(nb_edges_boundary));
    size_t nb_edges_intern;
    file.read(reinterpret_cast<char*>(&nb_edges_intern), sizeof(nb_edges_intern));
    size_t nb_node;
    file.read(reinterpret_cast<char*>(&nb_node), sizeof(nb_node));
    size_t nb_boundaries;
    file.read(reinterpret_cast<char*>(&nb_boundaries), sizeof(nb_boundaries));
    //== read number of elements

    //-- read cells
    for (size_t n = 0; n < nb_cell; ++n)
    {
        swCell* cell = new swCell();
        size_t idx;
        size_t sms_id;
        size_t nb_node;
        size_t nb_edge_bnd;
        size_t nb_edge_int;
        size_t nb_neigh;
        double area;
        double perimeter;
        double x;
        double y;
        double zb;
        file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
        file.read(reinterpret_cast<char*>(&sms_id), sizeof(sms_id));
        file.read(reinterpret_cast<char*>(&nb_node), sizeof(nb_node));
        file.read(reinterpret_cast<char*>(&nb_edge_bnd), sizeof(nb_edge_bnd));
        file.read(reinterpret_cast<char*>(&nb_edge_int), sizeof(nb_edge_int));
        file.read(reinterpret_cast<char*>(&nb_neigh), sizeof(nb_neigh));
        file.read(reinterpret_cast<char*>(&area), sizeof(area));
        file.read(reinterpret_cast<char*>(&perimeter), sizeof(perimeter));
        file.read(reinterpret_cast<char*>(&x), sizeof(x));
        file.read(reinterpret_cast<char*>(&y), sizeof(y));
        file.read(reinterpret_cast<char*>(&zb), sizeof(zb));
        cell->setProperties(idx, sms_id, x, y, zb, area, perimeter,
                            nb_node, nb_edge_bnd, nb_edge_int, nb_neigh);
        list_cells.push_back(cell);
    }
    logInfo() << "  cells ...done";
    //== read cells

    //-- read nodes
    std::vector<size_t> nb_edges_per_point;
    for (size_t p = 0; p < nb_node; ++p)
    {
        swNode* node = new swNode();
        size_t id;
        size_t sms_id;
        size_t nb_cell;
        size_t nb_edge;
        double x;
        double y;
        double z;
        file.read(reinterpret_cast<char*>(&id), sizeof(id));
        file.read(reinterpret_cast<char*>(&sms_id), sizeof(sms_id));
        file.read(reinterpret_cast<char*>(&x), sizeof(x));
        file.read(reinterpret_cast<char*>(&y), sizeof(y));
        file.read(reinterpret_cast<char*>(&z), sizeof(z));
        file.read(reinterpret_cast<char*>(&nb_cell), sizeof(nb_cell));
        file.read(reinterpret_cast<char*>(&nb_edge), sizeof(nb_edge));
        node->setCoord(sms_id, id, x, y, z);
        list_nodes.push_back(node);
        nb_edges_per_point.push_back(nb_edge);
    }
    logInfo() << "  nodes ...done";
    //== read nodes

    //-- read boundary edges
    for (size_t i = 0; i < nb_edges_boundary; ++i)
    {
        size_t idx;
        size_t node_1_idx;
        size_t node_2_idx;
        int bc_code;
        double x_m;
        double y_m;
        double z_m;
        double d_norm;
        double length;
        double norm_x;
        double norm_y;

        file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
        file.read(reinterpret_cast<char*>(&node_1_idx), sizeof(node_1_idx));
        file.read(reinterpret_cast<char*>(&node_2_idx), sizeof(node_2_idx));
        file.read(reinterpret_cast<char*>(&x_m), sizeof(x_m));
        file.read(reinterpret_cast<char*>(&y_m), sizeof(y_m));
        file.read(reinterpret_cast<char*>(&z_m), sizeof(z_m));
        file.read(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
        file.read(reinterpret_cast<char*>(&length), sizeof(length));
        file.read(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
        file.read(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));
        file.read(reinterpret_cast<char*>(&bc_code), sizeof(bc_code));

        swNode* node_1 = list_nodes[node_1_idx];
        swNode* node_2 = list_nodes[node_2_idx];
        swEdgeBoundary* edge_bnd = new swEdgeBoundary();
        edge_bnd->setProperties(idx, node_1, node_2,
                                x_m, y_m, z_m,
                                d_norm, length, norm_x, norm_y,
                                bc_code);

        size_t cell_idx;
        file.read(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx));
        edge_bnd->setCell(list_cells[cell_idx]);
        list_edges_boundary.push_back(edge_bnd);
    }
    logInfo() << "  boundary edges ...done";
    //== read boundary edges

    //-- read intern edges
    for (size_t i = 0; i < nb_edges_intern; ++i)
    {
        swEdgeIntern* edge_int = new swEdgeIntern();
        size_t idx;
        size_t node_1_idx;
        size_t node_2_idx;
        double x_m;
        double y_m;
        double z_m;
        double length;
        double d_norm;
        double norm_x;
        double norm_y;

        file.read(reinterpret_cast<char*>(&idx), sizeof(idx));
        file.read(reinterpret_cast<char*>(&node_1_idx), sizeof(node_1_idx));
        file.read(reinterpret_cast<char*>(&node_2_idx), sizeof(node_2_idx));
        file.read(reinterpret_cast<char*>(&x_m), sizeof(x_m));
        file.read(reinterpret_cast<char*>(&y_m), sizeof(y_m));
        file.read(reinterpret_cast<char*>(&z_m), sizeof(z_m));
        file.read(reinterpret_cast<char*>(&d_norm), sizeof(d_norm));
        file.read(reinterpret_cast<char*>(&length), sizeof(length));
        file.read(reinterpret_cast<char*>(&norm_x), sizeof(norm_x));
        file.read(reinterpret_cast<char*>(&norm_y), sizeof(norm_y));

        swNode* node_1 = list_nodes[node_1_idx];
        swNode* node_2 = list_nodes[node_2_idx];
        edge_int->setProperties(idx, node_1, node_2,
                                x_m, y_m, z_m,
                                d_norm, length, norm_x, norm_y);

        std::vector<swCell*> edge_cells_list;
        for (size_t n = 0; n < 2 ; ++n)
        {
            size_t cell_idx;
            file.read(reinterpret_cast<char*>(&cell_idx), sizeof(cell_idx));
            edge_cells_list.push_back(list_cells[cell_idx]);
        }
        edge_int->setCells(edge_cells_list);
        list_edges_intern.push_back(edge_int);
    }
    logInfo() << "  internal edges ...done";
    //== read intern edges

    //-- read boundaries
    for (size_t k = 0; k < nb_boundaries; k++)
    {
        size_t bnd_size;
        file.read(reinterpret_cast<char*>(&bnd_size), sizeof(bnd_size));

        std::vector<swEdgeBoundary *> boundary;
        for (size_t i = 0; i < bnd_size; i++)
        {
            size_t edge_idx;
            file.read(reinterpret_cast<char*>(&edge_idx), sizeof(edge_idx));
            for (swEdgeBoundary *edge : list_edges_boundary)
            {
                if (edge->getIdx() == edge_idx)
                {
                    boundary.push_back(edge);
                    break;
                }
            }
        }
        list_boundaries.push_back(boundary);
    }

    for (size_t k = 0; k < nb_boundaries; k++)
    {
        int code;
        file.read(reinterpret_cast<char*>(&code), sizeof(code));
        list_bc_codes.push_back(code);
    }
    //== read boundaries

    //-- read cell link with other elements
    for (swCell* cell : list_cells)
    {
        std::vector<swNode*> cell_nodes;
        for (size_t p = 0; p < cell->getNodes().size(); ++p)
        {
            size_t node_idx;
            file.read(reinterpret_cast<char*>(&node_idx), sizeof(node_idx));
            cell_nodes.push_back(list_nodes[node_idx]);
        }
        cell->setNodes(cell_nodes);

        std::vector<swEdgeBoundary*> cell_edges_bnd;
        for (size_t i = 0; i < cell->getEdgesBoundary().size(); ++i)
        {
            size_t edge_idx;
            file.read(reinterpret_cast<char*>(&edge_idx), sizeof(edge_idx));
            cell_edges_bnd.push_back(list_edges_boundary[edge_idx]);
        }
        cell->setEdgesBnd(cell_edges_bnd);

        std::vector<swEdgeIntern*> cell_edges_int;
        for (size_t i = 0; i < cell->getEdgesIntern().size(); ++i)
        {
            size_t edge_idx;
            file.read(reinterpret_cast<char*>(&edge_idx), sizeof(edge_idx));
            cell_edges_int.push_back(list_edges_intern[edge_idx]);
        }
        cell->setEdgesInt(cell_edges_int);

        std::vector<swCell*> cell_neighbouring_cells(cell->getNeigh().size());
        for (size_t n = 0; n < cell->getNeigh().size(); ++n)
        {
            size_t neight_cell_idx;
            file.read(reinterpret_cast<char*>(&neight_cell_idx), sizeof(neight_cell_idx));
            cell_neighbouring_cells[n] = list_cells[neight_cell_idx];
        }
        cell->setNeigh(cell_neighbouring_cells);
    }
    logInfo() << "  cells connectivities ...done";
    //== read cell link with other elements

    //-- read node link with other elements
    for (size_t p = 0; p < nb_node; p++)
    {
        swNode* node = list_nodes[p];
        std::vector<swCell*> node_cells;
        for (size_t n_p = 0; n_p < node->getCells().size(); ++n_p)
        {
            size_t cel_idx;
            file.read(reinterpret_cast<char*>(&cel_idx), sizeof(cel_idx));
            node_cells.push_back(list_cells[cel_idx]);
        }
        node->setCells(node_cells);
    }
    logInfo() << "  nodes connectivities ...done";
    //== read node link with other elements
    file.close();

    return EXIT_SUCCESS;
}

bool swMesh::buildCells(std::vector<size_t> cells_sms_id, std::vector<std::vector<size_t>> cells_nodes_idx)
/**
 * @brief Create cell objects and set list_cells attributes
 * 
 * @param[in] cells_sms_id list of the cells id (SMS numbering style)
 * @param[in] cells_nodes_idx list of the cells nodes
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    for (size_t n = 0; n < cells_sms_id.size(); ++n)
    {
        //-- Find node pointers
        std::vector<swNode *> nodes_ptr;
        std::vector<size_t> nodes_idx = cells_nodes_idx[n];
        for (size_t p = 0; p < nodes_idx.size(); p++)
        {
            bool found = false;
            for (swNode* node : this->list_nodes)
            {
                if (node->getSMSIdx() == nodes_idx[p])
                {
                    nodes_ptr.push_back(node);
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                logError() << "  Node" << nodes_idx[p] << "for cell" << cells_sms_id[n] << "does not exist";
                logFatal(__FILE__, __func__, __LINE__);
            }
        }
        //== Find node pointers
        
        // Create a new cell and add it to the list
        swCell *cell = new swCell(cells_sms_id[n], n, nodes_ptr);
        list_cells.push_back(cell);
    }

    return EXIT_SUCCESS;
}

bool swMesh::buildEdges(size_t default_bc)
/**
 * @brief Create edge objects and set list_edges_boundary and list_edges_intern attributes
 * 
 * @param[in] default_bc default value of bc number to assign to all edges
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    //-- Compute max number of edges
    size_t nb_edges_max = 0;
    for (swCell *cell : this->list_cells)
    {
        nb_edges_max += cell->getNodes().size();
    }
    //== Compute max number of edges

    //-- Identify the edges
    std::vector<std::vector<swNode*>> edge_nodes(nb_edges_max);
    std::vector<std::vector<swCell*>> edge_cells(nb_edges_max);

    std::vector<std::vector<size_t>> node_edges(this->list_nodes.size()); // list of edges associated to each node

    size_t edge_id = 0;
    for (swNode* node : this->list_nodes)
    {
        for (swNode* neigh : node->getNeighNodes())
        {
            // for each pair of nodes
            size_t id_node  = node->getIdx();
            size_t id_neigh = neigh->getIdx();
            
            if (id_node < id_neigh)
            {
                // check if an edge with these two nodes already exists
                bool found = false;
                for (size_t edge_number : node_edges[id_node])
                {
                    if ((node  == edge_nodes[edge_number][0]) &&
                        (neigh == edge_nodes[edge_number][1]))
                    {
                        found = true;
                    }
                }

                if (found == false)
                {
                    // if not, create edge
                    node_edges[id_node].push_back(edge_id);
                    node_edges[id_neigh].push_back(edge_id);
                    edge_nodes[edge_id].push_back(node);
                    edge_nodes[edge_id].push_back(neigh);
                }

                for (swCell* cell : node->getCells())
                {
                    // for each cell attached to first node
                    std::vector<swNode*> cell_nodes = cell->getNodes();

                    for (size_t cell_node_idx = 0; cell_node_idx < cell_nodes.size(); ++cell_node_idx)
                    {
                        swNode *source = nullptr;
                        swNode *target = nullptr;
                        if (cell_node_idx == 0)
                        {
                            source = cell_nodes.back();
                            target = cell_nodes[0];
                        }
                        else
                        {
                            source = cell_nodes[cell_node_idx - 1];
                            target = cell_nodes[cell_node_idx];
                        }
                        // if the edge is in this cell, add cell to edge_cells
                        if (((source->getIdx() == id_node) && (target->getIdx() == neigh->getIdx())) ||
                            ((source->getIdx() == neigh->getIdx()) && (target->getIdx() == id_node)))
                        {
                            edge_cells[edge_id].push_back(cell);
                        }
                    }
                }
                if (found == false)
                {
                    // if edge has been created, increment id
                    edge_id++;
                }
            }
        }
    }

    size_t nb_edges = edge_id;
    logInfo() << "   " << nb_edges << "edges identified";
    //== Identify the edges

    //-- Create the edges
    size_t edge_bnd_idx = 0;
    size_t edge_int_idx = 0;

    for (size_t i = 0; i < nb_edges; ++i)
    {
        if (edge_cells[i].size() == 1) // boundary interface
        {
            swEdgeBoundary* edge_boundary = new swEdgeBoundary(edge_bnd_idx, edge_nodes[i], edge_cells[i][0], default_bc);
            this->list_edges_boundary.push_back(edge_boundary);
            edge_bnd_idx++;
        }
        else // internal interface
        {
            swEdgeIntern* edge_intern = new swEdgeIntern(edge_int_idx, edge_nodes[i], edge_cells[i]);
            this->list_edges_intern.push_back(edge_intern);
            edge_int_idx++;
        }
    }
    //== Create the edges

    return EXIT_SUCCESS;
}

bool swMesh::checkQuality(double a_min)
/**
 * @brief check cells area and write smaller cells in a log file
 * 
 * @param[in] a_min area below which cells are considered small
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    if (a_min > 0)
    {
        // Write all cells log file
        if (writeLogCell(list_cells, m_io->absolutePath(swFile::log_cell_file)) != EXIT_SUCCESS)
        {
            logError() << "  Error in writeLogCell()";
            logFatal(__FILE__, __func__, __LINE__);
        }

        // Get small cells
        size_t nb_small_cells = 0;
        std::vector<swCell *> small_cells;
        for (swCell *cell : list_cells)
        {
            if (cell->getArea() < a_min)
            {
                nb_small_cells++;
                small_cells.push_back(cell);
            }
        }

        // Write small cells log file
        if (writeLogCell(small_cells, m_io->absolutePath(swFile::log_cell_small_file)) != EXIT_SUCCESS)
        {
            logError() << "  Error in writeLogCell()";
            logFatal(__FILE__, __func__, __LINE__);
        }

        // Feedback for user
        if (nb_small_cells > 0)
        {
            logWarn() << " " << nb_small_cells << "cells have their area smaller than" << a_min << "m2";
            logWarn() << "  Smaller cells listed in" << m_io->absolutePath(swFile::log_cell_small_file);
        }
    }

    return EXIT_SUCCESS;
}

std::vector<int> swMesh::getBcCodes()
/**
 * @brief return the list of bc codes
 */
{
    return list_bc_codes;
}

std::vector<std::vector<swEdgeBoundary *>> swMesh::getBoundaries()
/**
 * @brief return the list of boundaries
 */
{
    return list_boundaries;
}

swCell* swMesh::getCellFromCoordinates(double x, double y)
/**
 * @brief Returns the pointer to the cell containing the point with coordinates (x, y)
 * 
 * @param[in] x x-coordinate
 * @param[in] y y-coordinate
 * @return swCell* cell containing the point
 */
{
    swCell* found_cell = nullptr;

    for (swCell* cell : list_cells)
    {
        double theta = 0.;
        std::vector<swNode*> nodes = cell->getNodes();
        size_t nb_nodes = nodes.size();
        for (size_t p = 0; p < nb_nodes; ++p) // Parse the points of the cell
        { 
            double x1 = nodes[p]->getX();
            double y1 = nodes[p]->getY();
            double x2 = nodes[(p + 1)%nb_nodes]->getX(); // if j is past the last element, loop and take the first node
            double y2 = nodes[(p + 1)%nb_nodes]->getY();

            double cross_prod = (y2 - y) * (x1 - x) - (x2 - x) * (y1 - y);
            double dot_prod   = (x1 - x) * (x2 - x) + (y1 - y) * (y2 - y);

            double norm1 = std::sqrt((x1 - x) * (x1 - x) + (y1 - y) * (y1 - y));
            double norm2 = std::sqrt((x2 - x) * (x2 - x) + (y2 - y) * (y2 - y));
            norm1 = std::max(norm1, small_real);
            norm2 = std::max(norm2, small_real);

            double dtheta = (cross_prod / norm1 / norm2);
            double cosine = dot_prod / norm1 / norm2;
            if (cosine > 0.)
            {
                dtheta = std::asin(std::min(dtheta, 1.));
            }
            else
            {
                if (dtheta < 0.)
                {
                    dtheta = -4. * std::atan(1.) - std::asin(std::min(dtheta, 1.));
                }
                else
                {
                    dtheta = 4. * std::atan(1.) - std::asin(std::min(dtheta, 1.));
                }
            }
            theta = theta + dtheta;
        }

        if (std::fabs(theta) > 6.)
        {
            found_cell = cell;
            break;
        }
    }

    if (found_cell == nullptr)
    {
        logWarn() << __func__ << "Could not find cell for coordinates " << x << y;
    }

    return found_cell;
}

std::vector<swCell*> swMesh::getCellsList()
/**
 * @brief return the list of cells 
 */
{
    return list_cells;
}

std::vector<swEdgeBoundary*> swMesh::getEdgeBoundaryList()
/**
 * @brief return the list of boundary edges
 */
{
    return list_edges_boundary;
}

std::vector<swEdgeIntern*> swMesh::getEdgeInternList()
/**
 * @brief return the list of internal edges 
 */
{
    return list_edges_intern;
}

swIO* swMesh::getIO()
/**
 * @brief return the io class associated to the mesh
 */
{
    return m_io;
}

std::vector<swNode*> swMesh::getNodesList()
/**
 * @brief return the list of nodes
 */
{
    return list_nodes;
}

std::vector<std::vector<std::string>> swMesh::getQSections()
/**
 * @brief return the list of qsections
 */
{
    return list_qsections;
}

std::vector<std::string> swMesh::getQSectionNames(void)
/**
 * @brief return the list of qsections names
 */
{
    return list_qsections_names;
}

std::vector<std::vector<int>> swMesh::getQSectionSigns(void)
/**
 * @brief return the list of qsections signs
 */
{
    return list_qsections_signs;
}

bool swMesh::importGeometry(std::string filename_2dm, std::string filename_bc, size_t default_bc)
/**
 * @brief Read 2dm and bc files and set attributes values
 * @param[in] filename_2dm path to the 2dm file
 * @param[in] filename_bc path to the bc file
 * @param[in] default_bc default bc code
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    std::vector<std::vector<size_t>> cell_nodes_idx;
    std::vector<size_t> cell_sms_id;
    if (read2dm(filename_2dm, cell_nodes_idx, cell_sms_id) != EXIT_SUCCESS)
    {
        logError() << "  Error while reading 2dm file";
        logFatal(__FILE__, __func__, __LINE__);
    }

    //-- Build cells
    logInfo() << "  Building cells";
    buildCells(cell_sms_id, cell_nodes_idx);
    logInfo() << "  ...done";
    //== Build cells

    //-- Build edges
    logInfo() << "  Building edges";
    buildEdges(default_bc);
    logInfo() << "  ...done";
    //== Build edges


    if ((filename_bc != "") && (assignBC(filename_bc, default_bc) != EXIT_SUCCESS))
    {
        logError() << "  Error while reading bc file";
        logFatal(__FILE__, __func__, __LINE__);
    }

    return EXIT_SUCCESS;
}

bool swMesh::readQSections()
/**
 * @brief read the qsections input file and fill in the qsection attributes
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    // Open qsections file
    std::string filename = swFile::input_qsection_file;
    std::ifstream file(m_io->absolutePath(filename));
    if (file.is_open() == false)
    {
        logError() << "Unable to open file" << filename;
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }

    // Read number of lines
    int line_nb = 1;
    bool conversion_success = true;

    std::string line_txt = swIO::readLine(file);
    const size_t nb_line = swString::stringToUInt(line_txt, conversion_success);
    if (!conversion_success)
    {
        logError() << "line is misformatted: expecting a positive integer (number of lines to read)";
        swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
    }

    std::vector<std::vector<size_t>> read_node_pairs;
    std::vector<std::string> read_names;

    for (size_t k = 0; k < nb_line; ++k)
    {
        ++line_nb;
        //-- Read line
        line_txt = swIO::readLine(file);
        const std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());  
        if (line_split.size() != 3)
        {
            logError() << "line is misformatted : expecting 2 positive integers (nodes id) and a string (qsection name)";
            swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
        }

        size_t node_1 = swString::stringToUInt(line_split[0], conversion_success);
        size_t node_2 = swString::stringToUInt(line_split[1], conversion_success);
        std::string name = line_split[2];
        if (!conversion_success)
        {
            logError() << "line is misformatted : expecting 2 positive integers (nodes id) and a string (qsection name)";
            swIO::errorReadingLine(file, filename, line_txt, line_nb, __LINE__, __func__, __FILE__);
        }
        //== Read line

        // Fill lists with read elements
        read_node_pairs.push_back({node_1, node_2});
        read_names.push_back(name);

        // Fill the names list attribute
        if (std::find(list_qsections_names.begin(), list_qsections_names.end(), name) == list_qsections_names.end())
        {
            list_qsections_names.push_back(name);
        }
    }
    file.close();

    size_t nb_qsections = list_qsections_names.size();

    //-- Fill qsections list
    for (size_t k = 0; k < nb_qsections; k++)
    {
        std::vector<std::string> qsection;
        std::vector<int> qsection_signs;

        for (size_t l = 0; l < nb_line; l++)
        {
            if (read_names[l] == list_qsections_names[k])
            {
                
                bool found = false;
                for(swEdgeBoundary *edge_bnd : list_edges_boundary)
                {
                    std::vector<swNode *> nodes = edge_bnd->getNodes();
                    if ((nodes[0]->getIdx() == read_node_pairs[l][0] - 1) && 
                        (nodes[1]->getIdx() == read_node_pairs[l][1] - 1))
                    {
                        qsection.push_back("bnd_" + std::to_string(edge_bnd->getIdx()));
                        qsection_signs.push_back(1);
                        found = true;
                        break;
                    }
                    else if ((nodes[0]->getIdx() == read_node_pairs[l][1] - 1) && 
                             (nodes[1]->getIdx() == read_node_pairs[l][0] - 1))
                    {
                        qsection.push_back("bnd_" + std::to_string(edge_bnd->getIdx()));
                        qsection_signs.push_back(-1);
                        found = true;
                        break;
                    }
                }
                for(swEdgeIntern *edge_int : list_edges_intern)
                {
                    std::vector<swNode *> nodes = edge_int->getNodes();
                    if ((nodes[0]->getIdx() == read_node_pairs[l][0] - 1) && 
                        (nodes[1]->getIdx() == read_node_pairs[l][1] - 1))
                    {
                        qsection.push_back("int_" + std::to_string(edge_int->getIdx()));
                        qsection_signs.push_back(1);
                        found = true;
                        break;
                    }
                    else if ((nodes[0]->getIdx() == read_node_pairs[l][1] - 1) && 
                             (nodes[1]->getIdx() == read_node_pairs[l][0] - 1))
                    {
                        qsection.push_back("int_" + std::to_string(edge_int->getIdx()));
                        qsection_signs.push_back(-1);
                        found = true;
                        break;
                    }
                }
                if (!found)
                {
                    logError() << "No edge found with nodes" << read_node_pairs[k][0] << "and" << read_node_pairs[k][1];
                    logFatal(__FILE__, __func__, __LINE__);
                }
            }
        }
        list_qsections.push_back(qsection);
        list_qsections_signs.push_back(qsection_signs);
    }
    //== Fill qsections list

    return EXIT_SUCCESS;
}

bool swMesh::read2dm(std::string filename_2dm,
                     std::vector<std::vector<size_t>> &cell_nodes_idx,
                     std::vector<size_t> &cell_sms_id)
/**
 * @brief Read 2dm file and set geometry attributes
 * @param[in] filename_2dm name of the 2dm file
 * @param[in,out] cell_nodes_idx list of the nodes defining a cell
 * @param[in,out] cell_sms_id list of the cell id in SMS
 * @return bool EXIT_SUCCES/EXIT_FAILURE
 */
{
    // Open file
    std::ifstream file;
    file.open(m_io->absolutePath(filename_2dm));

    if (file.is_open() == false)
    {
        logError() << "  Could not open file" << m_io->absolutePath(filename_2dm);
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }

    //-- Read file (and build nodes)
    logInfo() << "  Reading 2dm file";
    std::string line_txt = swIO::readLine(file);
    size_t line_nb = 1;


    size_t node_idx = 0;
    while (!file.eof())
    {
        bool conversion_success = true;
        const std::vector<std::string> fields = swString::splitStr(line_txt, swString::regexWhiteSpace()); 

        if (fields[0] == "E3T")
        {
            size_t id = swString::stringToUInt(fields[1], conversion_success);
            std::vector<size_t> nodes(3);
            nodes[0] = swString::stringToUInt(fields[2], conversion_success);
            nodes[1] = swString::stringToUInt(fields[3], conversion_success);
            nodes[2] = swString::stringToUInt(fields[4], conversion_success);
            cell_sms_id.push_back(id);
            cell_nodes_idx.push_back(nodes);
        }
        else if (fields[0] == "E4Q")
        {
            size_t id = swString::stringToUInt(fields[1], conversion_success);
            std::vector<size_t> nodes(4);
            nodes[0] = swString::stringToUInt(fields[2], conversion_success);
            nodes[1] = swString::stringToUInt(fields[3], conversion_success);
            nodes[2] = swString::stringToUInt(fields[4], conversion_success);
            nodes[3] = swString::stringToUInt(fields[5], conversion_success);
            cell_sms_id.push_back(id);
            cell_nodes_idx.push_back(nodes);
        }
        else if (fields[0] == "ND")
        {
            size_t id = swString::stringToUInt(fields[1], conversion_success);
            double x = swString::stringToDouble(fields[2], conversion_success);
            double y = swString::stringToDouble(fields[3], conversion_success);
            double z = swString::stringToDouble(fields[4], conversion_success);
            swNode* node = new swNode(id, node_idx, x, y, z);
            node_idx++;
            list_nodes.push_back(node);
        }
        else if(fields[0] != "MESH2D")
        {
            logWarn() << "    Type not recognized while parsing 2dm file at line" << line_nb << ":" << line_txt.c_str();
        }

        if (!conversion_success)
        {
            logError() << "    Line" << line_nb << "is misformatted: format conversion failed";
            logFatal(__FILE__, __func__, __LINE__);
        }

        line_txt = swIO::readLine(file);
        line_nb++;
    }

    file.close();
    logInfo() << "  ...done";
    //== Read file (and build nodes)

    return EXIT_SUCCESS;
}

bool swMesh::topographyCorrection(const std::string& topo_filename,const std::string &param_corr)
/**
 * @brief read file containing topography data and compute the new topography of the cells
 * 
 * @param[in] topo_filename name of the file containing topography correction data
 * @param[in] param_corr value of the topography correction parameter
 * @return bool EXIT_SUCCES/EXIT_FAILURE
 */
{
    //-- Read file
    std::string path = m_io->absolutePath(topo_filename);
    std::ifstream file(path);
    if (file.is_open() == false)
    {
        logError() << "  Unable to open the topography correction file:" << path.c_str();
        logFatal(__FILE__, __func__, __LINE__);
    }
    logInfo() << "  Importing topography correction data from file:" << path.c_str();

    size_t line_nb = 0;

    // Skip the header
    std::string line_txt = swIO::readLine(file);
    ++line_nb;

    std::vector<std::vector<double>> list_cell_corr(list_cells.size());

    while (!file.eof())
    {
        line_txt = swIO::readLine(file);
        ++line_nb;
        if (!(line_txt.empty() || (line_txt[0] == '#') || (line_txt[0]== '=')))
        {
            //-- Read line
            std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
            if (line_split.size() != 3)
            {
                logError() << "    Line" << line_nb << "is misformatted:" << line_txt;
                logError() << "      expecting 3 numbers : x, y, z";
                logFatal(__FILE__, __func__, __LINE__);
            }
            bool conversion_success = true;
            double x = swString::stringToDouble(line_split[0], conversion_success);
            double y = swString::stringToDouble(line_split[1], conversion_success);
            double z = swString::stringToDouble(line_split[2], conversion_success);
            if (conversion_success == false)
            {
                logError() << "    Line" << line_nb << "is misformatted: format conversion";
                logFatal(__FILE__, __func__, __LINE__);
            }
            //== Read line

            // Search for cell containing point (x,y)
            swCell* cell = getCellFromCoordinates(x, y);
            if (cell != nullptr)
            {
                // Add zb to zb list of this cell
                size_t id_cell = cell->getIdx();
                list_cell_corr[id_cell].push_back(z);
            }
        }
    }

    logInfo() << "  ...done";
    // == read file

    //-- Open log files
    // for corrected topography (only the modified cells)
    std::string filename_corrected_cells = m_io->absolutePath(swFile::log_topo_corrected_cells_file);
    std::ofstream file_corrected(filename_corrected_cells);
    if (!file_corrected)
    {
        logError() << "  Unable to open topography corrected log file:" << filename_corrected_cells;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // for topography before correction (all cells)
    std::string filename_before_correction = m_io->absolutePath(swFile::log_topo_before_correction_file);
    std::ofstream file_before(filename_before_correction);
    if (!file_before)
    {
        logError() << "  Unable to open topography before correction log file:" << filename_before_correction;
        logFatal(__FILE__, __func__, __LINE__);
    }

    // for topography after correction (all cells)
    std::string filename_after_correction = m_io->absolutePath(swFile::log_topo_after_correction_file);
    std::ofstream file_after(filename_after_correction);
    if (!file_after)
    {
        logError() << "  Unable to open topography after correction log file:" << filename_after_correction;
        logFatal(__FILE__, __func__, __LINE__);
    }
    //== Open log files

    file_after     << "x\ty\tzb\n";
    file_before    << "x\ty\tzb\n";
    file_corrected << "x\ty\tzb\n";

    for (swCell* cell : list_cells)
    {
        // cell coordinates
        double x = cell->getXG();
        double y = cell->getYG();

        // write bottom elevation before correction
        double zb = cell->getZBottom();
        file_before << swString::numberToStr(x, 'e', 10) << "\t"
                    << swString::numberToStr(y, 'e', 10) << "\t"
                    << swString::numberToStr(zb, 'e', 10) << "\n";

        if (list_cell_corr[cell->getIdx()].size() > 0)
        {
            // correct bottom elevation according to param_corr option
            cell->correctTopography(list_cell_corr[cell->getIdx()], param_corr);

            // write bottom elevation after topography is corrected
            file_corrected << swString::numberToStr(x, 'e', 10) << "\t"
                           << swString::numberToStr(y, 'e', 10) << "\t"
                           << swString::numberToStr(cell->getZBottom(), 'e', 10) << "\n";
        }

        // write bottom elevation after correction
        file_after << swString::numberToStr(x, 'e', 10) << "\t"
                   << swString::numberToStr(y, 'e', 10) << "\t"
                   << swString::numberToStr(cell->getZBottom(), 'e', 10) << "\n";
    }

    file_after.close();
    file_before.close();
    file_corrected.close();

    return EXIT_SUCCESS;
}

bool swMesh::writeGeometryFile(const std::string& filename_geo,
                               const std::string& filename_input_conv,
                               const std::string& filename_2dm,
                               const std::string& filename_bc)
/**
 * @brief write the mesh into the binary file geometry.sw2d
 * @param[in] filename_geo name of the geometry file
 * @param[in] filename_input_conv name of the the converter input file
 * @param[in] filename_2dm name of the the 2dm file
 * @param[in] filename_bc name of the the bc file
 * @return bool EXIT_SUCCES/EXIT_FAILURE
 */
{
    // check that the input folder exists
    if (std::filesystem::exists(m_io->absolutePath(swFile::input_folder)) == false)
    {
        std::filesystem::create_directory(m_io->absolutePath(swFile::input_folder));
    }

    // Open file
    std::ofstream file(m_io->absolutePath(filename_geo), std::ios::binary);
    if (file.is_open() == false)
        {
            logError() << "    Unable to open file" << filename_geo.c_str();
            logFatal(__FILE__, __func__, __LINE__);
            return EXIT_FAILURE;
        }

    //-- write version number
    size_t version_lenght = m_version.size();
    file.write(reinterpret_cast<const char*>(&version_lenght), sizeof(version_lenght));
    file.write(m_version.c_str(), version_lenght);
    //== write version number
    
    //-- write number of elements
    size_t nb_cells = list_cells.size();
    size_t nb_edges_boundary = list_edges_boundary.size();
    size_t nb_edges_intern = list_edges_intern.size();
    size_t nb_nodes = list_nodes.size();
    size_t nb_boundaries = list_boundaries.size();
    file.write(reinterpret_cast<const char*>(&nb_cells), sizeof(nb_cells));
    file.write(reinterpret_cast<const char*>(&nb_edges_boundary), sizeof(nb_edges_boundary));
    file.write(reinterpret_cast<const char*>(&nb_edges_intern), sizeof(nb_edges_intern));
    file.write(reinterpret_cast<const char*>(&nb_nodes), sizeof(nb_nodes));
    file.write(reinterpret_cast<const char*>(&nb_boundaries), sizeof(nb_boundaries));
    //== write number of elements

    //-- write cells
    for (swCell* cell : list_cells)
    {
        size_t id = cell->getIdx(), sms_id = cell->getSMSIdx();
        size_t nodes_size = cell->getNodes().size();
        size_t edges_bnd_size = cell->getEdgesBoundary().size();
        size_t edges_int_size = cell->getEdgesIntern().size();
        size_t neigh_size = cell->getNeigh().size();
        double area = cell->getArea(), perimeter = cell->getPerimeter();
        double x_g = cell->getXG(), y_g = cell->getYG(), zBottom = cell->getZBottom();

        file.write(reinterpret_cast<const char*>(&id), sizeof(id));
        file.write(reinterpret_cast<const char*>(&sms_id), sizeof(sms_id));
        file.write(reinterpret_cast<const char*>(&nodes_size), sizeof(nodes_size));
        file.write(reinterpret_cast<const char*>(&edges_bnd_size), sizeof(edges_bnd_size));
        file.write(reinterpret_cast<const char*>(&edges_int_size), sizeof(edges_int_size));
        file.write(reinterpret_cast<const char*>(&neigh_size), sizeof(neigh_size));
        file.write(reinterpret_cast<const char*>(&area), sizeof(area));
        file.write(reinterpret_cast<const char*>(&perimeter), sizeof(perimeter));
        file.write(reinterpret_cast<const char*>(&x_g), sizeof(x_g));
        file.write(reinterpret_cast<const char*>(&y_g), sizeof(y_g));
        file.write(reinterpret_cast<const char*>(&zBottom), sizeof(zBottom));
    }
    logInfo() << "  cells ...done";
    //== write cells

    //-- write nodes
    for (swNode* node : list_nodes)
    {
        size_t id = node->getIdx(), sms_id= node->getSMSIdx();
        size_t cells_size = node->getCells().size();
        size_t edges_size = node->getEdgesBoundary().size() + node->getEdgesIntern().size();
        double x = node->getX(), y = node->getY(), z = node->getZ();

        file.write(reinterpret_cast<const char*>(&id), sizeof(id));
        file.write(reinterpret_cast<const char*>(&sms_id), sizeof(sms_id));
        file.write(reinterpret_cast<const char*>(&x), sizeof(x));
        file.write(reinterpret_cast<const char*>(&y), sizeof(y));
        file.write(reinterpret_cast<const char*>(&z), sizeof(z));
        file.write(reinterpret_cast<const char*>(&cells_size), sizeof(cells_size));
        file.write(reinterpret_cast<const char*>(&edges_size), sizeof(edges_size));
    }
    
    logInfo() << "  nodes ...done";
    //== write nodes

    //-- write edges
    for (swEdgeBoundary* edge_bnd : list_edges_boundary)
    {
        size_t id = edge_bnd->getIdx();
        size_t node_1_id = edge_bnd->getNodes()[0]->getIdx(), node_2_id = edge_bnd->getNodes()[1]->getIdx();
        double x_m = edge_bnd->getXM(), y_m = edge_bnd->getYM(), z_m = edge_bnd->getZM();
        double lenght = edge_bnd->getLength(), d_norm = edge_bnd->getDnorm();
        double norm_x = edge_bnd->getNormX(), norm_y = edge_bnd->getNormY();
        std::string type = "boundary";
        int bc_code = edge_bnd->getBCCode();

        file.write(reinterpret_cast<const char*>(&id), sizeof(id));
        file.write(reinterpret_cast<const char*>(&node_1_id), sizeof(node_1_id));
        file.write(reinterpret_cast<const char*>(&node_2_id), sizeof(node_2_id));
        file.write(reinterpret_cast<const char*>(&x_m), sizeof(x_m));
        file.write(reinterpret_cast<const char*>(&y_m), sizeof(y_m));
        file.write(reinterpret_cast<const char*>(&z_m), sizeof(z_m));
        file.write(reinterpret_cast<const char*>(&d_norm), sizeof(d_norm));
        file.write(reinterpret_cast<const char*>(&lenght), sizeof(lenght));
        file.write(reinterpret_cast<const char*>(&norm_x), sizeof(norm_x));
        file.write(reinterpret_cast<const char*>(&norm_y), sizeof(norm_y));
        file.write(reinterpret_cast<const char*>(&bc_code), sizeof(bc_code));

        size_t id_cell = edge_bnd->getCell()->getIdx();
        file.write(reinterpret_cast<const char*>(&id_cell), sizeof(id_cell));
    }

    for (swEdgeIntern* edge_int : list_edges_intern)
    {
        size_t id = edge_int->getIdx();
        size_t node_1_id = edge_int->getNodes()[0]->getIdx(), node_2_id = edge_int->getNodes()[1]->getIdx();
        double x_m = edge_int->getXM(), y_m = edge_int->getYM(), z_m = edge_int->getZM();
        double lenght = edge_int->getLength(), d_norm = edge_int->getDnorm();
        double norm_x = edge_int->getNormX(), norm_y = edge_int->getNormY();
        std::string type = "intern";

        file.write(reinterpret_cast<const char*>(&id), sizeof(id));
        file.write(reinterpret_cast<const char*>(&node_1_id), sizeof(node_1_id));
        file.write(reinterpret_cast<const char*>(&node_2_id), sizeof(node_2_id));
        file.write(reinterpret_cast<const char*>(&x_m), sizeof(x_m));
        file.write(reinterpret_cast<const char*>(&y_m), sizeof(y_m));
        file.write(reinterpret_cast<const char*>(&z_m), sizeof(z_m));
        file.write(reinterpret_cast<const char*>(&d_norm), sizeof(d_norm));
        file.write(reinterpret_cast<const char*>(&lenght), sizeof(lenght));
        file.write(reinterpret_cast<const char*>(&norm_x), sizeof(norm_x));
        file.write(reinterpret_cast<const char*>(&norm_y), sizeof(norm_y));

        for (swCell* cell : edge_int->getCells())
        {
            size_t id_cell = cell->getIdx();
            file.write(reinterpret_cast<const char*>(&id_cell), sizeof(id_cell));
        }
    }
    logInfo() << "  edges ...done";
    //== write edges

    //-- write boundaries
    for (size_t k = 0; k < nb_boundaries; k++)
    {
        size_t bnd_size = list_boundaries[k].size();
        file.write(reinterpret_cast<const char*>(&bnd_size), sizeof(bnd_size));
        for (size_t i = 0; i < bnd_size; i++)
        {
            size_t edge_idx = list_boundaries[k][i]->getIdx();
            file.write(reinterpret_cast<const char*>(&edge_idx), sizeof(edge_idx));
        }
    }
    
    for (size_t k = 0; k < nb_boundaries; k++)
    {
        int code = list_bc_codes[k];
        file.write(reinterpret_cast<const char*>(&code), sizeof(code));
    }
    //== write boundaries

    //-- write cell link with other elements
    for (swCell* cell : list_cells)
    {
        for (swNode* node : cell->getNodes())
        {
            size_t node_idx = node->getIdx();
            file.write(reinterpret_cast<const char*>(&node_idx), sizeof(node_idx));
        }

        for (swEdgeBoundary* edge_bnd : cell->getEdgesBoundary())
        {
            size_t edge_bnd_idx = edge_bnd->getIdx();
            file.write(reinterpret_cast<const char*>(&edge_bnd_idx), sizeof(edge_bnd_idx));
        }

        for (swEdgeIntern* edge_int : cell->getEdgesIntern())
        {
            size_t edge_int_idx = edge_int->getIdx();
            file.write(reinterpret_cast<const char*>(&edge_int_idx), sizeof(edge_int_idx));
        }

        for (swCell* neigh_cell : cell->getNeigh())
        {
            size_t neigh_cell_idx = neigh_cell->getIdx();
            file.write(reinterpret_cast<const char*>(&neigh_cell_idx), sizeof(neigh_cell_idx));
        }
    }
    logInfo() << "  cells connectivities ...done";
    //== write cell link with other elements

    //-- write node link with other elements
    for (swNode* node : list_nodes)
    {
        for (swCell* cell : node->getCells())
        {
            size_t id_cell = cell->getIdx();
            file.write(reinterpret_cast<const char*>(&id_cell), sizeof(id_cell));
        }
    }
    logInfo() << "  nodes connectivities ...done";
    //== write node link with other elements

    //-- copy content of the input files
    unsigned long long size_of_files = 0;

    //---- input.conv
    std::ifstream file_input;
    file_input.open(m_io->absolutePath(filename_input_conv));
    if (file_input.is_open() == false)
    {
        logError() << "  Unable to open:" << m_io->absolutePath(filename_input_conv);
        logFatal(__FILE__, __func__, __LINE__);
    }

    // read the file
    std::vector<std::string> lines_list;
    int nb_lines = 0;
    while (!file_input.eof())
    {
        std::string line_txt = swIO::readLine(file_input);
        lines_list.push_back(line_txt);
        nb_lines ++;
    }

    // copy to geometry.sw2d
    file.write(reinterpret_cast<const char*>(&nb_lines), sizeof(nb_lines));
    size_of_files += sizeof(nb_lines);
    for (std::string line_txt : lines_list)
    {
        size_t line_lenght = line_txt.size();
        file.write(reinterpret_cast<const char*>(&line_lenght), sizeof(line_lenght));
        file.write(line_txt.c_str(), line_lenght);
        size_of_files += line_lenght + sizeof(line_lenght);
    }
    //==== input.conv

    //---- .2dm
    std::ifstream file_2dm;
    file_2dm.open(m_io->absolutePath(filename_2dm));
    if (file_2dm.is_open() == false)
    {
        logError() << "  Unable to open:" << m_io->absolutePath(filename_2dm);
        logFatal(__FILE__, __func__, __LINE__);
    }

    // read the file
    lines_list.clear();
    nb_lines = 0;
    while (!file_2dm.eof())
    {
        std::string line_txt = swIO::readLine(file_2dm);
        lines_list.push_back(line_txt);
        nb_lines ++;
    }

    // copy to geometry.sw2d
    file.write(reinterpret_cast<const char*>(&nb_lines), sizeof(nb_lines));
    size_of_files += sizeof(nb_lines);
    for (std::string line_txt : lines_list)
    {
        size_t line_lenght = line_txt.size();
        file.write(reinterpret_cast<const char*>(&line_lenght), sizeof(line_lenght));
        file.write(line_txt.c_str(), line_lenght);
        size_of_files += line_lenght + sizeof(line_lenght);
    }
    //==== .2dm

    //---- .bc
    if (filename_bc == "")
    {
        nb_lines = 0;
        file.write(reinterpret_cast<const char*>(&nb_lines), sizeof(nb_lines));
        size_of_files += sizeof(nb_lines);
    }
    else
    {
        std::ifstream file_bc;
        file_bc.open(m_io->absolutePath(filename_bc));
        if (file_bc.is_open() == false)
        {
            logError() << "  Unable to open:" << m_io->absolutePath(filename_bc);
            logFatal(__FILE__, __func__, __LINE__);
        }

        // read the file
        lines_list.clear();
        nb_lines = 0;
        while (!file_bc.eof())
        {
            std::string line_txt = swIO::readLine(file_bc);
            lines_list.push_back(line_txt);
            nb_lines ++;
        }

        // copy to geometry.sw2d
        file.write(reinterpret_cast<const char*>(&nb_lines), sizeof(nb_lines));
        size_of_files += sizeof(nb_lines);
        for (std::string line_txt : lines_list)
        {
            size_t line_lenght = line_txt.size();
            file.write(reinterpret_cast<const char*>(&line_lenght), sizeof(line_lenght));
            file.write(line_txt.c_str(), line_lenght);
            size_of_files += line_lenght + sizeof(line_lenght);
        }
    }    
    //==== .bc

    file.write(reinterpret_cast<const char*>(&size_of_files), sizeof(size_of_files));
    //== copy content of the input files

    file.flush();
	file.close();

    return EXIT_SUCCESS;
}

bool swMesh::writeLog(void)
/**
 * @brief write the mesh logs
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    if (writeLogCell(list_cells, m_io->absolutePath(swFile::log_cell_file)) != EXIT_SUCCESS)
    {
        logError() << "Error in writeLogCell()";
        return EXIT_FAILURE;
    }
    if (writeLogEdgeBoundary() != EXIT_SUCCESS)
    {
        logError() << "Error in writeLogEdgeBoundary()";
        return EXIT_FAILURE;
    }
    if (writeLogEdgeIntern() != EXIT_SUCCESS)
    {
        logError() << "Error in writeLogEdgeIntern()";
        return EXIT_FAILURE;
    }
    if (writeLogNode() != EXIT_SUCCESS)
    {
        logError() << "Error in writeLogNode()";
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}

bool swMesh::writeLogCell(std::vector<swCell *> cells, std::string filename)
/**
 * @brief Write a cells log file
 * 
 * @param[in] cells list of cells to write in the file
 * @param[in] filename name of the log file
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
    // Open file
    std::ofstream log_file(filename);
    if (log_file.is_open() == false)
    {
        logError() << "  Unable to open the cell log file:" << filename;
        return EXIT_FAILURE;
    }
    logInfo() << "  Writing cell log file:" << filename;

    // Write header
    log_file << "CellNum\tSMSNum\tXg\tYg\tZb\tA\tNbNodes\t"
                "NodeNum\tNodeNum\tNodeNum\tNodeNum\t"
                "EdgeNum\tEdgeNum\tEdgeNum\tEdgeNum\t"
                "NeighborCellNum\tNeighborCellNum\tNeighborCellNum\tNeighborCellNum\n";

    //-- Write cells
    for (swCell *cell : list_cells)
    {
        double x = cell->getXG();
        double y = cell->getYG();
        double zb = cell->getZBottom();
        double a = cell->getArea();
        double p = cell->getPerimeter();
        size_t id = cell->getIdx();
        size_t sms_number = cell->getSMSIdx();
        std::vector<swNode *> nodes = cell->getNodes();
        std::vector<swEdgeBoundary *> edges_bnd = cell->getEdgesBoundary();
        std::vector<swEdgeIntern *> edges_int = cell->getEdgesIntern();
        std::vector<swCell *> neigh_cells = cell->getNeigh();
        size_t nb_nodes = nodes.size();
        log_file << id << "\t"
                 << sms_number << "\t"
                 << swString::numberToStr(x, 'e', 10) << "\t"
                 << swString::numberToStr(y, 'e', 10)<< "\t"
                 << swString::numberToStr(zb, 'e', 10) << "\t"
                 << swString::numberToStr(a, 'e', 10) << "\t"
                 << swString::numberToStr(p, 'e', 10) << "\t"
                 << nb_nodes << "\t";
        for (swNode *node : nodes)
        {
            log_file << node->getSMSIdx() << "\t";
        }
        for (size_t p = nb_nodes; p < 4; ++p)
        {
            log_file << "-999\t";
        }
        for (swEdgeBoundary *edge_bnd : edges_bnd)
        {
            log_file << edge_bnd->getIdx() << "\t";
        }
        for (swEdgeIntern *edge_int : edges_int)
        {
            log_file << edge_int->getIdx() << "\t";
        }
        for (size_t p = nb_nodes; p < 4; ++p)
        {
            log_file << "-999\t";
        }
        for (swCell *neigh_cell : neigh_cells)
        {
            log_file << neigh_cell->getSMSIdx() << "\t";
        }
        for (size_t p = nb_nodes; p < 4; ++p)
        {
            log_file << "-999\t";
        }
        log_file << "\n";
    }
    //== Write cells

    log_file.close();
    logInfo() << "  ...done";

    return EXIT_SUCCESS;
}

bool swMesh::writeLogEdgeBoundary(void)
/**
 * @brief write the boundary edge log
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
     
    // Open file
    std::string filepath = m_io->absolutePath(swFile::log_edge_bnd_file);
    std::ofstream log_file(filepath);
    if (log_file.is_open() == false)
    {
        logError() << "    Unable to open the boundary edge log file:" << filepath;
        return EXIT_FAILURE;
    }
    logInfo() << "  Writing boundary edge log file:" << filepath;

    // Write header
    log_file << "EdgeNum\t"
             << "Xint\tYint\tZint\t"
             << "Length\tNodeBeg\tNodeEnd\t"
             << "NormX\tNormY\t"
             << "Cell\tBC\n";

    //-- Write boundary edges
    for (swEdgeBoundary *edge_bnd : list_edges_boundary)
    {
        double x = edge_bnd->getXM();
        double y = edge_bnd->getYM();
        double z = edge_bnd->getZM();
        double length = edge_bnd->getLength();
        std::vector<swNode *> edge_nodes = edge_bnd->getNodes();
        double norm_x = edge_bnd->getNormX();
        double norm_y = edge_bnd->getNormY();
        swCell *cell = edge_bnd->getCell();
        int bc = edge_bnd->getBCCode();

        log_file << edge_bnd->getIdx() << "\t"
                 << swString::numberToStr(x, 'e', 10) << "\t"
                 << swString::numberToStr(y, 'e', 10) << "\t"
                 << swString::numberToStr(z, 'e', 10) << "\t"
                 << swString::numberToStr(length, 'e', 10) << "\t"
                 << edge_nodes[0]->getSMSIdx() << "\t" << edge_nodes[1]->getSMSIdx() << "\t"
                 << swString::numberToStr(norm_x, 'e', 10) << "\t"
                 << swString::numberToStr(norm_y, 'e', 10) << "\t"
                 << cell->getIdx() << "\t"
                 << bc << "\n";
    }
    //== Write boundary edges

    log_file.close();
    logInfo() << "  ...done";

    return EXIT_SUCCESS;
}

bool swMesh::writeLogEdgeIntern(void)
/**
 * @brief write the internal edge log
 * @return bool EXIT_SUCCES/EXIT_FAILURE 
 */
{
     
    // Open file
    std::string filepath = m_io->absolutePath(swFile::log_edge_int_file);
    std::ofstream log_file(filepath);
    if (log_file.is_open() == false)
    {
        logError() << "    Unable to open the internal edge log file:" << filepath;
        return EXIT_FAILURE;
    }
    logInfo() << "  Writing internal edge log file:" << filepath;

    // Write header
    log_file << "EdgeNum\t"
             << "Xint\tYint\tZint\t"
             << "Length\tNodeBeg\tNodeEnd\t"
             << "NormX\tNormY\t"
             << "CelL\tCelR\n";

    //-- Write internal edges
    for (swEdgeIntern *edge_int : list_edges_intern)
    {
        double x = edge_int->getXM();
        double y = edge_int->getYM();
        double z = edge_int->getZM();
        double length = edge_int->getLength();
        std::vector<swNode *> edge_nodes = edge_int->getNodes();
        double norm_x = edge_int->getNormX();
        double norm_y = edge_int->getNormY();
        std::vector<swCell *> cells = edge_int->getCells();

        log_file << edge_int->getIdx() << "\t"
                 << swString::numberToStr(x, 'e', 10) << "\t"
                 << swString::numberToStr(y, 'e', 10) << "\t"
                 << swString::numberToStr(z, 'e', 10) << "\t"
                 << swString::numberToStr(length, 'e', 10) << "\t"
                 << edge_nodes[0]->getSMSIdx() << "\t" << edge_nodes[1]->getSMSIdx() << "\t"
                 << swString::numberToStr(norm_x, 'e', 10) << "\t"
                 << swString::numberToStr(norm_y, 'e', 10) << "\t"
                 << cells[0]->getIdx() << "\t"
                 << cells[1]->getIdx() << "\n";
    }
    //== Write internal edges

    log_file.close();
    logInfo() << "  ...done";

    return EXIT_SUCCESS;
}

bool swMesh::writeLogNode(void)
/**
 * @brief write the nodes log
 * @return bool EXIT_SUCCES/EXIT_FAILURE
 */
{
    // Open file
    std::string filepath = m_io->absolutePath(swFile::log_node_file);
    std::ofstream log_file(filepath);
    if (log_file.is_open() == false)
    {
        logError() << "    Unable to open the edge log file:" << filepath;
        return EXIT_FAILURE;
    }
    logInfo() << "  Writing node log file:" << filepath;

    // Write header
    log_file << "NodeNum\tSMSNum\tX\tY\tZ\t"
             << "NbCell\tNbEdge\t"
             << "CellsSMSId\t"
             << "ListOfEdges\n";

    //-- Write nodes
    for (swNode *node : list_nodes)
    {
        double x = node->getX();
        double y = node->getY();
        double z = node->getZ();
        std::vector<swCell *> cells = node->getCells();
        std::vector<swEdgeBoundary *> edges_bnd = node->getEdgesBoundary();
        std::vector<swEdgeIntern *> edges_int = node->getEdgesIntern();
        size_t nb_cells = cells.size();
        size_t nb_edges = edges_bnd.size() + edges_int.size();
        log_file << node->getIdx() << "\t" << node->getSMSIdx() << "\t"
                 << swString::numberToStr(x, 'e', 10) << "\t"
                 << swString::numberToStr(y, 'e', 10) << "\t"
                 << swString::numberToStr(z, 'e', 10) << "\t"
                 << nb_cells << "\t" << nb_edges;
        for (swCell *cell : cells)
        {
            log_file << "\t" << cell->getSMSIdx();
        }
        for (swEdgeBoundary *edge_bnd : edges_bnd)
        {
            log_file << "\t" << edge_bnd->getIdx();
        }
        for (swEdgeIntern *edge_int : edges_int)
        {
            log_file << "\t" << edge_int->getIdx();
        }
        log_file << "\n";
    }
    //== Write nodes

    log_file.close();
    logInfo() << "  ...done";

    return EXIT_SUCCESS;
}
