#include "swPorosityMapperPipeline.h"

#include <cmath>
#include <fstream>
#include <string>

#include "../common/geometry/swCell.h"
#include "../common/geometry/swEdgeBoundary.h"
#include "../common/geometry/swEdgeIntern.h"
#include "../common/swFile.h"
#include "../common/swIO.h"
#include "../common/swGeometry.h"
#include "../common/swLog.h"
#include "../common/geometry/swMesh.h"
#include "../common/geometry/swNode.h"
#include "../common/swString.h"

swPorosityMapperPipeline::swPorosityMapperPipeline(swIO *io, swPorosityMapperParam *param, swMesh *mesh)
/**
 * @brief Construct a new swPorosityMapperPipeline::swPorosityMapperPipeline object
 * @param[in] param pointer to the param object
 * @param[in] mesh pointer to the mesh
 * @param[in] io pointer to the Input/Output handler
 */
{
    m_io = io;
    m_mesh = mesh;
    m_param = param;
}

double swPorosityMapperPipeline::cellPorosity(swCell *cell, double h)
/**
 * @brief Compute the porosity of a cell according to a heigh and the buildings contained in the cell
 * 
 * @param[in] cell cell whose porosity is computed
 * @param[in] h heigh considered to compute the porosity
 * @return phi_omega : the cell porosity 
 */
{
    // build vectors of the cell coordinates
    std::vector<double> x_cell, y_cell;
    for (swNode *node : cell->getNodes())
    {
        x_cell.push_back(node->getX());
        y_cell.push_back(node->getY());
    }
    
    double buildings_area = 0;

    if (m_param->getInput() == "footprint")
    {
        size_t nb_buildings = x_buildings.size();
        for (size_t k = 0; k < nb_buildings; k++)
        {
            if((m_param->getModel() != "ddp") || (h_buildings[k] > h))
            {
                std::vector<double> x_building = x_buildings[k];
                std::vector<double> y_building = y_buildings[k];

                buildings_area += swGeometry::computeIntersectionArea(x_cell, y_cell, x_building, y_building);
            }
        }
    }
    else if (m_param->getInput() == "raster")
    {
        // find the cell bounding box
        std::vector<size_t> bb = m_buildings->getBoundingBox(cell);

        for (size_t i=bb[0]; i<=bb[2]; i++)
        {
            for (size_t j=bb[1]; j<=bb[3]; j++)
            {
                // For each point of the bouding box, check if the point is a building
                if (m_buildings->getBuildings()[j][i] > h)
                {
                    // build vector of the building point coordinates
                    double dx = m_buildings->getDeltaX();
                    double dy = m_buildings->getDeltaY();
                    std::vector<double> x_buildings = {m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2, m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2};
                    std::vector<double> y_buildings = {m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] + dy/2, m_buildings->getY()[j] + dy/2};

                    // compute the area of intersection between the building point and the cell
                    buildings_area += swGeometry::computeIntersectionArea(x_cell, y_cell, x_buildings, y_buildings);
                }
            }
        }
    }

    return(1 - buildings_area/cell->getArea());
}

bool swPorosityMapperPipeline::computeDDP()
/**
 * @brief Compute the porosities (with ddp model) of each cell and edge of the mesh and store the values is **_ddp_array
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    const std::vector<swCell *> cells = m_mesh->getCellsList();
    const std::vector<swEdgeBoundary *> edges_bnd = m_mesh->getEdgeBoundaryList();
    const std::vector<swEdgeIntern *> edges_int = m_mesh->getEdgeInternList();

    //-- Compute phi_omega
    logInfo() << "    computing phi_omega";
    for (swCell *cell : cells)
    {
        // compute the porosity for each heigh
        std::vector<double> cell_porosity_line;
        for (size_t idx_h=0; idx_h<topo_list.size(); idx_h++)
        {
            double phi_omega = cellPorosity(cell, topo_list[idx_h]);
            cell_porosity_line.push_back(phi_omega);
        }

        cell_porosity_array.push_back(cell_porosity_line);
        cell_porosity_line.clear();
    }
    logInfo() << "    ...done";
    //== Compute phi_omega

    //-- Compute phi_gamma
    logInfo() << "    computing phi_gamma";
    for (swEdgeBoundary *edge_bnd : edges_bnd)
    {
        // compute the porosity for each heigh
        std::vector<double> edge_porosity_line;
        for (size_t idx_h=0; idx_h<topo_list.size(); idx_h++)
        {
            double phi_gamma = edgePorosity(edge_bnd, topo_list[idx_h]);
            edge_porosity_line.push_back(phi_gamma);

            // check phiG <= phiW
            swCell *cell = edge_bnd->getCell();
            if (phi_gamma > cell_porosity_array[cell->getIdx()][idx_h])
            {
                logWarn() << "Boundary edge" << edge_bnd->getIdx() << "has greater porosity (" << phi_gamma << ") than its cell" << cell->getIdx() << "(" << cell_porosity_array[cell->getIdx()][idx_h] << ") for topography h = " << topo_list[idx_h];
            }   
        }
        edge_bnd_porosity_array.push_back(edge_porosity_line);
        edge_porosity_line.clear();
    }
    logInfo() << "    ...done";
    
    for (swEdgeIntern *edge_int : edges_int)
    {
        // compute the porosity for each heigh
        std::vector<double> edge_porosity_line;
        for (size_t idx_h=0; idx_h<topo_list.size(); idx_h++)
        {
            double phi_gamma = edgePorosity(edge_int, topo_list[idx_h]);
            edge_porosity_line.push_back(phi_gamma);

            // check phiG <= phiW
            swCell *cellL = edge_int->getCells()[0];
            if (phi_gamma > cell_porosity_array[cellL->getIdx()][idx_h])
            {
                logWarn() << "Internal edge" << edge_int->getIdx() << "has greater porosity (" << phi_gamma << ") than its left cell" << cellL->getIdx() << "(" << cell_porosity_array[cellL->getIdx()][idx_h] << ") for topography h = " << topo_list[idx_h];
            }
            swCell *cellR = edge_int->getCells()[1];
            if (phi_gamma > cell_porosity_array[cellR->getIdx()][idx_h])
            {
                logWarn() << "Internal edge" << edge_int->getIdx() << "has greater porosity (" << phi_gamma << ") than its right cell" << cellR->getIdx() << "(" << cell_porosity_array[cellR->getIdx()][idx_h] << ") for topography h = " << topo_list[idx_h];
            }
        }
        edge_int_porosity_array.push_back(edge_porosity_line);
        edge_porosity_line.clear();
    }
    logInfo() << "    ...done";
    //== Compute phi_gamma

    return EXIT_SUCCESS;
}

bool swPorosityMapperPipeline::computeDIP()
/**
 * @brief Compute the porosity (with dip model) of each cell and edge of the mesh and store the values is **_dip_array
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    const std::vector<swCell *> cells = m_mesh->getCellsList();
    const std::vector<swEdgeBoundary *> edges_bnd = m_mesh->getEdgeBoundaryList();
    const std::vector<swEdgeIntern *> edges_int = m_mesh->getEdgeInternList();

    //-- Compute phi_omega
    logInfo() << "    computing phi_omega";
    for (swCell *cell : cells)
    {
        // compute the porosity
        double phi_omega = cellPorosity(cell, 0);
        cell_porosity_list.push_back(phi_omega);
    }
    logInfo() << "    ...done";
    //== Compute phi_omega

    //-- Compute phi_gamma
    logInfo() << "    computing phi_gamma";
    for (swEdgeBoundary *edge_bnd : edges_bnd)
    {
        // compute the porosity
        double phi_gamma = edgePorosity(edge_bnd, 0);
        edge_bnd_porosity_list.push_back(phi_gamma);

        // check phiG <= phiW
        swCell *cell = edge_bnd->getCell();
        if (phi_gamma > cell_porosity_list[cell->getIdx()])
        {
            logWarn() << "Boundary edge" << edge_bnd->getIdx() << "has greater porosity (" << phi_gamma << ") than its cell" << cell->getIdx() << "(" << cell_porosity_list[cell->getIdx()] << ")";
        }        
    }
    logInfo() << "    ...done";
    for (swEdgeIntern *edge_int : edges_int)
    {
        // compute the porosity
        double phi_gamma = edgePorosity(edge_int, 0);
        edge_int_porosity_list.push_back(phi_gamma);

        // check phiG <= phiW
        swCell *cellL = edge_int->getCells()[0];
        if (phi_gamma > cell_porosity_list[cellL->getIdx()])
        {
            logWarn() << "Internal edge" << edge_int->getIdx() << "has greater porosity (" << phi_gamma << ") than its left cell" << cellL->getIdx() << "(" << cell_porosity_list[cellL->getIdx()] << ")";
        }
        swCell *cellR = edge_int->getCells()[1];
        if (phi_gamma > cell_porosity_list[cellR->getIdx()])
        {
            logWarn() << "Internal edge" << edge_int->getIdx() << "has greater porosity (" << phi_gamma << ") than its right cell" << cellR->getIdx() << "(" << cell_porosity_list[cellR->getIdx()] << ")";
        }
    }
    logInfo() << "    ...done";
    //== Compute phi_gamma

    return EXIT_SUCCESS;
}

bool swPorosityMapperPipeline::computeSP()
/**
 * @brief Compute the porosity (with sp model) of each cell of the mesh and store the values is **_sp_array
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    const std::vector<swCell *> cells = m_mesh->getCellsList();

    //-- Compute phi_omega
    logInfo() << "    computing phi_omega";
    for (swCell *cell : cells)
    {
        // compute the porosity
        double phi_omega = cellPorosity(cell, 0);
        cell_porosity_list.push_back(phi_omega);
    }
    logInfo() << "    ...done";
    //== Compute phi_omega

    return EXIT_SUCCESS;
}

double swPorosityMapperPipeline::edgePorosity(swEdgeBoundary *edge, double h)
/**
 * @brief Compute the porosity of a edge according to a heigh and the buildings contained on the edge
 * 
 * @param[in] edge edge whose porosity is computed
 * @param[in] h heigh considered to compute the porosity
 * @return phi_omega : the edge porosity 
 */
{
    // build vectors of the edge coordinates
    std::vector<double> x_edge = {edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX()};
    std::vector<double> y_edge = {edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY()};

    double buildings_lenght = 0;

    if(m_param->getInput() == "footprint")
    {
        size_t nb_buildings = x_buildings.size();
        for (size_t k = 0; k < nb_buildings; k++)
        {
            if((m_param->getModel() != "ddp") || (h_buildings[k] > h))
            {
                std::vector<double> x_building = x_buildings[k];
                std::vector<double> y_building = y_buildings[k];

                buildings_lenght += swGeometry::lenghtSegmentInPolygon(x_building, y_building, x_edge, y_edge);
            }
        }
    }
    else if(m_param->getInput() == "raster")
    {
        // find the edge bounding box
        std::vector<size_t> bb = m_buildings->getBoundingBox(edge);

        // Check if the edge is on the side of a building point to avoid counting twice the lenght
        bool horizontal = (edge->getNodes()[0]->getY() == edge->getNodes()[1]->getY()) && (edge->getNodes()[0]->getY() == m_buildings->getY()[bb[3]]-m_buildings->getDeltaY()/2) && bb[1]!=bb[3];
        bool vertical = (edge->getNodes()[0]->getX() == edge->getNodes()[1]->getX()) && (edge->getNodes()[0]->getX() == m_buildings->getX()[bb[2]]-m_buildings->getDeltaX()/2) && bb[0]!=bb[2];

        for (size_t i=bb[0]; i<=bb[2]; i++)
        {
            for (size_t j=bb[1]; j<=bb[3]; j++)
            {
                // For each point of the bouding box, check if the point is a building
                if (m_buildings->getBuildings()[j][i] > h)
                {                
                    // build vectors of the building point coordinates
                    double dx = m_buildings->getDeltaX();
                    double dy = m_buildings->getDeltaY();
                    std::vector<double> x_buildings = {m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2, m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2};
                    std::vector<double> y_buildings = {m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] + dy/2, m_buildings->getY()[j] + dy/2};

                    // compute the lenght of intersection between the building point ad the edge
                    double int_lenght = swGeometry::lenghtSegmentInPolygon(x_buildings, y_buildings, x_edge, y_edge);
                    if ((horizontal && (m_buildings->getBuildings()[bb[1]][i] > h) && (m_buildings->getBuildings()[bb[3]][i] > h)) ||
                        (vertical && (m_buildings->getBuildings()[j][bb[0]] > h) && (m_buildings->getBuildings()[j][bb[2]] > h)))
                    {
                        int_lenght /= 2;                    
                    }
                    buildings_lenght += int_lenght;
                }
            }
        }
    }
    
    return(1 - buildings_lenght/edge->getLength());
}

double swPorosityMapperPipeline::edgePorosity(swEdgeIntern *edge, double h)
/**
 * @brief Compute the porosity of a edge according to a heigh and the buildings contained on the edge
 * 
 * @param[in] edge edge whose porosity is computed
 * @param[in] h heigh considered to compute the porosity
 * @return phi_omega : the edge porosity 
 */
{
    // build vectors of the edge coordinates
    std::vector<double> x_edge = {edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX()};
    std::vector<double> y_edge = {edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY()};

    double buildings_lenght = 0;

    if(m_param->getInput() == "footprint")
    {
        size_t nb_buildings = x_buildings.size();
        for (size_t k = 0; k < nb_buildings; k++)
        {
            if((m_param->getModel() != "ddp") || (h_buildings[k] > h))
            {
                std::vector<double> x_building = x_buildings[k];
                std::vector<double> y_building = y_buildings[k];

                buildings_lenght += swGeometry::lenghtSegmentInPolygon(x_building, y_building, x_edge, y_edge);
            }
        }
    }
    else if(m_param->getInput() == "raster")
    {
        // find the edge bounding box
        std::vector<size_t> bb = m_buildings->getBoundingBox(edge);

        // Check if the edge is on the side of a building point to avoid counting twice the lenght
        bool horizontal = (edge->getNodes()[0]->getY() == edge->getNodes()[1]->getY()) && (edge->getNodes()[0]->getY() == m_buildings->getY()[bb[3]]-m_buildings->getDeltaY()/2) && bb[1]!=bb[3];
        bool vertical = (edge->getNodes()[0]->getX() == edge->getNodes()[1]->getX()) && (edge->getNodes()[0]->getX() == m_buildings->getX()[bb[2]]-m_buildings->getDeltaX()/2) && bb[0]!=bb[2];

        for (size_t i=bb[0]; i<=bb[2]; i++)
        {
            for (size_t j=bb[1]; j<=bb[3]; j++)
            {
                // For each point of the bouding box, check if the point is a building
                if (m_buildings->getBuildings()[j][i] > h)
                {                
                    // build vectors of the building point coordinates
                    double dx = m_buildings->getDeltaX();
                    double dy = m_buildings->getDeltaY();
                    std::vector<double> x_buildings = {m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2, m_buildings->getX()[i] - dx/2, m_buildings->getX()[i] + dx/2};
                    std::vector<double> y_buildings = {m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] - dy/2, m_buildings->getY()[j] + dy/2, m_buildings->getY()[j] + dy/2};

                    // compute the lenght of intersection between the building point ad the edge
                    double int_lenght = swGeometry::lenghtSegmentInPolygon(x_buildings, y_buildings, x_edge, y_edge);
                    if ((horizontal && (m_buildings->getBuildings()[bb[1]][i] > h) && (m_buildings->getBuildings()[bb[3]][i] > h)) ||
                        (vertical && (m_buildings->getBuildings()[j][bb[0]] > h) && (m_buildings->getBuildings()[j][bb[2]] > h)))
                    {
                        int_lenght /= 2;                    
                    }
                    buildings_lenght += int_lenght;
                }
            }
        }
    }
    
    return(1 - buildings_lenght/edge->getLength());
}

void swPorosityMapperPipeline::init(void)
/**
 * @brief initialize the pipeline
 */
{
    //-- Read the buildings
    std::string filename = m_io->absolutePath(swFile::input_porosity_mapper_buildings_file);
    logInfo() << "Read the buildings position from :" << filename;

    if (m_param->getInput() == "raster")
    {
        m_buildings = new swRasterBuildings();
        m_buildings->build(filename);
        topo_list = m_buildings->getTopoList();
    }
    else if (m_param->getInput() == "footprint")
    {
        topo_list.push_back(0);

        std::ifstream file(filename);
        if (file.is_open() == false)
        {
            logError() << "    Unable to open file :" << filename;
            logFatal(__FILE__, __func__, __LINE__);
        }

        // read the header
        std::string line_txt = swIO::readLine(file);

        // read coordinates
        while(!file.eof())
        {
            std::vector<double> x_line;
            std::vector<double> y_line;

            line_txt = swIO::readLine(file);
            std::vector<std::string> line_split = swString::splitStr(line_txt, swString::regexWhiteSpace());
            bool conversion_success = true;

            size_t nb_points;
            double last_elt = swString::stringToDouble(line_split.back(), conversion_success);
            if (last_elt < 0)
            {
                nb_points = (line_split.size()-1)/2;
                h_buildings.push_back(-last_elt);
                if (std::find(topo_list.begin(), topo_list.end(), -last_elt) == topo_list.end()) // add the building heigh to topo_list if it's not already in
                {
                    topo_list.push_back(-last_elt);
                }
            }
            else
            {
                nb_points = line_split.size()/2;
            }

            for (size_t k = 0; k < nb_points*2; k += 2)
            {
                double x_k = swString::stringToDouble(line_split[k], conversion_success);
                double y_k = swString::stringToDouble(line_split[k+1], conversion_success);
                x_line.push_back(x_k);
                y_line.push_back(y_k);
            }

            x_buildings.push_back(x_line);
            y_buildings.push_back(y_line);
        }

        file.close();

        // Sort topo_list in ascending order
        std::sort(topo_list.begin(), topo_list.end());
    }

    logInfo() << "...done";
    //== Read the buildings
}

bool swPorosityMapperPipeline::mapPorosity(std::string model)
/**
 * @brief compute the porosities according to the model over the mesh and write map files
 * @param[in] model the porosity model
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    if (model == "ddp")
    {
        logInfo() << "Compute DDP";
        if (computeDDP() != EXIT_SUCCESS)
        {
            logError() << "Error in computeDDP()";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";

        logInfo() << "Write porosity maps";
        if (writeDdp(-2, topo_list.size()) != EXIT_SUCCESS)
        {
            logError() << "Error while writing DDP maps";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";
    }
    else if (model == "dip")
    {
        logInfo() << "Compute DIP";
        if (computeDIP() != EXIT_SUCCESS)
        {
            logError() << "Error in computeDIP()";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";

        logInfo() << "Write porosity maps";
        if (writeDip() != EXIT_SUCCESS)
        {
            logError() << "Error while writing DIP maps";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";
    }
    else if (model == "sp")
    {
        logInfo() << "Compute SP";
        if (computeSP() != EXIT_SUCCESS)
        {
            logError() << "Error in computeSP()";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";

        logInfo() << "Write porosity map";
        if (writeSp() != EXIT_SUCCESS)
        {
            logError() << "Error while writing SP map";
            logFatal(__FILE__, __func__, __LINE__);
        }
        logInfo() << "...done";
    }
    else
    {
        logError() << "unknown model:" << model;
        logFatal(__FILE__, __func__, __LINE__);
    }

    return EXIT_SUCCESS;
}

bool swPorosityMapperPipeline::writeDdp(int law_type, int n_tab_max)
/**
 * @brief write the DDP porosity maps
 * 
 * @param[in] law_type the law type
 * @param[in] n_tab_max the number of vertical discretization levels to use for the law
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    size_t n_topo = topo_list.size();

    //-- Cells
    std::ofstream file;
    std::string filename = m_io->absolutePath(swFile::input_ddp_porosity_cell_file);
    file.open(filename);
    if (file.is_open() == false)
    {
        logError() << "    Unable to open" << (filename).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing cell porosity map :" << swFile::input_ddp_porosity_cell_file;

    file << "Unif\t0" << std::endl;
    file << "NtabMax\t" << n_tab_max << std::endl;
    file << "==== Default Param" << std::endl;
    file << "Type\tParameters" << std::endl;
    file << law_type << "\t" << n_topo << "\t";
    if ((law_type == -2) || (law_type == -1))
    {
        for (size_t k = 0; k < n_topo-1; k++)
        {
            file << cell_porosity_array[0][k] << "\t" << topo_list[k] << "\t";
        }
        file << cell_porosity_array[0][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;
    
        file << "==== Distrib" << std::endl;
        file << "Type\tParameters" << std::endl;

        for (size_t n = 0; n < cell_porosity_array.size(); n++)
        {
            file << law_type << "\t" << n_topo << "\t";
            for (size_t k = 0; k < n_topo-1; k++)
            {
                file << cell_porosity_array[n][k] << "\t" << topo_list[k] << "\t";
            }
            file << cell_porosity_array[n][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;
        }
    }
    else
    {
        logError() << "Porosity mapper not implemeted for law type" << law_type;
        logFatal(__FILE__, __func__, __LINE__);
    }

    file.close();
    logInfo() << "    ...done";
    //== Cells

    //-- Boundary edges
    std::ofstream file_bnd;
    std::string filename_bnd = m_io->absolutePath(swFile::input_ddp_porosity_edge_bnd_file);
    file_bnd.open(filename_bnd);
    if (file_bnd.is_open() == false)
    {
        logError() << "    Unable to open" << (swFile::input_ddp_porosity_edge_bnd_file).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing edge porosity map :" << swFile::input_ddp_porosity_edge_bnd_file;

    file_bnd << "Unif\t0" << std::endl;
    file_bnd << "NtabMax\t" << n_tab_max << std::endl;
    file_bnd << "==== Default Param" << std::endl;
    file_bnd << "Type\tParameters" << std::endl;
    if ((law_type == -2) || (law_type == -1))
    {
        file_bnd << law_type << "\t" << n_topo << "\t";
        for (size_t k = 0; k < n_topo-1; k++)
        {
            file_bnd << edge_bnd_porosity_array[0][k] << "\t" << topo_list[k] << "\t";
        }
        file_bnd << edge_bnd_porosity_array[0][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;

        file_bnd << "==== Distrib" << std::endl;
        file_bnd << "Type\tParameters" << std::endl;

        for (size_t n = 0; n < edge_bnd_porosity_array.size(); n++)
        {
            file_bnd << law_type << "\t" << n_topo << "\t";
            for (size_t k = 0; k < n_topo-1; k++)
            {
                file_bnd << edge_bnd_porosity_array[n][k] << "\t" << topo_list[k] << "\t";
            }
            file_bnd << edge_bnd_porosity_array[n][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;
        }
    }
    else
    {
        logError() << "Porosity mapper not implemeted for law type" << law_type;
        logFatal(__FILE__, __func__, __LINE__);
    }

    file_bnd.close();
    logInfo() << "    ...done";
    //== Boundary edges

    //-- Internal egdes
    std::ofstream file_int;
    std::string filename_int = m_io->absolutePath(swFile::input_ddp_porosity_edge_int_file);
    file_int.open(filename_int);
    if (file_int.is_open() == false)
    {
        logError() << "    Unable to open" << (swFile::input_ddp_porosity_edge_int_file).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing edge porosity map :" << swFile::input_ddp_porosity_edge_int_file;

    file_int << "Unif\t0" << std::endl;
    file_int << "NtabMax\t" << n_tab_max << std::endl;
    file_int << "==== Default Param" << std::endl;
    file_int << "Type\tParameters" << std::endl;
    if ((law_type == -2) || (law_type == -1))
    {
        file_int << law_type << "\t" << n_topo << "\t";
        for (size_t k = 0; k < n_topo-1; k++)
        {
            file_int << edge_int_porosity_array[0][k] << "\t" << topo_list[k] << "\t";
        }
        file_int << edge_int_porosity_array[0][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;

        file_int << "==== Distrib" << std::endl;
        file_int << "Type\tParameters" << std::endl;

        for (size_t n = 0; n < edge_int_porosity_array.size(); n++)
        {
            file_int << law_type << "\t" << n_topo << "\t";
            for (size_t k = 0; k < n_topo-1; k++)
            {
                file_int << edge_int_porosity_array[n][k] << "\t" << topo_list[k] << "\t";
            }
            file_int << edge_int_porosity_array[n][n_topo-1] << "\t" << topo_list[n_topo-1] << std::endl;
        }
    }
    else
    {
        logError() << "Porosity mapper not implemeted for law type" << law_type;
        logFatal(__FILE__, __func__, __LINE__);
    }

    file_int.close();
    logInfo() << "    ...done";
    //== Internal edges

    return EXIT_SUCCESS;
}

bool swPorosityMapperPipeline::writeDip()
/**
 * @brief write the DIP porosity maps
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- Cells
    std::ofstream file;
    std::string filename = m_io->absolutePath(swFile::input_dip_porosity_cell_file);
    file.open(filename);
    if (file.is_open() == false)
    {
        logError() << "    Unable to open" << (filename).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing cell porosity map :" << filename;

    file << "Unif\t0" << std::endl;
    file << "==== Default Param" << std::endl;
    file << "PhiW" << std::endl;
    file << cell_porosity_list[0] << std::endl;
    file << "==== Distrib" << std::endl;

    for (size_t n = 0; n < cell_porosity_list.size(); n++)
    {
        file << cell_porosity_list[n] << std::endl;
    }

    file.close();
    logInfo() << "    ...done";
    //== Cells

    //-- Boundary edges
    std::ofstream file_bnd;
    std::string filename_bnd = m_io->absolutePath(swFile::input_dip_porosity_edge_bnd_file);
    file_bnd.open(filename_bnd);
    if (file_bnd.is_open() == false)
    {
        logError() << "    Unable to open" << (filename_bnd).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing edge porosity map :" << filename_bnd;

    file_bnd << "Unif\t0" << std::endl;
    file_bnd << "==== Default Param" << std::endl;
    file_bnd << "PhiG" << std::endl;
    file_bnd << edge_bnd_porosity_list[0] << std::endl;;
    file_bnd << "==== Distrib" << std::endl;

    for (size_t n = 0; n < edge_bnd_porosity_list.size(); n++)
    {
        file_bnd << edge_bnd_porosity_list[n] << std::endl;
    }

    file_bnd.close();
    logInfo() << "    ...done";
    //== Boundary edges

    //-- Internal edges
    std::ofstream file_int;
    std::string filename_int = m_io->absolutePath(swFile::input_dip_porosity_edge_int_file);
    file_int.open(filename_int);
    if (file_int.is_open() == false)
    {
        logError() << "    Unable to open" << (filename_int).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing edge porosity map :" << filename_int;

    file_int << "Unif\t0" << std::endl;
    file_int << "==== Default Param" << std::endl;
    file_int << "PhiG" << std::endl;
    file_int << edge_int_porosity_list[0] << std::endl;;
    file_int << "==== Distrib" << std::endl;

    for (size_t n = 0; n < edge_int_porosity_list.size(); n++)
    {
        file_int << edge_int_porosity_list[n] << std::endl;
    }

    file_int.close();
    logInfo() << "    ...done";
    //== Internal edges

    return EXIT_SUCCESS;
}

bool swPorosityMapperPipeline::writeSp()
/**
 * @brief write the SP or porosity map
 * 
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
{
    //-- Cells
    std::ofstream file;
    std::string filename = m_io->absolutePath(swFile::input_sp_porosity_cell_file);
    file.open(filename);
    if (file.is_open() == false)
    {
        logError() << "    Unable to open" << (filename).c_str();
        logFatal(__FILE__, __func__, __LINE__);
        return EXIT_FAILURE;
    }
    logInfo() << "    writing cell porosity map :" << filename;

    file << "Unif\t0" << std::endl;
    file << "==== Default Param" << std::endl;
    file << "PhiW" << std::endl;
    file << cell_porosity_list[0] << std::endl;
    file << "==== Distrib" << std::endl;

    for (size_t n = 0; n < cell_porosity_list.size(); n++)
    {
        file << cell_porosity_list[n] << std::endl;
    }

    file.close();
    logInfo() << "    ...done";
    //== Cells

    return EXIT_SUCCESS;
}

