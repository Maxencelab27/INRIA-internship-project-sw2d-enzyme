#include "swRasterBuildings.h"

#include <cmath>
#include <filesystem>

#include "../common/geometry/swCell.h"
#include "../common/geometry/swEdgeBoundary.h"
#include "../common/geometry/swEdgeIntern.h"
#include "../common/swIO.h"
#include "../common/swLog.h"
#include "../common/geometry/swNode.h"
#include "../common/swString.h"

/**
 * @brief Construct a new sw Raster Buildings::sw Raster Buildings object
 * 
 */
swRasterBuildings::swRasterBuildings(void)
{
}

/**
 * @brief Destroy the sw Raster Buildings::sw Raster Buildings object
 * 
 */
swRasterBuildings::~swRasterBuildings(void)
{
}

/**
 * @brief Read the buildings file to assign values to attributes
 * @param[in] filename path + name of the buildings file
 * @return bool EXIT_SUCCESS/EXIT_FAILURE
 */
bool swRasterBuildings::build(const std::string &filename)
{
    
    std::ifstream file;
    file.open(filename);
    if (file.is_open() == false)
    {
        logError() << "    Unable to open file :" << filename;
        logFatal(__FILE__, __func__, __LINE__);
    }
    else
    {
        std::vector<std::vector<double>> t_buildings;

        // read the header
        std::string line_txt = swIO::readLine(file);
        
        //read the first line
        line_txt = swIO::readLine(file);
        std::vector<std::string> line_txt_list = swString::splitStr(line_txt, swString::regexWhiteSpace());

        bool conversion_success = true;
        double x_i = swString::stringToDouble(line_txt_list[0], conversion_success);
        double y_i = swString::stringToDouble(line_txt_list[1], conversion_success);
        double building_i = swString::stringToDouble(line_txt_list[2], conversion_success);

        double x_prev = x_i;
        double y_prev = y_i;
        double x_0 = x_i;
        double y_0 = y_i;

        // add first elements to the vectors
        std::vector<double> buildings_line;
        buildings_line.push_back(building_i);
        topo_list.push_back(building_i);
        x_coords.push_back(x_i);
        y_coords.push_back(y_i);

        // read the following lines
        line_txt = swIO::readLine(file);
        while (!file.eof())
        {
            line_txt_list = swString::splitStr(line_txt, swString::regexWhiteSpace());

            double x_i = swString::stringToDouble(line_txt_list[0], conversion_success);
            double y_i = swString::stringToDouble(line_txt_list[1], conversion_success);
            double building_i = swString::stringToDouble(line_txt_list[2], conversion_success);
            
            if (y_i != y_prev)
            {
                if (y_i > y_prev)
                {
                    logError() << "File" << filename << "is misformated :";
                    logError() << "  y values should be in descending order";
                    logFatal(__FILE__, __func__, __LINE__);
                }
                t_buildings.push_back(buildings_line);
                buildings_line.clear();
                buildings_line.push_back(building_i);
                y_coords.push_back(y_i);
                y_prev = y_i;
                x_prev = x_0;
            }
            else
            {
                if (x_i < x_prev)
                {
                    logError() << "File" << filename << "is misformated :";
                    logError() << "  x values should be in ascending order for each value of y.";
                    logError() << "x_i" << x_i << "x_prev" << x_prev;
                    logFatal(__FILE__, __func__, __LINE__);
                }
                buildings_line.push_back(building_i);
                if (y_i == y_0)
                {
                    x_coords.push_back(x_i);
                }
                x_prev = x_i;
            }

            // add the building heigh to topo_list if it's not already in
            if (std::find(topo_list.begin(), topo_list.end(), building_i) == topo_list.end())
            {
                topo_list.push_back(building_i);
            }

            line_txt = swIO::readLine(file);
        }
        t_buildings.push_back(buildings_line);

        // Transpose buildings so y is in ascending order
        buildings = t_buildings;
        std::reverse(buildings.begin(), buildings.end());

        // Transpose y so it is in ascending order
        std::reverse(y_coords.begin(), y_coords.end());

        // Sort topo_list in ascending order
        std::sort(topo_list.begin(), topo_list.end());

        delta_x = x_coords[1]-x_coords[0];
        delta_y = y_coords[1]-y_coords[0];
        }
    
    return EXIT_SUCCESS;
}

/**
 * @brief find the indices of x and y coordinates of the cell bounding box in the buildings raster file
 * @param[in] cell pointer to the cell
 * @return vector containing minimal x, minimal y, maximal x, maximal y of the bounding box
 */
std::vector<size_t> swRasterBuildings::getBoundingBox(swCell *cell)
{
    // Find min and max x and y for the cell
    double min_x_cell = cell->getNodes()[0]->getX();
    double min_y_cell = cell->getNodes()[0]->getY();
    double max_x_cell = cell->getNodes()[0]->getX();
    double max_y_cell = cell->getNodes()[0]->getY();

    for (swNode *node : cell->getNodes())
    {
        if (node->getX() < min_x_cell) {min_x_cell = node->getX();}
        if (node->getY() < min_y_cell) {min_y_cell = node->getY();}
        if (node->getX() > max_x_cell) {max_x_cell = node->getX();}
        if (node->getY() > max_y_cell) {max_y_cell = node->getY();}
    }

    // Return corresponding indices in the buildings tab
    size_t i_x_min = getIndexX(min_x_cell);
    size_t i_y_min = getIndexY(min_y_cell);
    size_t i_x_max = getIndexX(max_x_cell) + 1;
    size_t i_y_max = getIndexY(max_y_cell) + 1;

    // Special cases
    if (i_x_max > x_coords.size() - 1) {i_x_max -= 1;}
    if (i_y_max > y_coords.size() - 1) {i_y_max -= 1;}

    std::vector<size_t> result = {i_x_min, i_y_min, i_x_max, i_y_max};

    return(result);
}

/**
 * @brief find the indices of x and y coordinates of the edge bounding box in the buildings raster file
 * @param[in] edge pointer to the edge
 * @return vector containing minimal x, minimal y, maximal x, maximal y of the bounding box
 */
std::vector<size_t> swRasterBuildings::getBoundingBox(swEdgeBoundary *edge)
{
    // Find min and max x and y for the edge
    double min_x_edge = std::min(edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX());
    double min_y_edge = std::min(edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY());
    double max_x_edge = std::max(edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX());
    double max_y_edge = std::max(edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY());

    // Return corresponding indices in the buildings tab
    size_t i_x_min = getIndexX(min_x_edge);
    size_t i_y_min = getIndexY(min_y_edge);
    size_t i_x_max = getIndexX(max_x_edge) + 1;
    size_t i_y_max = getIndexY(max_y_edge) + 1;

    // Special cases
    if (max_x_edge == x_coords[0] - delta_x/2) {i_x_max -= 1;}
    if (max_y_edge == y_coords[0] - delta_y/2) {i_y_max -= 1;}
    if (i_x_max > x_coords.size() - 1) {i_x_max -= 1;}
    if (i_y_max > y_coords.size() - 1) {i_y_max -= 1;}

    std::vector<size_t> result = {i_x_min, i_y_min, i_x_max, i_y_max};

    return(result);
}

/**
 * @brief find the indices of x and y coordinates of the edge bounding box in the buildings raster file
 * @param[in] edge pointer to the edge
 * @return vector containing minimal x, minimal y, maximal x, maximal y of the bounding box
 */
std::vector<size_t> swRasterBuildings::getBoundingBox(swEdgeIntern *edge)
{
    // Find min and max x and y for the edge
    double min_x_edge = std::min(edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX());
    double min_y_edge = std::min(edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY());
    double max_x_edge = std::max(edge->getNodes()[0]->getX(), edge->getNodes()[1]->getX());
    double max_y_edge = std::max(edge->getNodes()[0]->getY(), edge->getNodes()[1]->getY());

    // Return corresponding indices in the buildings tab
    size_t i_x_min = getIndexX(min_x_edge);
    size_t i_y_min = getIndexY(min_y_edge);
    size_t i_x_max = getIndexX(max_x_edge) + 1;
    size_t i_y_max = getIndexY(max_y_edge) + 1;

    // Special cases
    if (max_x_edge == x_coords[0] - delta_x/2) {i_x_max -= 1;}
    if (max_y_edge == y_coords[0] - delta_y/2) {i_y_max -= 1;}
    if (i_x_max > x_coords.size() - 1) {i_x_max -= 1;}
    if (i_y_max > y_coords.size() - 1) {i_y_max -= 1;}

    std::vector<size_t> result = {i_x_min, i_y_min, i_x_max, i_y_max};

    return(result);
}

/**
 * @brief Return the array containing the buildings heighs
 * @return array buildings
 */
std::vector<std::vector<double>> swRasterBuildings::getBuildings()
{
    return(buildings);
}

/**
 * @brief Return the numerical difference between 2 consecutive x coordinates
 * @return delta_x value 
 */
double swRasterBuildings::getDeltaX(){
    return(delta_x);
}

/**
 * @brief Return the numerical difference between 2 consecutive y coordinates
 * @return delta_y value 
 */
double swRasterBuildings::getDeltaY(){
    return(delta_y);
}

/**
 * @brief Find the index of greater x[i] coordinate smaller than x parameter
 * @param[in] x coordinate whose index is searched
 * @return index of greater x[i] smaller than x
 */
size_t swRasterBuildings::getIndexX(double x)
{
    size_t i = 0;

    if (x < x_coords[0] - delta_x/2)
    {
        logError() << "Argument" << swString::numberToStr(x, 'e', 10) << "is smaller than minimal x :" << swString::numberToStr(x_coords[0] - delta_x/2, 'e', 10);
        logFatal(__FILE__, __func__, __LINE__);
    }
    if (x > x_coords[x_coords.size()-1] + delta_x/2)
    {
        logError() << "Argument" << swString::numberToStr(x, 'e', 10) << "is greater than maximal x :" << swString::numberToStr(x_coords[x_coords.size()-1] + delta_x/2, 'e', 10);
        logFatal(__FILE__, __func__, __LINE__);
    }

    if (x < x_coords[0] + delta_x/2)
    {
        return(0);
    }

    while ((i < x_coords.size() - 1) && (x > x_coords[i+1] - delta_x/2))
    {
        i++;
    }
    
    return(i);
}

/**
 * @brief Find the index of greater y[i] coordinate smaller than y parameter
 * @param[in] y coordinate whose index is searched
 * @return index of greater y[i] smaller than y
 */
size_t swRasterBuildings::getIndexY(double y)
{
    size_t i = 0;

    if (y < y_coords[0] - delta_y/2)
    {
        logError() << "Argument" << swString::numberToStr(y, 'e', 10) << "is smaller than minimal y :" << swString::numberToStr(y_coords[0] - delta_y/2, 'e', 10);
        logFatal(__FILE__, __func__, __LINE__);
    }
    if (y > y_coords[y_coords.size()-1] + delta_y/2)
    {
        logError() << "Argument" << swString::numberToStr(y, 'e', 10) << "is greater than maximal y :" << swString::numberToStr(y_coords[y_coords.size()-1] + delta_y/2, 'e', 10);
        logFatal(__FILE__, __func__, __LINE__);
    }

    if (y < y_coords[1] + delta_y/2)
    {
        return(0);
    }

    while ((i < y_coords.size() - 1) &&(y > y_coords[i+1] - delta_y/2))
    {
        i++;
    }
    
    return(i);
}

/**
 * @brief Return the vector containing the differrent buildings heigh (in ascending order)
 * @return topo_list vector
 */
std::vector<double> swRasterBuildings::getTopoList(void)
{
    return(topo_list);
}

/**
 * @brief Return a vector of x coordinates of the raster file
 * @return x_coord vector
 */
std::vector<double> swRasterBuildings::getX(){
    return(x_coords);
}

/**
 * @brief Return a vector of y coordinates of the raster file
 * @return y_coord vector
 */
std::vector<double> swRasterBuildings::getY(){
    return(y_coords);
}