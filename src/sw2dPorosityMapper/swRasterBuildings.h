#include <map>
#include <string>
#include <vector>

class swCell;
class swEdgeBoundary;
class swEdgeIntern;
class swNode;
class swParam;

class swRasterBuildings
{
public:
    //-- General methods
    swRasterBuildings(void);
    ~swRasterBuildings(void);

    bool build(const std::string &filename);
    std::vector<size_t> getBoundingBox(swCell *cell);
    std::vector<size_t> getBoundingBox(swEdgeBoundary *edge);
    std::vector<size_t> getBoundingBox(swEdgeIntern *edge);
    //== General methods

    //-- Getters
    std::vector<std::vector<double>> getBuildings(void);
    double getDeltaX(void);
    double getDeltaY(void);
    size_t getIndexX(double x);
    size_t getIndexY(double y);
    std::vector<double> getTopoList(void);
    std::vector<double> getX(void);
    std::vector<double> getY(void);
    //== Getters

private :
    //-- Attributes
    std::vector<std::vector<double>> buildings;    /**< array of topography : buildings[j][i] is the height of buildings at (x[i],y[j]) */
    double delta_x;                                /**< difference between two consecutive x coordinates */
    double delta_y;                                /**< difference between two consecutive y coordinates */
    std::vector<double> topo_list;                 /**< list of the buildings height */
    std::vector<double> x_coords;                  /**< list of x coordinates in asecending order */
    std::vector<double> y_coords;                  /**< list of y coordinates in asecending order */
    //== Attributes
};
