#include "../common/geometry/swCell.h"
#include "../common/geometry/swMesh.h"
#include "swPorosityMapperParam.h"
#include "swRasterBuildings.h"
#include "../common/swString.h"

class swPorosityMapperPipeline
{
public:
    //-- General methods
    swPorosityMapperPipeline(swIO *io, swPorosityMapperParam *param, swMesh *mesh);
    ~swPorosityMapperPipeline() = default;

    void init();
    bool mapPorosity(std::string model);
    //== General methods
    
private:
    //-- General methods
    bool computeDDP();
    bool computeDIP();
    bool computeSP();

    double cellPorosity(swCell *cell, double h);
    double edgePorosity(swEdgeBoundary *edge, double h);
    double edgePorosity(swEdgeIntern *edge, double h);

    bool writeDdp(int law_type, int n_tab_max);
    bool writeDip();
    bool writeSp();
    //== General methods

    //-- Attributes
    std::vector<std::vector<double>> cell_porosity_array;          /**< table of the cell porosities for each heighs (DDP) */
    std::vector<double> cell_porosity_list;                        /**< vector of the cell porosities (SP + DIP) */
    std::vector<std::vector<double>> edge_bnd_porosity_array;      /**< table of the boundary edge porosities for each heighs (DDP) */
    std::vector<std::vector<double>> edge_int_porosity_array;      /**< table of the internal edge porosities for each heighs (DDP) */
    std::vector<double> edge_bnd_porosity_list;                    /**< vector of the boundary cell porosities (DIP) */
    std::vector<double> edge_int_porosity_list;                    /**< vector of the internal cell porosities (DIP) */
    std::vector<double> topo_list;                                 /**< vector of the building heighs (in ascending order) */

    // raster
    swRasterBuildings *m_buildings;                                /**< table of spatial raster containing the buildings heigh */

    // footprints
    std::vector<std::vector<double>> x_buildings;                  /**< list of x-coordinates for each building */
    std::vector<std::vector<double>> y_buildings;                  /**< list of y-coordinates for each building */
    std::vector<double> h_buildings;                               /**< list of the heigh of each building */

    swIO *m_io = nullptr;
    swMesh *m_mesh = nullptr;
    swPorosityMapperParam *m_param = nullptr;
    //== Attributes
};
