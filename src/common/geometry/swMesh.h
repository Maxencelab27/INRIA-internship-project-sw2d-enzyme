#pragma once

#include <string>
#include <vector>

#include "swCell.h"
#include "swEdgeBoundary.h"
#include "swEdgeIntern.h"
#include "../../common/swIO.h"
#include "swNode.h"

/**
 * @class swMesh
 * @brief Describes the mesh properties (cells, edges, nodes)
 */
class swMesh
{
public:
    //-- General methods
    swMesh(swIO *io);
    ~swMesh() = default;

    bool build();
    bool checkQuality(double a_min);
    swCell* getCellFromCoordinates(double, double);
    bool importGeometry(std::string filename_2dm, std::string filename_bc, size_t default_bc);
    bool readQSections();
    bool topographyCorrection(const std::string& topo_filename,const std::string &param_corr);
    bool writeGeometryFile(const std::string& filename_geo,
                           const std::string& filename_input_conv,
                           const std::string& filename_2dm,
                           const std::string& filename_bc);
    bool writeLog();
    //== General methods

    //-- Getters
    std::vector<int> getBcCodes();
    std::vector<std::vector<swEdgeBoundary *>> getBoundaries();
    std::vector<swCell*> getCellsList();
    std::vector<swEdgeBoundary*> getEdgeBoundaryList();
    std::vector<swEdgeIntern*> getEdgeInternList();
    swIO* getIO();
    std::vector<swNode*> getNodesList();
    std::vector<double> getCellsAreaList();
    std::vector<std::vector<std::string>> getQSections(void);
    std::vector<std::string> getQSectionNames(void);
    std::vector<std::vector<int>> getQSectionSigns(void);
    //== Getters

private:
    //-- General methods
    bool assignBC(std::string filename_bc, size_t default_bc);
    bool buildCells(std::vector<size_t> cells_sms_id, std::vector<std::vector<size_t>> cells_nodes_idx);
    bool buildEdges(size_t default_bc);
    bool read2dm(std::string filename_2dm,
                 std::vector<std::vector<size_t> > &cell_nodes_idx,
                 std::vector<size_t> &cell_sms_id);
    bool readGeometryFile();
    bool readGeometryFileV2(std::ifstream &file);

    // log methods
    bool writeLogCell(std::vector<swCell *> cells, std::string filename);
    bool writeLogEdgeBoundary();
    bool writeLogEdgeIntern();
    bool writeLogNode();
    //== General methods

    //-- Attributes
    swIO *m_io = nullptr;

    //geometry
    std::vector<swCell*> list_cells;                           /**< list of the pointers to cell **/
    std::vector<swEdgeBoundary*> list_edges_boundary;          /**< list of the boundary edges **/
    std::vector<swEdgeIntern*> list_edges_intern;              /**< list of the internal edges **/
    std::vector<swNode*> list_nodes;                           /**< list of the nodes **/
    std::string m_version = "2";                               /**< version number of geometry file **/

    // boundary conditions
    std::vector<int> list_bc_codes;                             /**< list of the bc used **/
    std::vector<std::vector<swEdgeBoundary*>> list_boundaries; /**< list of the boundary associated to each bc code **/

    // q sections
    std::vector<std::vector<std::string>> list_qsections;      /**< list of the q sections associated to each name **/
    std::vector<std::string> list_qsections_names;             /**< list of the q section names **/
    std::vector<std::vector<int>> list_qsections_signs;        /**< list of the qsections signs **/
    //== Attributes
};
