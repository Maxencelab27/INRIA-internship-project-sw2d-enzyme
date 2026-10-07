#pragma once
#include <fstream>
#include <string>
#include <vector>

#include "../../common/swIO.h"

class swBoundaryCondition;
class swBoundaryConditionRating;
class swBoundaryConditionTimeSeries;
class swEdgeBoundary;
class swMesh;
class swSolverParam;

class swParent
{
	public:
        swParent(void);
        ~swParent(void);
        static swParent* create(swSolverParam *solver_param, swIO *io);

        void checkStore(const double& t, const double& dt);
        void computeBC(double const& t);   
        void displayBiblio(void);
        virtual bool execute() = 0;
        virtual void importAttributesFromSolverParam(swSolverParam *solver_param) = 0;
        virtual bool init(swMesh *mesh) = 0;
        bool initBoundaryConditions();
        double interpolateTimeSeries(std::vector<double> const& times_vector, std::vector<double> const& vals_vector, 
                                     double const& t); 

        void setMeshTables(swMesh *mesh); /**< Implemented here, will be inherited by all child (model) classes */
        void setIOPointer(swIO *io);
        void setPrecipitationStationCodes(void);
        void setPrecipitationTimeSeriesTable(void);
        void setWindStationCodes(void);
        void setWindVelocityTimeSeriesTables(void);

        std::vector<std::string> biblio_reference_list;
        swIO *m_io = nullptr;                       /**< pointer to the swIO object */

        //-- BC attributes
        std::vector<swBoundaryConditionRating*> bc_ratingcurve_list = {};      /**< list of the different boundary conditions */
        std::vector<swBoundaryConditionTimeSeries*> bc_timeseries_list = {};   /**< list of the different boundary conditions */
        std::vector<double> bc_times_vector = {};                              /**< [n]: the time that correponds to the n-th record in the BC file */
        std::vector<std::string> bc_type = {};                                 /**< [i]: type for BC with code (index) i */
        std::vector<double> bc_val_interpolated = {};                          /**< [i]: numerical value for BC with code (index) i */
        std::vector<std::vector<double>> bc_vals_table = {};                   /**< [i,n]: BC value for the i-th BC code at the n-th time record */
        //== BC attributes

        ///-- Files attributes
        std::ofstream bnd_fluxes_output_stream;
        std::string probes_filename; /** @todo to be removed ? **/
        std::ofstream probes_output_stream;
        std::ofstream qsection_output_stream;
        //==  Files attributes

        //---- Parabolic step attributes
        std::vector<double> parab_edge_int_flux_factors = {};    /**< [i] The quantity D * L / d_norm for interface i */
        //==== Parabolic step attributes


        //----Precipitation attributes (@todo move to swSource.h)
        std::vector<double> stationwise_precipitation_rate_interpolated = {};         /**< [i]: rainfall rate interpolated at a given time fir station i */
        std::vector<size_t> precipitation_station_code = {};                          /**< [k]: rainfall station code for cell k */
        std::vector<double> precipitation_times_vector = {};                          /**< [n]: the time that correponds to the n-th record in the precipitation time series file */
        std::vector<std::vector<double>> stationwise_precipitation_rate_table = {};   /**< [i,n]: rainfall rate for staion i at time n */
        //====Precipitation attributes


        //---- Wind attributes (@todo move to swSource.h)
        std::vector<double> stationwise_wind_velocity_x_interpolated = {};   /**< [i,n]: x-velocity interpolated at a given time */
        std::vector<double> stationwise_wind_velocity_y_interpolated = {};   /**< [i,n]: x-velocity interpolated at a given time */
        std::vector<size_t> wind_station_code = {};                          /**< [i]:   wind station code for cell i*/
        std::vector<double> wind_times_vector = {};                          /**< [n]: the time that correponds to the n-th record in the wind time series file */
        std::vector<std::vector<double>> wind_velocity_x_table = {};         /**< [i,n]: x-velocity at the i-th wind station at time n */
        std::vector<std::vector<double>> wind_velocity_y_table = {};         /**< [i,n]: x-velocity at the i-th wind station at time n */
        //==== Wind attributes
        
        //---- Result storage
        double dt_map;
        double dt_time_series;
        bool store_map = false;
        bool store_time_series = false;
        //==== Result storage

		//---- Mesh tables for faster CPU 
        size_t bc_type_nb = 0;							/**< number of BC types **/
        size_t edge_total_nb = 0;                       /**< total number of edges **/
        
        //------ boundaries
        std::vector<int> list_bc_codes;                             /**< list of the boundary codes **/
        std::vector<std::vector<swEdgeBoundary*>> list_boundaries; /**< list of the boundary associated to each bc code **/
        //====== boundaries

        //------ q sections
        std::vector<std::string> list_qsections_names;             /**< list of the q section names **/
        std::vector<std::vector<std::string>> list_qsections;      /**< list of the q sections associated to each name **/
        std::vector<std::vector<int>> list_qsections_signs;        /**< list of the qsections signs **/
        //====== q sections
        
		//------ Boundary edges
        size_t edge_bnd_nb = 0;                          /**< number of boundary edges **/
        std::vector<int> edge_bnd_bc_code = {};          /**< [i]: code of the boundary condition for boundary edge #i **/
        std::vector<swBoundaryCondition*> edge_bnd_bc_pointer = {};       /**< [i]: Pointer to the boundary condition object for boundary edge #i **/
        std::vector<size_t> edge_bnd_cell_idx = {};      /**< [i]: Index of the cell on the RHS of boundary edge #i **/
        std::vector<size_t> edge_bnd_idx = {};           /**< [i]: Idx of boundary edge #i **/
        std::vector<double> edge_bnd_length = {};        /**< [i]: length of boundary edge #i **/
        std::vector<double> edge_bnd_midpoint_x = {};    /**< [i]: x-coordinate of midpoint for boundary edge #i **/
        std::vector<double> edge_bnd_midpoint_y = {};    /**< [i]: y-coordinate of midpoint for boundary edge #i **/
        std::vector<double> edge_bnd_midpoint_z = {};    /**< [i]: z-coordinate of midpoint for boundary edge #i **/
        std::vector<double> edge_bnd_node_begin_x = {};  /**< [i]: x-coordinate of begin node for boundary edge #i **/
        std::vector<double> edge_bnd_node_begin_y = {};  /**< [i]: y-coordinate of begin node for boundary edge #i **/
        std::vector<double> edge_bnd_node_end_x = {};    /**< [i]: x-coordinate of end node for boundary edge #i **/
        std::vector<double> edge_bnd_node_end_y = {};    /**< [i]: y-coordinate of end node for boundary edge #i **/
        std::vector<double> edge_bnd_norm_x = {};        /**< [i]: x-component of normal unit vector for boundary edge #i **/
        std::vector<double> edge_bnd_norm_y = {};        /**< [i]: y-component of normal unit vector for boundary edge #i **/
		//====== Boundary edges
		
		//------ Internal edges
        size_t edge_int_nb = 0;                          /**< number of internal edges **/
        std::vector<size_t> edge_int_celL_idx = {};      /**< [i]: Index of the cell on the LHS of internal edge #i **/
        std::vector<size_t> edge_int_celR_idx = {};      /**< [i]: Index of the cell on the RHS of internal edge #i **/
        std::vector<size_t> edge_int_idx = {};           /**< [i]: Idx of internal edge #i **/
        std::vector<double> edge_int_length = {};        /**< [i]: length of internal edge #i **/
        std::vector<double> edge_int_midpoint_x = {};    /**< [i]: x-coordinate of midpoint for internal edge #i **/
        std::vector<double> edge_int_midpoint_y = {};    /**< [i]: y-coordinate of midpoint for internal edge #i **/
        std::vector<double> edge_int_midpoint_z = {};    /**< [i]: z-coordinate of midpoint for internal edge #i **/
        std::vector<double> edge_int_node_begin_x = {};  /**< [i]: x-coordinate of begin node for internal edge #i **/
        std::vector<double> edge_int_node_begin_y = {};  /**< [i]: y-coordinate of begin node for internal edge #i **/
        std::vector<double> edge_int_node_end_x = {};    /**< [i]: x-coordinate of end node for internal edge #i **/
        std::vector<double> edge_int_node_end_y = {};    /**< [i]: y-coordinate of end node for internal edge #i **/
        std::vector<double> edge_int_norm_x = {};        /**< [i]: x-component of normal unit vector for internal edge #i **/
        std::vector<double> edge_int_norm_y = {};        /**< [i]: y-component of normal unit vector for internal edge #i **/
		//====== Internal edges

		
		//------ Cells
        size_t cell_nb = 0;                                       /**< number of cells **/
        std::vector<double> cell_area = {};                       /**< [i]: area of cell #i **/
        std::vector<double> cell_centre_x = {};                   /**< [i]: x-coordinate of cell #i **/
        std::vector<double> cell_centre_y = {};                   /**< [i]: y-coordinate of cell #i **/
        std::vector<std::vector<size_t>> cell_edges_bnd_idx = {}; /**< [i][j]: Index of the j-th boundary edge of the i-th cell */
        std::vector<std::vector<size_t>> cell_edges_int_idx = {}; /**< [i][j]: Index of the j-th internal edge of the i-th cell */
        std::vector<size_t> cell_id = {};                         /**< [i]: Id of cell #i **/
        std::vector<std::vector<size_t>> cell_neigh_cell_idx = {}; /**< [i][j]: Idx of the j-th neighbour cell of the i-th cell */
        std::vector<double> cell_perimeter = {};                  /**< [i]: perimeter of cell #i **/
        std::vector<double> cell_z_b = {};                        /**< [i]: bottom level of cell #i **/
		//====== Cells
	 
        //------ Nodes
        size_t node_nb = 0;                             /**< total number of nodes **/
        std::vector<size_t> node_idx = {};               /**< [i]: Idx of node #i **/
        std::vector<double> node_coord_x = {};     /**< [i]: x-coordinate of node #i **/
        std::vector<double> node_coord_y = {};     /**< [i]: y-coordinate of node #i **/
        std::vector<double> node_coord_z = {};     /**< [i]: z-coordinate of node #i **/
        //====== Nodes	 
		
		//==== Mesh tables for faster CPU

};

