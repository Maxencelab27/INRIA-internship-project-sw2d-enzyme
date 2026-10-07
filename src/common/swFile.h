#pragma once

#include "../common/swString.h"

namespace swFile
/**
 * @brief set the simulation files structure
 * example : swFile::primary_folder return the path of the primary_files folder (relatively to the project folder)
 */
{
    /// Folders
    const std::string primary_folder = "primary_files" + swString::dirSeparator();
    const std::string input_folder = "input" + swString::dirSeparator();
    const std::string log_folder = "log" + swString::dirSeparator();
    const std::string output_folder = "output" + swString::dirSeparator();

    /// Input files
    const std::string input_boussinesq_file = input_folder + "Boussinesq_map.txt";
    const std::string input_converter_file = primary_folder + "input.conv";
    const std::string input_ddp_porosity_cell_file = input_folder + "DDP_cell_porosity_map.txt";
    const std::string input_ddp_porosity_edge_bnd_file = input_folder + "DDP_edge_bnd_porosity_map.txt";
    const std::string input_ddp_porosity_edge_int_file = input_folder + "DDP_edge_int_porosity_map.txt";
    const std::string input_dip_porosity_cell_file = input_folder + "DIP_cell_porosity_map.txt";
    const std::string input_dip_porosity_edge_bnd_file = input_folder + "DIP_edge_bnd_porosity_map.txt";
    const std::string input_dip_porosity_edge_int_file = input_folder + "DIP_edge_int_porosity_map.txt";
    const std::string input_extraction_time_file = input_folder + "extraction_times_list.txt";
    const std::string input_friction_chezy_map = input_folder + "friction_chezy_map.txt";
    const std::string input_friction_manning_map = input_folder + "friction_manning_map.txt";
    const std::string input_friction_strickler_map = input_folder + "friction_strickler_map.txt";
    const std::string input_hydro_boundary_timeseries = input_folder + "hydro_boundary_time_series.txt";
    const std::string input_hydro_boundary_rating_curve = input_folder + "hydro_boundary_rating_curve.txt";
    const std::string input_infiltration_h_proportional_parameter_map = input_folder + "infiltration_h_proportional_parameter_map.txt";
    const std::string input_infiltration_horton_parameter_map = input_folder + "infiltration_horton_parameter_map.txt";
    const std::string input_hotstart_file = input_folder + "ResultHydro_";
    const std::string input_geo_file = input_folder + "geometry.sw2d";
    const std::string input_initial_conditions_file = input_folder + "initial_conditions_map.txt";
    const std::string input_momentum_dissipation_parameters_map = input_folder + "cell_momentum_dissipation_map.txt";
    const std::string input_mp_porosity_cell_file = input_folder + "MP_cell_porosity_map.txt";
    const std::string input_mp_exchange_cell_file = input_folder + "MP_cell_exchange_coeffs_map.txt";
    const std::string input_postprocess_file = input_folder + "input.postproc";
    const std::string input_porosity_mapper_file = primary_folder + "input.pmap";
    const std::string input_porosity_mapper_buildings_file = primary_folder + "buildings.txt";
    const std::string input_precipitation_stations = input_folder + "precipitation_station_codes_map.txt";
    const std::string input_precipitation_time_series = input_folder + "precipitation_time_series.txt";
    const std::string input_probe_locations_file = input_folder + "probe_locations.txt";
    const std::string input_qsection_file = input_folder + "qsection.txt";
    const std::string input_runoff_infiltration_coefficient_map = input_folder + "infiltration_coefficient_map.txt";
    const std::string input_runoff_runoff_coefficient_map = input_folder + "runoff_coefficient_map.txt";
    const std::string input_runoff_scs_parameter_map = input_folder + "scs_parameter_map.txt";
    const std::string input_singular_head_loss_map = input_folder + "singular_head_loss_map.txt";
    const std::string input_solver_file = input_folder + "input.solv";
    const std::string input_sp_porosity_cell_file = input_folder + "SP_cell_porosity_map.txt";
    const std::string input_wind_stations = input_folder + "wind_station_codes_map.txt";
    const std::string input_wind_time_series = input_folder + "wind_time_series.txt";
    const std::string input_zone_file = primary_folder + "input.zone";

    /// Log files
    const std::string log_boussinesq_file = log_folder + "boussinesq.txt";
    const std::string log_cell_file = log_folder + "cells.txt";
    const std::string log_cell_small_file = log_folder + "cells_small.txt";
    const std::string log_converter_file = log_folder + "console_Converter.log";
    const std::string log_ddp_porosity_cell_file = log_folder + "DDP_porosity_cells.txt";
    const std::string log_ddp_porosity_edge_bnd_file = log_folder + "DDP_porosity_edge_bnd.txt";
    const std::string log_ddp_porosity_edge_int_file = log_folder + "DDP_porosity_edge_int.txt";
    const std::string log_dip_porosity_cell_file = log_folder + "DIP_porosity_cells.txt";
    const std::string log_dip_porosity_edge_bnd_file = log_folder + "DIP_porosity_edge_bnd.txt";
    const std::string log_dip_porosity_edge_int_file = log_folder + "DIP_porosity_edge_int.txt";
    const std::string log_edge_int_file = log_folder + "edges_intern.txt";
    const std::string log_edge_bnd_file = log_folder + "edges_boundary.txt";
    const std::string log_edge_bnd_factors_file = log_folder + "edge_bnd_factors.txt";
    const std::string log_edge_int_factors_file = log_folder + "edge_int_factors.txt";
    const std::string log_friction_file = log_folder + "friction.txt";
    const std::string log_infiltration_file = log_folder + "infiltration.txt";
    const std::string log_input_conv_file = log_folder + "input.conv";
    const std::string log_input_postproc_file = log_folder + "input.postproc";
    const std::string log_input_porosity_mapper_file = log_folder + "input.pmap";
    const std::string log_mp_porosity_cell_file = log_folder + "MP_porosity_cell_parameters.txt";
    const std::string log_node_file = log_folder + "nodes.txt";
    const std::string log_postprocess_file = log_folder + "console_PostProcess.log";
    const std::string log_porositymapper_file = log_folder + "console_PorosityMapper.log";
    const std::string log_runoff_file = log_folder + "runoff.txt";
    const std::string log_singular_head_loss_file = log_folder + "singular_head_loss.txt";
    const std::string log_solver_file = log_folder + "console_Solver.log";
    const std::string log_sp_porosity_cell_file = log_folder + "SP_porosity_cells.txt";
    const std::string log_topo_after_correction_file = log_folder + "topography_after_correction.txt";
    const std::string log_topo_before_correction_file = log_folder + "topography_before_correction.txt";
    const std::string log_topo_corrected_cells_file = log_folder + "topography_corrected_cells.txt";
    const std::string log_wind_file = log_folder + "wind.txt";
    const std::string log_zone_file = log_folder + "console_Zone.log";

    /// Output files
    const std::string output_boundary_fluxes_file = output_folder + "ResultBoundaryFluxes.txt";
    const std::string output_extraction_time_file_root = output_folder + "ResultTime_p";
    const std::string output_hydro_map_file_root = output_folder + "ResultHydro_p";
    const std::string output_transport_map_file_root = output_folder + "ResultTransport_p";
    const std::string output_map_binary_file = output_folder + "ResultHydro.bin";
    const std::string output_map_max_date_file = output_folder + "ResultHydro_DateOfMax.txt";
    const std::string output_map_max_file = output_folder + "ResultHydro_Max.txt";
    const std::string output_map_submersion_durations = output_folder + "SubmersionDurations.txt";
    const std::string output_probes_file = output_folder + "ProbesTimeSeries.txt";
    const std::string output_qsection_file = output_folder + "qsection.txt";
};
