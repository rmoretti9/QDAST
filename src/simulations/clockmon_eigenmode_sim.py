import logging
import sys
import ast

from qdast import defaults
from kqcircuits import defaults
from kqcircuits.pya_resolver import pya
from kqcircuits.simulations.export.ansys.ansys_export import export_ansys
from kqcircuits.simulations.export.elmer.elmer_export import export_elmer
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
)
from kqcircuits.simulations.export.simulation_export import export_simulation_oas
from kqcircuits.simulations.post_process import PostProcess
from kqcircuits.util.export_helper import (
    create_or_empty_tmp_directory,
    get_active_or_new_layout,
    open_with_klayout_or_default_application,
    get_simulation_directory,
)
from kqcircuits.simulations.simulation import Simulation
from kqcircuits.simulations.single_element_simulation import (
    get_single_element_sim_class,
)
from qdast.qubits.clockmon import Clockmon
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
    visualise_xsection_cut_on_original_layout,
)
from kqcircuits.simulations.export.xsection.xsection_export import (
    create_xsections_from_simulations,
    separate_signal_layer_shapes,
    visualise_xsection_cut_on_original_layout,
)

sim_tools = ["elmer", "eigenmode"]
# sim_tools = ["eigenmode"]
Lj = 1.4348926834765362e-08

flip_chip = False
etch_opposite_face = False


mer_dims_default = {
    "gap": 7.0,
    "vacuum": 3.0,
    "metal": 7.0,
    "substrate": 3.0,
}
center = (1000, 1000)

# Simulation parameters
sim_class = get_single_element_sim_class(
    Clockmon,
    ignore_ports=["port_drive", "port_1", "port_2", "port_3", "port_4", "port_5", "port_island1_signal", "port_island2_signal"],
)
sim_class.junction_inductance = Lj # Manually adjusting Lj
sim_class.sim_tool = "eig"
sim_parameters = {
    "name": f"clockmon",
    # 'use_internal_ports': True,
    # 'use_ports': True,
    "qubit_face": ["1t1"],
    "face_stack": ["1t1"],
    "with_squid": False,
    "junction_inductance": Lj,
    "box": pya.DBox(pya.DPoint(0, 0), pya.DPoint(2000, 2000)),
    "tls_layer_thickness": [
        4.8e-9 * 1e6,
        0.3e-9 * 1e6,
        2.4e-9 * 1e6,
    ],  # MA, MS, and SA,
    "tls_sheet_approximation": True,
    "metal_height": 0.2,
    "tls_layer_material": ["oxideMA", "oxideMS", "oxideSA"],
    "material_dict": {
        **ast.literal_eval(Simulation.material_dict),
        "oxideMA": {"permittivity": 8},
        "oxideMS": {"permittivity": 11.4},
        "oxideSA": {"permittivity": 4},
    },
    "ground_gap": [630, 610],
    "a": 10,
    "b": 6,
    "coupler_widths": [150, 0, 0, 0, 0, 0],
    "island_extent": [535, 200], 
    "island_to_island_distance": 50,
    "coupler_offsets": [255, 0, 0, 0, 0, 0],
    "clock_diameter": 90,
    "sim_tool": "eig",
    "bending_angle": 0,
}

mer_correction_path = get_simulation_directory(
    f"{sim_parameters['name']}_elmer_metal_edge"
)

for sim_tool in sim_tools:
    is_elmer = sim_tool == "elmer"
    dir_path = create_or_empty_tmp_directory(f"{sim_parameters['name']}_{sim_tool}")

    #######################
    # SIMULATION SETTINGS #
    #######################

    ansys_export_parameters = {
        "path": dir_path,
        "exit_after_run": False,
    }

    ansys_export_parameters.update(
        {
            "ansys_tool": "eigenmode",
            "max_delta_f": 0.05,  # quite tight
            "maximum_passes": 25,
            "minimum_passes": 1,
            "minimum_converged_passes": 2,
            "n_modes": 1,
            "min_frequency": 1,  # minimum allowed frequency
            "mesh_size": {
                "1t1_substratemer": 20,
                "1t1_vacuummer": 20,
            },
            # mer_correction_path is not portable! TODO: make relative path
            "post_process": PostProcess(
                "produce_epr_table.py",
                sheet_approximations={
                    "MA": {"thickness": 4.8e-9, "eps_r": 8, "background_eps_r": 1.0},
                    "MS": {
                        "thickness": 0.3e-9,
                        "eps_r": 11.4,
                        "background_eps_r": 11.45,
                    },
                    "SA": {"thickness": 2.4e-9, "eps_r": 4, "background_eps_r": 11.45},
                },
                groups=["MA", "MS", "SA", "substrate", "vacuum"],
                mer_correction_paths=[
                    ("default", str(mer_correction_path)),
                ],
            ),
            "integrate_energies": True,
        }
    )

    quiet = True
    dielectric_surfaces = {
        "b_substrate": {
            "tan_delta_surf": 5e-7,
        },
        "ma_layer": {
            "tan_delta_surf": 9.9e-3,
            "th": 4.8e-9 * 1e6,
            "eps_r": 8,
        },
        "ms_layer": {
            "tan_delta_surf": 2.6e-3,
            "th": 0.3e-9 * 1e6,
            "eps_r": 11.4,
        },
        "sa_layer": {
            "tan_delta_surf": 2.1e-3,
            "th": 2.4e-9 * 1e6,
            "eps_r": 4,
        },
    }

    elmer_export_parameters = {
        "path": mer_correction_path,
        "tool": "cross-section",
        "linear_system_method": "mg",
        "p_element_order": 3,
        "is_axisymmetric": False,
        "workflow": {
            "run_klayout_gui": not quiet,  # this is visual view of the results which can be
            # removed to speed up the process
            "run_gmsh": True,
            "run_gmsh_gui": not quiet,  # this is visual view of the results which can be
            # removed to speed up the process
            "run_elmergrid": True,
            "run_elmer": True,
            "run_paraview": not quiet,  # this is visual view of the results which can be
            # removed to speed up the process
            "python_executable": "python",  # use 'kqclib' when using singularity image (you can also put a full path)
            "elmer_n_processes": 1,  # -1 means all the physical cores
            "gmsh_n_threads": 1,  # -1 means all the physical cores
        },
        "mesh_size": {
            "vacuum": 1000,
            "b_substrate": 1000,
            "b_signal_1": 1,
            "b_ground": 1,
            "ma_layer": 0.1,
            "ms_layer": 0.1,
            "sa_layer": 0.1,
            "sa_layer&ma_layer": [0.001, 0.010, 0.4],
            #'gap': [4, 8, 200],
            #'port_min_mesh_size': [4., 8, 200],
        },
        "integrate_energies": True,
        "post_process": [
            PostProcess(
                "produce_q_factor_table.py",
                **{k: v["tan_delta_surf"] for k, v in dielectric_surfaces.items()},
            ),
            PostProcess(
                "produce_epr_table.py", groups=["ma", "ms", "sa", "substrate", "vacuum"]
            ),
        ],
    }

    #####################
    # GEOMETRY SETTINGS #
    #####################

    # Get layout
    logging.basicConfig(level=logging.INFO, stream=sys.stdout)
    layout = get_active_or_new_layout()
    # layout.dbu = 1e-5  # need finer DBU for SA-interface

    metal_edge_region_x = 40.0
    sweep_parameters_list = [
        {
            "name": "clockmon",  # required when specifying sims manually
            "partition_regions": [
                {
                    "name": "mer",
                    "face": "1t1",
                    "metal_edge_dimensions": [
                        mer_dims_default["gap"],
                        mer_dims_default["metal"],
                    ],
                    "vertical_dimensions": [
                        mer_dims_default["substrate"],
                        mer_dims_default["vacuum"],
                    ],
                },
            ],
        }
    ]

    # Sweep simulations
    simulations = []  # [sim_class(layout, **sim_parameters)]
    simulations += [
        sim_class(layout, **{**sim_parameters, **sweep_parameters})
        for sweep_parameters in sweep_parameters_list
    ]

    # Export files
    if is_elmer:
        for simulation in simulations:
            separate_signal_layer_shapes(simulation)

        elmer_export_parameters["mesh_size"] = {
            "ma_layer_straight": 0.1,
            "ms_layer_straight": 0.1,
            "sa_layer_straight": 0.1,
            "ma_layer_bend": 0.05,
            "ms_layer_bend": 0.05,
            "sa_layer_bend": 0.05,
            "ma_layer_straight_mer": 0.025,
            "ms_layer_straight_mer": 0.025,
            "sa_layer_straight_mer": 0.025,
            "ma_layer_bend_mer": 0.025,
            "ms_layer_bend_mer": 0.025,
            "sa_layer_bend_mer": 0.025,
            "ma_layer": 0.02,
            "ms_layer": 0.02,
            "sa_layer": 0.02,
            "ma_layer_mer": 0.0025,
            "ms_layer_mer": 0.0025,
            "sa_layer_mer": 0.0025,
            "vacuum_mer": 0.08,
            "b_substrate_mer": 0.08,
        }

        # Create cross sections using xsection tool
        elmer_export_parameters_default = elmer_export_parameters.copy()
        elmer_export_parameters_default["boundary_conditions"] = {
            "xmax": {"potential": 0}
        }
        cords_list = []
        for sweep_parameters in sweep_parameters_list:
            cord1 = pya.DPoint(1248, 1060)
            cord2 = pya.DPoint(1248 + metal_edge_region_x, 1060)
            cords_list.append((cord1, cord2))
        xsection_simulations_metal_edge = create_xsections_from_simulations(
            simulations,
            mer_correction_path,
            cords_list,
            ma_permittivity=dielectric_surfaces["ma_layer"]["eps_r"],
            ms_permittivity=dielectric_surfaces["ms_layer"]["eps_r"],
            sa_permittivity=dielectric_surfaces["sa_layer"]["eps_r"],
            ma_thickness=dielectric_surfaces["ma_layer"]["th"],
            ms_thickness=dielectric_surfaces["ms_layer"]["th"],
            sa_thickness=dielectric_surfaces["sa_layer"]["th"],
            magnification_order=1,
            process_path=defaults.XSECTION_PROCESS_PATH,
            vertical_cull=(
                -mer_dims_default["substrate"] * 6,
                mer_dims_default["vacuum"] * 6,
            ),
            mer_box=pya.DBox(
                metal_edge_region_x / 2.0 - mer_dims_default["metal"],
                -mer_dims_default["substrate"],
                metal_edge_region_x / 2.0 + mer_dims_default["gap"],
                +mer_dims_default["vacuum"],
            ),
        )

        export_elmer(xsection_simulations_metal_edge, **elmer_export_parameters_default)
    else:
        export_ansys(simulations, **ansys_export_parameters, skip_errors=True)


if "elmer" in sim_tools:
    visualise_xsection_cut_on_original_layout(simulations, cords_list)
    open_with_klayout_or_default_application(
        export_simulation_oas(
            simulations, mer_correction_path, file_prefix="xsection_cut_preview"
        )
    )
    oas_file_metal_edge = export_simulation_oas(
        xsection_simulations_metal_edge, dir_path
    )
    open_with_klayout_or_default_application(oas_file_metal_edge)
