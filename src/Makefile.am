bin_PROGRAMS = PGFem3D BPA_driver unit_test_43_driver IVD_driver J2D_driver CH_layer_driver comm_hints_driver
all_sources = \
ALM.cc \
Arc_length.cc \
boomerAMGInterface.cc \
GRedist_node.cc \
LINE.cc \
Newton_Raphson.cc \
PLC.cc \
PLoc_Sparse.cc \
Psparse_ApAi.cc \
Re1_Re2_Re3.cc \
TA_GA.cc \
addresses.cc \
build_distribution.cc \
def_grad.cc \
differentiation_of_UU.cc \
elem3d.cc \
fd_increment.cc \
fd_residuals.cc \
homogen.cc \
hypre_interface.cc \
in.cc \
initial_guess.cc \
integration.cc \
load.cc \
localizat.cc \
macroscopic_load_AL.cc \
matice.cc \
matrix_printing.cc \
metis_printing.cc \
null.cc \
out.cc \
press_theta.cc \
pressu_shape.cc \
Printing.cc \
print_dist.cc \
read_cryst_plast.cc \
renumber_ID.cc \
res_fini_def.cc \
resice.cc \
resid_on_elem.cc \
rn_skyline.cc \
rowlength.cc \
set_fini_def.cc \
SetGlobalNodeNumbers.cc \
skyline.cc \
stabilized.cc \
stiffmat_fd.cc \
stiffmatel_fd.cc \
stress_strain.cc \
subdivision.cc \
tensors.cc \
utils.cc \
incl.cc \
profiler.cc \
potential.cc \
computeMacroF.cc \
computeMacroS.cc \
vtk_output.cc \
MINI_element.cc \
MINI_3f_element.cc \
two_field_element.cc \
Hu_Washizu_element.cc \
new_potentials.cc \
displacement_based_element.cc \
PGFem3D_options.cc \
gen_path.cc \
volumetric_damage.cc \
initialize_damage.cc \
vol_damage_int_alg.cc \
interface_macro.cc \
compute_reactions.cc \
bounding_element.cc \
bounding_element_utils.cc \
get_dof_ids_on_elem.cc \
get_ndof_on_elem.cc \
transform_coordinates.cc \
generate_dof_ids.cc \
node.cc \
element.cc \
solve_system.cc \
PGFEM_HYPRE_preconditioners.cc \
allocation.cc \
cohesive_element_utils.cc \
cohesive_element.cc \
cohesive_potentials.cc \
quadrature_rules.cc \
integrate_surface.cc \
applied_traction.cc \
pgfem_comm.cc \
tfa.cc \
hommat.cc \
sig.cc \
supp.cc \
mesh_load.cc \
material.cc \
eps.cc \
matgeom.cc \
crpl.cc \
ensight.cc \
read_input_file.cc \
microscale_information.cc \
ms_cohe_job_info.cc \
ms_cohe_job_list.cc \
PGFEM_par_matvec.cc \
compute_ms_cohe_job.cc \
PGFEM_io.cc \
PGFEM_mpi.cc \
ms_job_intercomm.cc \
macro_micro_functions.cc \
pgf_fe2_micro_server.cc \
pgf_fe2_job.cc \
pgf_fe2_server_rebalance.cc \
pgf_fe2_rebalancer.cc \
pgf_fe2_macro_client.cc \
pgf_fe2_restart.cc \
pgf_fe2_compute_max_n_jobs.cc \
solver_file.cc \
condense.cc \
three_field_element.cc \
femlib.cc \
macro_micro_main.cc \
main.cc \
dynamics.cc \
plasticity_model.cc \
plasticity_model_none.cc \
plasticity_model_BPA.cc \
cm_iso_viscous_damage.cc \
cm_j2_plasticity.cc \
cm_placeholder_functions.cc \
cm_uqcm.cc \
state_variables.cc \
constitutive_model.cc \
post_processing.cc \
comm_hints.cc \
restart.cc \
PGFem3D_data_structure.cc

comm_hints_driver_SOURCES = $(all_sources) comm_hints.cc comm_hints_driver.cc

BPA_driver_SOURCES = $(all_sources) plasticity_model_BPA_driver.cc

IVD_driver_SOURCES = $(all_sources) cm_iso_viscous_damage_driver.cc

J2D_driver_SOURCES = $(all_sources) cm_j2_plasticity_driver.cc

CH_layer_driver_SOURCES = $(all_sources) CH_layer_driver.cc

unit_test_43_driver_SOURCES = $(all_sources) unit_test_43_driver.cc

PGFem3D_SOURCES = $(all_sources) super_main.cc

VTK_IO_include = -I$(top_srcdir)/lib/VTK_IO/include
VTK_IO_lib = $(top_builddir)/lib/VTK_IO/src/libPGFem3D_to_VTK.a

energy_equation_include = -I$(top_srcdir)/energy_equation/include
energy_equation_lib = $(top_builddir)/energy_equation/src/libenergy_equation.a

# Universal flags for all binaries
AM_CXXFLAGS = -fpermissive
AM_CPPFLAGS = -Wall -Werror -Wno-misleading-indentation -I$(top_srcdir)/include -I$(top_builddir)/include $(VTK_IO_include) $(energy_equation_include) $(hypre_include) $(suitesparse_include) $(cnstvm_include) $(mkl_include)
LDADD       = $(VTK_IO_lib) $(energy_equation_lib) $(suitesparse_lib) $(hypre_lib) $(vtk_lib) $(cnstvm_lib) $(mkl_lib)
AM_LDFLAGS  = -Wl,--allow-multiple-definition
