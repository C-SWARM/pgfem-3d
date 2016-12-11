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
differentiation_of_UU.c \
elem3d.cc \
fd_increment.cc \
fd_residuals.cc \
homogen.cc \
hypre_interface.cc \
in.cc \
initial_guess.cc \
integration.c \
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
utils.c \
incl.cc \
profiler.cc \
potential.cc \
computeMacroF.cc \
computeMacroS.cc \
vtk_output.c \
MINI_element.cc \
MINI_3f_element.cc \
two_field_element.cc \
Hu_Washizu_element.cc \
new_potentials.cc \
displacement_based_element.cc \
PGFem3D_options.cc \
gen_path.cc \
volumetric_damage.c \
initialize_damage.cc \
vol_damage_int_alg.cc \
interface_macro.cc \
compute_reactions.cc \
bounding_element.cc \
bounding_element_utils.c \
get_dof_ids_on_elem.cc \
get_ndof_on_elem.cc \
transform_coordinates.cc \
generate_dof_ids.cc \
node.cc \
element.cc \
solve_system.cc \
PGFEM_HYPRE_preconditioners.c \
allocation.cc \
cohesive_element_utils.cc \
cohesive_element.cc \
cohesive_potentials.cc \
quadrature_rules.cc \
integrate_surface.cc \
applied_traction.c \
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
ms_cohe_job_list.c \
PGFEM_par_matvec.cc \
compute_ms_cohe_job.cc \
PGFEM_io.cc \
PGFEM_mpi.cc \
ms_job_intercomm.cc \
macro_micro_functions.c \
pgf_fe2_micro_server.c \
pgf_fe2_job.cc \
pgf_fe2_server_rebalance.cc \
pgf_fe2_rebalancer.c \
pgf_fe2_macro_client.c \
pgf_fe2_restart.cc \
pgf_fe2_compute_max_n_jobs.cc \
solver_file.c \
condense.cc \
three_field_element.c \
femlib.cc \
macro_micro_main.c \
main.c \
dynamics.cc \
plasticity_model.c \
plasticity_model_none.c \
plasticity_model_BPA.c \
cm_iso_viscous_damage.c \
cm_j2_plasticity.c \
cm_placeholder_functions.cc \
cm_uqcm.c \
state_variables.cc \
constitutive_model.c \
post_processing.c \
comm_hints.cc \
restart.c

comm_hints_driver_SOURCES = comm_hints.cc comm_hints_driver.cc

BPA_driver_SOURCES = \
$(all_sources) \
plasticity_model_BPA_driver.c

IVD_driver_SOURCES = \
$(all_sources) \
cm_iso_viscous_damage_driver.c

J2D_driver_SOURCES = \
$(all_sources) \
cm_j2_plasticity_driver.c

CH_layer_driver_SOURCES = \
$(all_sources) \
CH_layer_driver.c

unit_test_43_driver_SOURCES = \
$(all_sources) \
unit_test_43_driver.cc

PGFem3D_SOURCES = \
$(all_sources) \
super_main.c

VTK_IO_include = -I$(top_srcdir)/lib/VTK_IO/include
VTK_IO_lib = -L$(top_srcdir)/lib/VTK_IO/src -lPGFem3D_to_VTK

AM_INCLUDES = $(VTK_IO_include) \
	      $(hypre_include) \
	      $(suitesparse_include) \
	      $(cnstvm_include) \
	      $(mkl_include)

all_libs = $(VTK_IO_lib) $(suitesparse_lib) $(hypre_lib) $(vtk_lib) $(cnstvm_lib) $(mkl_lib)

PGFem3D_LDADD = $(all_libs)
BPA_driver_LDADD = $(all_libs)
IVD_driver_LDADD = $(all_libs)
J2D_driver_LDADD = $(all_libs)
CH_layer_driver_LDADD = $(all_libs)
unit_test_43_driver_LDADD = $(all_libs)

AM_CFLAGS = $(AM_INCLUDES)

PGFem3D_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
BPA_driver_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
IVD_driver_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
J2D_driver_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
CH_layer_driver_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
unit_test_43_driver_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
