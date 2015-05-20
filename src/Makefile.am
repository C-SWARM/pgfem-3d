bin_PROGRAMS = PGFem3D
PGFem3D_SOURCES = \
ALM.c \
Arc_length.c \
boomerAMGInterface.c \
GRedist_node.c \
LINE.c \
Newton_Raphson.c \
PLC.c \
PLoc_Sparse.c \
Psparse_ApAi.c \
Re1_Re2_Re3.c \
TA_GA.c \
addresses.c \
build_distribution.c \
def_grad.c \
differentiation_of_UU.c \
elem3d.c \
fd_increment.c \
fd_residuals.c \
homogen.c \
hypre_interface.c \
in.c \
initial_guess.c \
integration.c \
load.c \
localizat.c \
macroscopic_load_AL.c \
matice.c \
matrix_printing.c \
metis_printing.c \
null.c \
out.c \
press_theta.c \
pressu_shape.c \
Printing.c \
print_dist.c \
read_cryst_plast.c \
renumber_ID.c \
res_fini_def.c \
resice.c \
resid_on_elem.c \
rn_skyline.c \
rowlength.c \
set_fini_def.c \
SetGlobalNodeNumbers.c \
skyline.c \
stabilized.c \
stiffmat_fd.c \
stiffmatel_fd.c \
stress_strain.c \
subdivision.c \
tensors.c \
utils.c \
incl.c \
profiler.c \
potential.c \
computeMacroF.c \
computeMacroS.c \
vtk_output.c \
MINI_element.c \
MINI_3f_element.c \
two_field_element.c \
Hu_Washizu_element.c \
new_potentials.c \
displacement_based_element.c \
PGFem3D_options.c \
gen_path.c \
volumetric_damage.c \
initialize_damage.c \
vol_damage_int_alg.c \
interface_macro.c \
compute_reactions.c \
bounding_element.c \
bounding_element_utils.c \
get_dof_ids_on_elem.c \
get_ndof_on_elem.c \
transform_coordinates.c \
generate_dof_ids.c \
node.c \
element.c \
solve_system.c \
PGFEM_HYPRE_preconditioners.c \
allocation.c \
cohesive_element_utils.c \
cohesive_element.c \
cohesive_potentials.c \
quadrature_rules.c \
integrate_surface.c \
applied_traction.c \
pgfem_comm.c \
tfa.c \
hommat.c \
sig.c \
supp.c \
mesh_load.c \
material.c \
eps.c \
matgeom.c \
crpl.c \
ensight.c \
read_input_file.c \
microscale_information.c \
ms_cohe_job_info.c \
ms_cohe_job_list.c \
PGFEM_par_matvec.c \
compute_ms_cohe_job.c \
PGFEM_io.c \
PGFEM_mpi.c \
ms_job_intercomm.c \
macro_micro_functions.c \
plasticity.c \
plasticity_model_none.c \
pgf_fe2_micro_server.c \
pgf_fe2_job.c \
pgf_fe2_server_rebalance.c \
pgf_fe2_rebalancer.c \
pgf_fe2_macro_client.c \
pgf_fe2_restart.c \
pgf_fe2_compute_max_n_jobs.c \
solver_file.c \
condense.c \
three_field_element.c \
femlib.c \
macro_micro_main.c \
main.c \
dynamics.c \
super_main.c

renumbering_include = -I$(top_srcdir)/lib/Renumbering/include
renumbering_lib = -L$(top_srcdir)/lib/Renumbering/src -lrenumbering

VTK_IO_include = -I$(top_srcdir)/lib/VTK_IO/include
VTK_IO_lib = -L$(top_srcdir)/lib/VTK_IO/src -lPGFem3D_to_VTK

AM_INCLUDES = $(renumbering_include) \
	      $(VTK_IO_include) \
	      $(hypre_include) \
	      $(suitesparse_include) \
	      $(mkl_include)

PGFem3D_LDADD = $(renumbering_lib) $(VTK_IO_lib) $(parmetis_lib) $(suitesparse_lib) $(hypre_lib) $(vtk_lib) $(mkl_lib)

AM_CFLAGS = $(AM_INCLUDES)

PGFem3D_LINK = $(PGFem3D_linker) $(AM_CFLAGS) $(CFLAGS) $(AM_LDFLAGS) $(LDFLAGS) -o $@
