SUBDIRS = solvers

# There are two targets in this source directory, the PGFem3D application
# itself, and a convenience library called libgpfem3d.la. The convenience
# library exists because elsewhere in the project there are executables that
# want to use the functionality, but want to write their own main.
bin_PROGRAMS = PGFem3D
noinst_LTLIBRARIES = libpgfem3d.la

PGFem3D_LDFLAGS = -Wl,--allow-multiple-definition # https://gitlab-cswarm.crc.nd.edu/Constitutive_Model/Generalizsed_constitutive_model/issues/9
PGFem3D_LDADD   = libpgfem3d.la
PGFem3D_SOURCES = super_main.cc

# Universal flags for all targets
AM_CPPFLAGS = -Wall -Werror -Wno-misleading-indentation \
              -I$(top_srcdir)/include \
              -I$(top_srcdir)/lib/VTK_IO/include \
              -I$(top_srcdir)/energy_equation/include \
              -I$(top_builddir)/include \
              $(VTK_IO_include) \
              $(energy_equation_include) \
              $(suitesparse_include) \
              $(ttl_include) \
              $(cnstvm_include) \
              $(mkl_include)

libpgfem3d_la_LIBADD = $(builddir)/solvers/libsolvers.la \
                       $(top_builddir)/lib/VTK_IO/src/libPGFem3D_to_VTK.la \
                       $(top_builddir)/energy_equation/src/libenergy_equation.la \
                       $(suitesparse_lib) \
                       $(VTK_LDFLAGS) $(VTK_LIBS) \
                       $(cnstvm_lib) \
                       $(mkl_lib)

libpgfem3d_la_SOURCES = ALM.cc \
                        Arc_length.cc \
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
