SUBDIRS = 
noinst_LTLIBRARIES        = libmultiscale.la

libmultiscale_la_SOURCES  = PGFEM_par_matvec.cc            \
			    MultiscaleSolution.cc          \
			    MultiscaleCommon.cc            \
			    pgf_fe2_macro_client.cc        \
			    pgf_fe2_compute_max_n_jobs.cc  \
			    pgf_fe2_restart.cc             \
			    pgf_fe2_rebalancer.cc          \
			    pgf_fe2_server_rebalance.cc    \
			    pgf_fe2_micro_server.cc        \
			    pgf_fe2_job.cc                 \
			    ms_cohe_job_info.cc            \
			    ms_cohe_job_list.cc            \
			    compute_ms_cohe_job.cc         \
			    macro_micro_functions.cc       \
			    macro_micro_main.cc

libmultiscale_la_CPPFLAGS = -Wall -Werror                  \
                            -I$(top_srcdir)/include        \
                            -I$(top_builddir)/include      \
	                    -I$(top_srcdir)/FEMLIB/include \
                             $(LIBGCM_CPPFLAGS)
